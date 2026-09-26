use cosmolkit_model::{
    Atom, AtomId, AtomRangeBounds, AtomRangeDataFunction, AtomSpec, Bond, BondId, BondSpec,
    Element, QueryAtom, QueryAtomIdentity, TopologyBlock,
};
use cosmolkit_search::{
    AtomQueryPredicate, BondQueryPredicate, QueryGraph, QueryNode, SmartsParseError,
    SmartsParseParams, SmartsWriteParams, match_query, parse_smarts, write_smarts,
};
use cosmolkit_types::{BondDirection, BondOrder, BondStereo, ChiralTag, Hybridization};

fn parse_source_case(smarts: &str) -> QueryGraph {
    parse_smarts(smarts, &SmartsParseParams::default())
        .unwrap_or_else(|error| panic!("pinned RDKit SMARTS {smarts:?}: {error}"))
}

fn assert_atom(
    atom: &QueryAtom,
    input: &str,
    identity: QueryAtomIdentity,
    predicate: QueryNode<AtomQueryPredicate>,
) {
    assert_eq!(atom.identity(), identity, "{input}");
    assert_eq!(atom.atomic_number(), identity.atomic_number(), "{input}");
    assert!(atom.is_aromatic(), "{input}");
    assert_eq!(atom.predicate(), &predicate, "{input}");
    assert!(!atom.predicate_is_carrier_derived(), "{input}");
}

#[test]
fn generic_a_and_as_follow_bracket_lexical_state_and_carrier_identity() {
    // Pinned smarts.ll emits generic `a` as SIMPLE_ATOM_QUERY_TOKEN in both
    // states, but its two-character `as` token exists only in atom state.
    let generic_a = QueryNode::predicate(AtomQueryPredicate::IsAromatic(true));
    for smarts in ["a", "[a]"] {
        let graph = parse_source_case(smarts);
        assert_eq!(graph.num_atoms(), 1, "{smarts}");
        assert_atom(
            graph.atom(0).expect("generic aromatic query atom"),
            smarts,
            QueryAtomIdentity::Element(Element::DUMMY),
            generic_a.clone(),
        );
    }

    let arsenic_query = QueryNode::predicate(AtomQueryPredicate::AtomType {
        atomic_number: 33,
        aromatic: true,
    });
    let bracketed_as = parse_source_case("[as]");
    assert_eq!(bracketed_as.num_atoms(), 1);
    assert_atom(
        bracketed_as.atom(0).expect("single aromatic arsenic atom"),
        "[as]",
        QueryAtomIdentity::Element(Element::AS),
        arsenic_query,
    );

    // Outside IN_ATOM_STATE, the source matches `a` then `s`, not `[as]`.
    let unbracketed_as = parse_source_case("as");
    assert_eq!(unbracketed_as.num_atoms(), 2);
    assert_atom(
        unbracketed_as.atom(0).expect("generic aromatic prefix"),
        "as",
        QueryAtomIdentity::Element(Element::DUMMY),
        generic_a,
    );
    assert_atom(
        unbracketed_as.atom(1).expect("aromatic sulfur suffix"),
        "as",
        QueryAtomIdentity::Element(Element::S),
        QueryNode::predicate(AtomQueryPredicate::AtomType {
            atomic_number: 16,
            aromatic: true,
        }),
    );
}

#[test]
fn q07a_element_numeric_and_wildcard_atoms_preserve_identity_and_origin() {
    // Pinned smarts.yy sends organic symbols through simple_atom/AtomType,
    // ATOM_TOKEN and #N through numeric identity, and `*` through AtomNull.
    // The scanner's legacy Uut/Uup aliases resolve to atomic numbers 113/115.
    let cases = [
        (
            "[C]",
            6,
            QueryNode::predicate(AtomQueryPredicate::AtomType {
                atomic_number: 6,
                aromatic: false,
            }),
        ),
        (
            "[Cl]",
            17,
            QueryNode::predicate(AtomQueryPredicate::AtomType {
                atomic_number: 17,
                aromatic: false,
            }),
        ),
        (
            "[H]",
            1,
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(1)),
        ),
        (
            "[He]",
            2,
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(2)),
        ),
        (
            "[Hf]",
            72,
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(72)),
        ),
        (
            "[Hg]",
            80,
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(80)),
        ),
        (
            "[Ho]",
            67,
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(67)),
        ),
        (
            "[Hs]",
            108,
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(108)),
        ),
        (
            "[Uut]",
            113,
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(113)),
        ),
        (
            "[Uup]",
            115,
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(115)),
        ),
        (
            "[#0]",
            0,
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(0)),
        ),
        (
            "[#119]",
            119,
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(119)),
        ),
        ("*", 0, QueryNode::predicate(AtomQueryPredicate::Any)),
        ("[*]", 0, QueryNode::predicate(AtomQueryPredicate::Any)),
    ];

    for (smarts, atomic_number, predicate) in cases {
        let graph = parse_source_case(smarts);
        assert_eq!(graph.num_atoms(), 1, "{smarts}");
        let atom = graph.atom(0).expect("one source atom");
        assert_eq!(
            atom.identity(),
            QueryAtomIdentity::from_atomic_number(atomic_number),
            "{smarts}"
        );
        assert_eq!(atom.atomic_number(), atomic_number, "{smarts}");
        assert_eq!(atom.predicate(), &predicate, "{smarts}");
        assert!(
            !atom.predicate_is_carrier_derived(),
            "SMARTS query origin is explicit: {smarts}"
        );
    }

    let numeric_zero = parse_source_case("[#0]");
    let wildcard = parse_source_case("*");
    assert_eq!(
        numeric_zero.atom(0).unwrap().identity(),
        wildcard.atom(0).unwrap().identity()
    );
    assert_ne!(
        numeric_zero.atom(0).unwrap().predicate(),
        wildcard.atom(0).unwrap().predicate()
    );
}

#[test]
fn scanner_setup_trims_control_bytes_only_at_input_edges() {
    let empty = parse_source_case("");
    assert_eq!(empty.num_atoms(), 0);

    // toMol accepts a truly empty input, while nonempty edge controls reach
    // the parser after setup_smarts_string trims them to its terminal NUL.
    assert!(parse_smarts(" \t\r\n", &SmartsParseParams::default()).is_err());

    for control in [" ", "\t", "\r", "\n"] {
        let smarts = format!("{control}C{control}");
        let graph = parse_source_case(&smarts);
        assert_eq!(graph.num_atoms(), 1, "{smarts:?}");
    }

    let leading_newline = parse_source_case("\nC");
    assert_eq!(leading_newline.num_atoms(), 1);

    let raw_smarts = SmartsParseParams {
        allow_cxsmiles: false,
        parse_name: false,
        ..SmartsParseParams::default()
    };
    for control in [" ", "\t", "\r"] {
        let smarts = format!("C{control}C");
        assert!(
            parse_smarts(&smarts, &raw_smarts).is_err(),
            "internal {control:?} must remain a BAD_CHARACTER"
        );
    }

    // An internal LF remains the source EOS token and terminates after C.
    let internal_newline = parse_smarts("C\nC", &raw_smarts).expect("EOS accepts prefix");
    assert_eq!(internal_newline.num_atoms(), 1);
}

#[test]
fn q03_parser_errors_use_trimmed_parser_byte_positions() {
    let error = parse_smarts("\t  CCC☃C \t", &SmartsParseParams::default())
        .expect_err("snowman is a BAD_CHARACTER after the valid prefix");
    assert_eq!(
        error,
        SmartsParseError::UnexpectedCharacter {
            position: 4,
            character: '☃',
            context: "unexpected character in SMARTS string".to_owned(),
        }
    );

    let eof_error = parse_smarts("C(", &SmartsParseParams::default())
        .expect_err("a branch must begin with atomd or bond_expr atomd before EOS");
    assert_eq!(
        eof_error,
        SmartsParseError::UnexpectedCharacter {
            position: 2,
            character: '?',
            context: "expected atom expression".to_owned(),
        }
    );

    let missing_close = parse_smarts("C(O", &SmartsParseParams::default())
        .expect_err("a populated branch reaches EOS without GROUP_CLOSE_TOKEN");
    assert_eq!(
        missing_close,
        SmartsParseError::UnexpectedCharacter {
            position: 3,
            character: 'E',
            context: "expected close parenthesis".to_owned(),
        }
    );
}

#[test]
fn q13_branch_stack_restores_attachment_and_rejects_invalid_branch_starts() {
    // The pinned mol productions attach each branch's first atom to the
    // active atom, then GROUP_CLOSE restores the branch attachment point.
    for (smarts, expected_endpoints) in [
        ("C(O)(N)F", vec![(0, 1), (0, 2), (0, 3)]),
        ("C(O(N))F", vec![(0, 1), (1, 2), (0, 3)]),
    ] {
        let graph = parse_source_case(smarts);
        assert_eq!((graph.num_atoms(), graph.num_bonds()), (4, 3), "{smarts}");
        let endpoints = (0..graph.num_bonds())
            .map(|index| graph.bond(index).expect("branch bond").endpoints())
            .collect::<Vec<_>>();
        assert_eq!(endpoints, expected_endpoints, "{smarts}");
    }

    let explicit_branch_bond = parse_source_case("C(=O)N");
    assert_eq!(
        (
            explicit_branch_bond.num_atoms(),
            explicit_branch_bond.num_bonds()
        ),
        (3, 2)
    );
    assert_eq!(
        explicit_branch_bond
            .bond(0)
            .expect("explicit branch bond")
            .endpoints(),
        (0, 1)
    );
    assert_eq!(
        explicit_branch_bond.bond(0).unwrap().bond().order(),
        BondOrder::Double
    );
    assert_eq!(
        explicit_branch_bond
            .bond(1)
            .expect("post-branch bond")
            .endpoints(),
        (0, 2)
    );

    // smarts.yy requires atomd or bond_expr atomd immediately after the
    // opening parenthesis. YY_USER_ACTION reports the ending byte of the
    // offending token, including EOS for the missing first branch atom.
    for (smarts, position) in [
        ("C(", 2),
        ("C()", 3),
        ("C(.O)", 3),
        ("C(1)1", 3),
        ("C((O))", 3),
    ] {
        assert_eq!(
            parse_smarts(smarts, &SmartsParseParams::default())
                .expect_err("invalid branch start must fail at its source token"),
            SmartsParseError::UnexpectedCharacter {
                position,
                character: '?',
                context: "expected atom expression".to_owned(),
            },
            "{smarts}"
        );
    }
}

#[test]
fn q14_disconnected_components_keep_source_order_and_reject_dangling_tokens() {
    // The pinned separator production appends a new source atom without a
    // bond, then makes it active for subsequent bond/atom productions.
    let disconnected = parse_source_case("O.C.N");
    assert_eq!((disconnected.num_atoms(), disconnected.num_bonds()), (3, 0));
    let source_order = (0..disconnected.num_atoms())
        .map(|index| {
            disconnected
                .atom(index)
                .expect("source ordered component atom")
                .atomic_number()
        })
        .collect::<Vec<_>>();
    assert_eq!(source_order, [8, 6, 7]);

    let continued_component = parse_source_case("C.N-O");
    assert_eq!(
        (
            continued_component.num_atoms(),
            continued_component.num_bonds()
        ),
        (3, 1)
    );
    assert_eq!(
        continued_component
            .bond(0)
            .expect("bond remains in the component after the separator")
            .endpoints(),
        (1, 2)
    );
    assert_eq!(
        continued_component.bond(0).unwrap().predicate(),
        &QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single))
    );

    // toMol returns an empty value only for an exactly empty input. The
    // START_MOL grammar requires an atom for nonempty scanner input.
    let params = SmartsParseParams {
        allow_cxsmiles: false,
        parse_name: false,
        ..SmartsParseParams::default()
    };
    let empty = parse_smarts("", &params).expect("exact empty input maps to an empty graph");
    assert_eq!((empty.num_atoms(), empty.num_bonds()), (0, 0));
    assert_eq!(
        parse_smarts(" \t\r\n", &params).expect_err("nonempty whitespace is parsed"),
        SmartsParseError::UnexpectedEnd("expected atom but reached end".to_owned())
    );

    // A separator or explicit bond cannot finish without its required atomd.
    for smarts in ["C.", "C-", "C.N-"] {
        assert_eq!(
            parse_smarts(smarts, &params).expect_err("dangling source token must fail"),
            SmartsParseError::UnexpectedEnd("expected atom but reached end".to_owned()),
            "{smarts}"
        );
    }
}

#[test]
fn q15_ring_closures_follow_bookmark_order_priority_and_source_errors() {
    // Pinned CloseMolRings walks sorted atom-bookmark labels, pairs each
    // label's atoms in insertion order, and gives the first bond without
    // `_unspecifiedOrder` priority. Query identity remains separate from the
    // final carrier type projected by SetUnspecifiedBondTypes.
    let single_or_aromatic = QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
        BondOrder::Single,
        BondOrder::Aromatic,
    ]));
    for (smarts, endpoints, order, predicate) in [
        (
            "C=1CC#1",
            (0, 2),
            BondOrder::Double,
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Double)),
        ),
        (
            "C1CC#1",
            (2, 0),
            BondOrder::Triple,
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Triple)),
        ),
        (
            "C~1CC1",
            (0, 2),
            BondOrder::Unspecified,
            QueryNode::predicate(BondQueryPredicate::Any),
        ),
        (
            "C1CC~1",
            (2, 0),
            BondOrder::Unspecified,
            QueryNode::predicate(BondQueryPredicate::Any),
        ),
    ] {
        let graph = parse_source_case(smarts);
        assert_eq!((graph.num_atoms(), graph.num_bonds()), (3, 3), "{smarts}");
        let closure = graph.bond(2).expect("ring closure follows chain bonds");
        assert_eq!(closure.endpoints(), endpoints, "{smarts}");
        assert_eq!(closure.bond().order(), order, "{smarts}");
        assert_eq!(closure.predicate(), &predicate, "{smarts}");
        assert!(!closure.predicate_is_carrier_derived(), "{smarts}");
    }

    let aliphatic_default = parse_source_case("C1CC1");
    let aliphatic_closure = aliphatic_default.bond(2).expect("aliphatic closure");
    assert_eq!(aliphatic_closure.endpoints(), (2, 0));
    assert_eq!(aliphatic_closure.bond().order(), BondOrder::Single);
    assert_eq!(aliphatic_closure.predicate(), &single_or_aromatic);

    let aromatic_default = parse_source_case("c1cc1");
    let aromatic_closure = aromatic_default.bond(2).expect("aromatic closure");
    assert_eq!(aromatic_closure.endpoints(), (2, 0));
    assert_eq!(aromatic_closure.bond().order(), BondOrder::Aromatic);
    assert_eq!(aromatic_closure.predicate(), &single_or_aromatic);

    // Source closure rows are appended after ordinary bonds by sorted label,
    // even when the final token closes label 2 before label 1.
    let sorted_labels = parse_source_case("C1CC2CC21");
    let endpoints = sorted_labels
        .bonds()
        .iter()
        .map(|bond| bond.endpoints())
        .collect::<Vec<_>>();
    assert_eq!(endpoints, [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0), (4, 2)]);

    // The bookmark list pairs repeated uses of one label from left to right.
    let repeated_label = parse_source_case("C1CC1C1CC1");
    let closure_endpoints = repeated_label.bonds()[5..]
        .iter()
        .map(|bond| bond.endpoints())
        .collect::<Vec<_>>();
    assert_eq!(closure_endpoints, [(2, 0), (5, 3)]);

    for (smarts, expected) in [
        (
            "C11",
            SmartsParseError::Parse("duplicated ring closure 1 bonds atom 0 to itself".to_owned()),
        ),
        (
            "C1C1",
            SmartsParseError::Parse(
                "ring closure 1 duplicates bond between atom 0 and atom 1".to_owned(),
            ),
        ),
        ("C1CC", SmartsParseError::Parse("unclosed ring".to_owned())),
    ] {
        assert_eq!(
            parse_smarts(smarts, &SmartsParseParams::default())
                .expect_err("pinned ring closure failure"),
            expected,
            "{smarts}"
        );
    }
}

#[test]
fn q16_dative_endpoint_orientation_and_source_row_order() {
    // The pinned direct and branch bond_expr atomd actions normalize DATIVER
    // and DATIVEL endpoints, then append the bond at that source reduction.
    let right = QueryNode::predicate(BondQueryPredicate::Order(BondOrder::DativeRight));
    let left = QueryNode::predicate(BondQueryPredicate::Order(BondOrder::DativeLeft));
    let single = QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single));
    let implicit_single_or_aromatic = QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
        BondOrder::Single,
        BondOrder::Aromatic,
    ]));

    let direct = parse_source_case("C->N-C<-O");
    assert_eq!((direct.num_atoms(), direct.num_bonds()), (4, 3));
    for (index, endpoints, order, predicate) in [
        (0, (0, 1), BondOrder::Dative, &right),
        (1, (1, 2), BondOrder::Single, &single),
        (2, (3, 2), BondOrder::Dative, &left),
    ] {
        let bond = direct.bond(index).expect("source-ordered direct bond");
        assert_eq!(bond.endpoints(), endpoints, "direct bond {index}");
        assert_eq!(bond.bond().order(), order, "direct bond {index}");
        assert_eq!(bond.predicate(), predicate, "direct bond {index}");
        assert!(!bond.predicate_is_carrier_derived(), "direct bond {index}");
    }

    for (smarts, endpoints, predicate) in [("C(->N)O", (0, 1), &right), ("C(<-N)O", (1, 0), &left)]
    {
        let graph = parse_source_case(smarts);
        assert_eq!((graph.num_atoms(), graph.num_bonds()), (3, 2), "{smarts}");
        let branch_bond = graph
            .bond(0)
            .expect("branch dative row precedes continuation");
        assert_eq!(branch_bond.endpoints(), endpoints, "{smarts}");
        assert_eq!(branch_bond.bond().order(), BondOrder::Dative, "{smarts}");
        assert_eq!(branch_bond.predicate(), predicate, "{smarts}");
        assert!(!branch_bond.predicate_is_carrier_derived(), "{smarts}");

        let continuation = graph.bond(1).expect("post-branch continuation row");
        assert_eq!(continuation.endpoints(), (0, 2), "{smarts}");
        assert_eq!(continuation.bond().order(), BondOrder::Single, "{smarts}");
        assert_eq!(
            continuation.predicate(),
            &implicit_single_or_aromatic,
            "{smarts}"
        );
    }
}

#[test]
fn q17_recursive_labels_follow_source_order_and_preserve_explicit_duplicates() {
    fn collect_recursive_serials(graph: &QueryGraph, descend: bool, serials: &mut Vec<u32>) {
        fn visit(node: &QueryNode<AtomQueryPredicate>, descend: bool, serials: &mut Vec<u32>) {
            match node {
                QueryNode::Predicate(AtomQueryPredicate::RecursiveSmarts(recursive)) => {
                    serials.push(recursive.serial_number());
                    if descend {
                        if let Some(nested) = recursive.query_graph() {
                            collect_recursive_serials(nested, true, serials);
                        }
                    }
                }
                QueryNode::Predicate(_) => {}
                QueryNode::And(children) | QueryNode::Or(children) | QueryNode::Xor(children) => {
                    for child in children {
                        visit(child, descend, serials);
                    }
                }
                QueryNode::Not(child) => visit(child, descend, serials),
            }
        }

        for atom_index in 0..graph.num_atoms() {
            let atom = graph.atom(atom_index).expect("query atom index");
            visit(atom.predicate(), descend, serials);
        }
    }

    // Pinned labelRecursivePatterns allocates IDs as recursive expressions
    // close, and uses exact pattern text to reuse an earlier automatic label.
    let source_order = parse_source_case("[$(N)]-[$(C)]-[$(N)]");
    let mut serials = Vec::new();
    collect_recursive_serials(&source_order, false, &mut serials);
    assert_eq!(serials, [100, 101, 100]);

    // Nested closures receive their label before the containing expression;
    // each repeated outer and inner pattern retains its own source label.
    let nested = parse_source_case("[$(C[$(N)])]-[$(C[$(N)])]");
    serials.clear();
    collect_recursive_serials(&nested, true, &mut serials);
    assert_eq!(serials, [101, 100, 101, 100]);

    // Explicit repeated labels bypass automatic registration, so the later
    // automatic N pattern still receives the first generated serial.
    let explicit = parse_source_case("[$(N)_7]-[$(C)_7]-[$(N)]");
    serials.clear();
    collect_recursive_serials(&explicit, false, &mut serials);
    assert_eq!(serials, [7, 7, 100]);
}

#[test]
fn q18_recursive_graph_ownership_root_errors_and_source_serial_limit() {
    fn first_recursive_query(node: &QueryNode<AtomQueryPredicate>) -> Option<(u32, &QueryGraph)> {
        match node {
            QueryNode::Predicate(AtomQueryPredicate::RecursiveSmarts(recursive)) => recursive
                .query_graph()
                .map(|graph| (recursive.serial_number(), graph)),
            QueryNode::Predicate(_) => None,
            QueryNode::And(children) | QueryNode::Or(children) | QueryNode::Xor(children) => {
                children.iter().find_map(first_recursive_query)
            }
            QueryNode::Not(child) => first_recursive_query(child),
        }
    }

    // The recursive predicate owns its compiled inner graph. Source atom
    // insertion order supplies the graph root at index zero at each nesting.
    let graph = parse_source_case("[$([O;$([N])]-S)]");
    let (outer_serial, outer_query) =
        first_recursive_query(graph.atom(0).expect("recursive host atom").predicate())
            .expect("owned outer recursive query graph");
    assert_eq!(outer_serial, 101);
    assert_eq!(outer_query.num_atoms(), 2);
    assert_eq!(
        outer_query.atom(0).expect("source root atom").identity(),
        QueryAtomIdentity::Element(Element::O)
    );
    assert_eq!(
        outer_query.atom(1).expect("source second atom").identity(),
        QueryAtomIdentity::Element(Element::S)
    );
    assert_eq!(
        outer_query.bond(0).expect("source bond").endpoints(),
        (0, 1)
    );

    let (inner_serial, inner_query) =
        first_recursive_query(outer_query.atom(0).expect("inner root atom").predicate())
            .expect("owned nested recursive query graph");
    assert_eq!(inner_serial, 100);
    assert_eq!(inner_query.num_atoms(), 1);
    assert_eq!(
        inner_query.atom(0).expect("nested source root").identity(),
        QueryAtomIdentity::Element(Element::N)
    );

    // Empty recursive molecules, missing close delimiters, and malformed
    // inner bond syntax remain parser errors instead of empty query graphs.
    for invalid in ["[$()]", "[$(O]", "[$(O-)]"] {
        assert!(
            parse_smarts(invalid, &SmartsParseParams::default()).is_err(),
            "{invalid}"
        );
    }

    let source_max = parse_source_case("[$(N)_2147483639]");
    let (max_serial, _) = first_recursive_query(
        source_max
            .atom(0)
            .expect("max serial host atom")
            .predicate(),
    )
    .expect("max serial recursive graph");
    assert_eq!(max_serial, 2_147_483_639);
    assert_eq!(
        parse_smarts("[$(N)_2147483640]", &SmartsParseParams::default())
            .expect_err("source nonzero_number limit"),
        SmartsParseError::InvalidAtomPrimitive {
            position: 5,
            detail: "number too large".to_string(),
        }
    );
}

#[test]
fn q19_preprocessing_defaults_replacements_name_delimiter_and_malformed_input() {
    // Defaults mirror the pinned SmartsParserParams declaration.
    let defaults = SmartsParseParams::default();
    assert!(defaults.allow_cxsmiles);
    assert!(defaults.strict_cxsmiles);
    assert!(defaults.parse_name);
    assert!(!defaults.merge_hs);
    assert!(!defaults.skip_cleanup);
    assert!(!defaults.debug_parse);
    assert!(defaults.replacements.is_empty());

    // preprocessSmiles splits only at the first ASCII space/tab after byte 0,
    // applies ordered global replacements to the SMARTS prefix only, and
    // repeats the ordered map until a later-key replacement can reach an
    // earlier key. Boost trim_copy retains trailing non-ASCII NBSP in the name.
    let mut params = SmartsParseParams {
        allow_cxsmiles: false,
        ..SmartsParseParams::default()
    };
    params.replacements.insert("N".to_string(), "C".to_string());
    params.replacements.insert("O".to_string(), "N".to_string());
    let graph = parse_smarts("O-O\tOriginal O\u{00a0}", &params)
        .expect("source replacement closure and name split");
    assert_eq!((graph.num_atoms(), graph.num_bonds()), (2, 1));
    for index in 0..2 {
        assert_eq!(
            graph.atom(index).expect("replacement atom").identity(),
            QueryAtomIdentity::Element(Element::C)
        );
    }
    assert_eq!(graph.name(), Some("Original O\u{00a0}"));

    // A delimiter at offset zero is not a split point; its remaining embedded
    // whitespace stays malformed SMARTS, and replacement output is not
    // recovered through a guessed fallback after the parser rejects it.
    assert!(parse_smarts("\tC name", &params).is_err());
    let mut malformed_replacement = params;
    malformed_replacement.replacements.clear();
    malformed_replacement
        .replacements
        .insert("C".to_string(), "?".to_string());
    assert!(parse_smarts("C", &malformed_replacement).is_err());
}

#[test]
fn q20_wrapper_cx_name_flags_and_error_order() {
    let valid_cx = parse_source_case("C |$label$| note");
    assert_eq!(
        valid_cx.atom(0).and_then(|atom| atom.prop("atomLabel")),
        Some("label")
    );
    assert_eq!(valid_cx.prop("_CXSMILES_Data"), Some("|$label$|"));
    assert_eq!(valid_cx.name(), Some("note"));

    let no_name_params = SmartsParseParams {
        parse_name: false,
        ..SmartsParseParams::default()
    };
    let no_name = parse_smarts("C |$label$| note", &no_name_params)
        .expect("parseName=false still applies valid CX records");
    assert_eq!(
        no_name.atom(0).and_then(|atom| atom.prop("atomLabel")),
        Some("label")
    );
    assert_eq!(no_name.prop("_CXSMILES_Data"), Some("|$label$|"));
    assert_eq!(no_name.name(), None);

    let no_cx_params = SmartsParseParams {
        allow_cxsmiles: false,
        ..SmartsParseParams::default()
    };
    let cx_text_is_name = parse_smarts("C |$label$| note", &no_cx_params)
        .expect("allowCXSMILES=false treats the full suffix as the name");
    assert_eq!(cx_text_is_name.name(), Some("|$label$| note"));
    assert_eq!(cx_text_is_name.prop("_CXSMILES_Data"), None);

    let retained_non_ascii =
        parse_smarts("C |$label$| note\u{00a0}", &SmartsParseParams::default())
            .expect("source-trimmed CX suffix");
    assert_eq!(retained_non_ascii.name(), Some("note\u{00a0}"));

    let strict_no_name = SmartsParseParams {
        parse_name: false,
        ..SmartsParseParams::default()
    };
    assert!(matches!(
        parse_smarts("C note", &strict_no_name),
        Err(SmartsParseError::CxSmiles(_))
    ));
    let lenient_no_name = SmartsParseParams {
        strict_cxsmiles: false,
        ..strict_no_name.clone()
    };
    let parsed_without_name = parse_smarts("C note", &lenient_no_name)
        .expect("strictCXSMILES=false accepts the source name suffix");
    assert_eq!(parsed_without_name.name(), None);

    let strict_malformed_cx = parse_smarts("C |sense|", &SmartsParseParams::default());
    assert!(matches!(
        strict_malformed_cx,
        Err(SmartsParseError::CxSmiles(_))
    ));

    // MolFromSmarts calls toMol before handleCXPartAndName: a malformed SMARTS
    // prefix therefore wins over a malformed CX suffix.
    let parser_error_first = parse_smarts("C- |sense|", &SmartsParseParams::default());
    assert!(!matches!(
        parser_error_first,
        Err(SmartsParseError::CxSmiles(_))
    ));
}

#[test]
fn q20_lenient_cx_failure_retains_source_cursor_property() {
    let params = SmartsParseParams {
        strict_cxsmiles: false,
        ..SmartsParseParams::default()
    };
    let graph = parse_smarts("C |sense| ignored", &params)
        .expect("lenient CX failure preserves the successfully parsed SMARTS");
    assert_eq!(graph.name(), None);
    // Pinned handleCXPartAndName stores [cxPart.begin(), parser iterator).
    // parse_substitution rejects `sense` at its initial `s`, so the property is `|`.
    assert_eq!(graph.prop("_CXSMILES_Data"), Some("|"));
}

#[test]
fn q20_lenient_cx_failure_retains_prior_record_effects() {
    let params = SmartsParseParams {
        strict_cxsmiles: false,
        ..SmartsParseParams::default()
    };
    let graph = parse_smarts("C |$label$ s:0:x| ignored", &params)
        .expect("lenient CX failure preserves the successfully parsed SMARTS");
    // The pinned parser installs the label before parse_substitution fails;
    // its iterator remains at `x`, so the CX data prefix also includes `s:0:`.
    assert_eq!(
        graph.atom(0).and_then(|atom| atom.prop("atomLabel")),
        Some("label")
    );
    assert_eq!(graph.prop("_CXSMILES_Data"), Some("|$label$ s:0:"));
    assert_eq!(graph.name(), None);
}

#[test]
fn q20_lenient_cx_lowering_failure_retains_record_effects_and_cursor() {
    let params = SmartsParseParams {
        strict_cxsmiles: false,
        ..SmartsParseParams::default()
    };
    // The second source pair fails the duplicate-wedge check after its bond
    // index is read, so the prefix excludes the closing pipe and later name.
    let graph = parse_smarts("C-C |$label$ wU:0.0,1.0| ignored", &params)
        .expect("lenient lowering failure preserves earlier source mutations");
    assert_eq!(
        graph.atom(0).and_then(|atom| atom.prop("atomLabel")),
        Some("label")
    );
    let bond = graph.bond(0).expect("first wedge bond remains");
    assert_eq!(bond.endpoints(), (0, 1));
    assert_eq!(bond.bond().direction(), BondDirection::BeginWedge);
    assert_eq!(bond.bond().prop("_MolFileBondCfg"), Some("1"));
    assert_eq!(graph.prop("_CXSMILES_Data"), Some("|$label$ wU:0.0,1.0"));
    assert_eq!(graph.name(), None);
}

#[test]
fn q20_strict_cx_lowering_failure_returns_cx_error() {
    assert!(matches!(
        parse_smarts("C-C |wU:0.0,1.0| ignored", &SmartsParseParams::default()),
        Err(SmartsParseError::CxSmiles(_))
    ));
}

#[test]
fn q03_bad_character_dispatch_preserves_byte_position_and_parser_priority() {
    let raw_smarts = SmartsParseParams {
        allow_cxsmiles: false,
        parse_name: false,
        ..SmartsParseParams::default()
    };

    for (input, character, position) in [
        ("C C", ' ', 2),
        ("C\tC", '\t', 2),
        ("C\rC", '\r', 2),
        ("C?C", '?', 2),
        ("CC☃C", '☃', 3),
        ("[C?]", '?', 3),
        ("C(C?C)", '?', 4),
    ] {
        assert_eq!(
            parse_smarts(input, &raw_smarts).expect_err("BAD_CHARACTER dispatch"),
            SmartsParseError::UnexpectedCharacter {
                position,
                character,
                context: "unexpected character in SMARTS string".to_owned(),
            },
            "{input:?}"
        );
    }

    // The pinned parser reduces mol GROUP_CLOSE_TOKEN and aborts at `)`
    // before asking the scanner for the following snowman token.
    assert_eq!(
        parse_smarts("C)☃", &raw_smarts)
            .expect_err("extra close parenthesis precedes later BAD_CHARACTER"),
        SmartsParseError::UnexpectedCharacter {
            position: 2,
            character: ')',
            context: "unexpected trailing token in molecule SMARTS".to_owned(),
        }
    );

    let newline_prefix =
        parse_smarts("C\nC", &raw_smarts).expect("internal LF is the scanner's consumed EOS token");
    assert_eq!(newline_prefix.num_atoms(), 1);
}

#[test]
fn q03_common_bond_tokens_keep_source_order_and_arrow_endpoints() {
    let cases: [(&str, BondOrder, (usize, usize)); 11] = [
        ("C-C", BondOrder::Single, (0, 1)),
        ("C=C", BondOrder::Double, (0, 1)),
        ("C#N", BondOrder::Triple, (0, 1)),
        ("C:N", BondOrder::Aromatic, (0, 1)),
        ("C~N", BondOrder::Unspecified, (0, 1)),
        ("C$N", BondOrder::Quadruple, (0, 1)),
        ("C/N", BondOrder::Single, (0, 1)),
        ("C\\N", BondOrder::Single, (0, 1)),
        ("C\\\\N", BondOrder::Single, (0, 1)),
        ("C->N", BondOrder::Dative, (0, 1)),
        ("C<-N", BondOrder::Dative, (1, 0)),
    ];

    for (smarts, expected_order, expected_endpoints) in cases {
        let graph = parse_source_case(smarts);
        assert_eq!((graph.num_atoms(), graph.num_bonds()), (2, 1), "{smarts}");
        let bond = graph.bond(0).expect("one source bond token");
        assert_eq!(bond.bond().order(), expected_order, "{smarts}");
        assert_eq!(bond.endpoints(), expected_endpoints, "{smarts}");
    }

    // Flex's [\\]{1,2} rule emits one BOND_TOKEN for either spelling.
    // The one-token query carrier and direction are therefore identical.
    let one_backslash = parse_source_case(r"C\N");
    let two_backslashes = parse_source_case(r"C\\N");
    assert_eq!(one_backslash.bond(0), two_backslashes.bond(0));
}

#[test]
fn q11_bond_boolean_precedence_preserves_left_carrier_and_direction() {
    // Pinned smarts.yy orders SEMI < OR < AND and keeps the left QueryBond at
    // each expansion. The parser's adjacent bond_query reduction is also AND.
    let graph = parse_source_case("C-,:&=;@N");
    assert_eq!((graph.num_atoms(), graph.num_bonds()), (2, 1));
    let bond = graph.bond(0).expect("one Boolean bond query");
    assert_eq!(
        bond.predicate(),
        &QueryNode::and(vec![
            QueryNode::or(vec![
                QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
                QueryNode::and(vec![
                    QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Aromatic)),
                    QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Double)),
                ]),
            ]),
            QueryNode::predicate(BondQueryPredicate::IsInRing(true)),
        ])
    );
    assert_eq!(bond.bond().order(), BondOrder::Single);
    assert_eq!(bond.bond().direction(), BondDirection::None);
    assert_eq!(bond.endpoints(), (0, 1));
    assert!(!bond.predicate_is_carrier_derived());

    // The slash lexer action supplies a single-or-aromatic predicate plus a
    // direction carrier. OR expansion keeps that first carrier when the
    // right-hand aromatic primitive is added.
    let directional = parse_source_case("C/,:N");
    let directional_bond = directional.bond(0).expect("one directional OR bond");
    assert_eq!(
        directional_bond.predicate(),
        &QueryNode::or(vec![
            QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
                BondOrder::Single,
                BondOrder::Aromatic,
            ])),
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Aromatic)),
        ])
    );
    assert_eq!(directional_bond.bond().order(), BondOrder::Single);
    assert_eq!(
        directional_bond.bond().direction(),
        BondDirection::EndUpRight
    );
    assert_eq!(directional_bond.endpoints(), (0, 1));

    for (smarts, expected) in [
        (
            "C!@N",
            QueryNode::not(QueryNode::predicate(BondQueryPredicate::IsInRing(true))),
        ),
        (
            "C!!@N",
            QueryNode::predicate(BondQueryPredicate::IsInRing(true)),
        ),
    ] {
        let graph = parse_source_case(smarts);
        let bond = graph.bond(0).expect("one negated bond query");
        assert_eq!(bond.predicate(), &expected, "{smarts}");
        assert_eq!(bond.bond().order(), BondOrder::Unspecified, "{smarts}");
        assert_eq!(bond.bond().direction(), BondDirection::None, "{smarts}");
        assert!(!bond.predicate_is_carrier_derived(), "{smarts}");
    }
}

#[test]
fn q12_bond_primitive_predicates_keep_source_carrier_defaults() {
    // smarts.yy constructs the single/triple/aromatic/ring primitives;
    // smarts.ll supplies the remaining BOND_TOKEN query and carrier actions.
    // Dative parsing changes the carrier type and endpoints after construction,
    // while the QueryBond predicate retains its source directional target.
    let cases = [
        (
            "C-C",
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
            BondOrder::Single,
            BondDirection::None,
            (0, 1),
        ),
        (
            "C=C",
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Double)),
            BondOrder::Double,
            BondDirection::None,
            (0, 1),
        ),
        (
            "C#N",
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Triple)),
            BondOrder::Triple,
            BondDirection::None,
            (0, 1),
        ),
        (
            "C:C",
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Aromatic)),
            BondOrder::Aromatic,
            BondDirection::None,
            (0, 1),
        ),
        (
            "C~N",
            QueryNode::predicate(BondQueryPredicate::Any),
            BondOrder::Unspecified,
            BondDirection::None,
            (0, 1),
        ),
        (
            "C@N",
            QueryNode::predicate(BondQueryPredicate::IsInRing(true)),
            BondOrder::Unspecified,
            BondDirection::None,
            (0, 1),
        ),
        (
            "C$N",
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Quadruple)),
            BondOrder::Quadruple,
            BondDirection::None,
            (0, 1),
        ),
        (
            "C->N",
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::DativeRight)),
            BondOrder::Dative,
            BondDirection::None,
            (0, 1),
        ),
        (
            "C<-N",
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::DativeLeft)),
            BondOrder::Dative,
            BondDirection::None,
            (1, 0),
        ),
        (
            "C/N",
            QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
                BondOrder::Single,
                BondOrder::Aromatic,
            ])),
            BondOrder::Single,
            BondDirection::EndUpRight,
            (0, 1),
        ),
        (
            "C\\N",
            QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
                BondOrder::Single,
                BondOrder::Aromatic,
            ])),
            BondOrder::Single,
            BondDirection::EndDownRight,
            (0, 1),
        ),
    ];

    for (smarts, expected_predicate, order, direction, endpoints) in cases {
        let graph = parse_source_case(smarts);
        assert_eq!((graph.num_atoms(), graph.num_bonds()), (2, 1), "{smarts}");
        let bond = graph.bond(0).expect("one source bond primitive");
        assert_eq!(bond.predicate(), &expected_predicate, "{smarts}");
        assert_eq!(bond.bond().order(), order, "{smarts}");
        assert_eq!(bond.bond().direction(), direction, "{smarts}");
        assert_eq!(bond.endpoints(), endpoints, "{smarts}");
        assert!(!bond.predicate_is_carrier_derived(), "{smarts}");
    }
}

#[test]
fn q03_common_bond_token_width_sets_following_bad_character_byte_position() {
    let raw_smarts = SmartsParseParams {
        allow_cxsmiles: false,
        parse_name: false,
        ..SmartsParseParams::default()
    };

    // setup_smarts_string trims trailing high-bit UTF-8 bytes on the pinned
    // signed-char build. The scanner therefore sees `C-` and reports EOS
    // while expecting the atom after the bond.
    assert_eq!(
        parse_smarts("C-☃", &raw_smarts)
            .expect_err("terminal snowman bytes are trimmed before lexing"),
        SmartsParseError::UnexpectedEnd("expected atom but reached end".to_owned())
    );

    for (input, position) in [
        ("C-☃C", 3),
        ("C=☃C", 3),
        ("C#☃C", 3),
        ("C:☃C", 3),
        ("C~☃C", 3),
        ("C$☃C", 3),
        ("C@☃C", 3),
        ("C/☃C", 3),
        ("C\\☃C", 3),
        ("C\\\\☃C", 4),
        ("C->☃C", 4),
        ("C<-☃C", 4),
    ] {
        assert_eq!(
            parse_smarts(input, &raw_smarts).expect_err("internal BAD_CHARACTER after bond token"),
            SmartsParseError::UnexpectedCharacter {
                position,
                character: '☃',
                context: "unexpected character in SMARTS string".to_owned(),
            },
            "{input:?}"
        );
    }
}

#[test]
fn q04_percent_ring_numbers_follow_pinned_grammar_and_error_position() {
    // RDKit 2026.03.1 smarts.yy accepts exactly two shorthand digits with a
    // nonzero first digit, or one through five grouped digits (including 0).
    // The grouped grammar's largest value is 99999, so its six-digit failure
    // is a syntax error at digit six rather than an integer overflow.
    for smarts in [
        "C%10CCCCC%10",
        "C%12CCCCC%12",
        "C%(0)CCCCC%(0)",
        "C%(00001)CCCCC%(00001)",
        "C%(99999)CCCCC%(99999)",
    ] {
        let graph = parse_source_case(smarts);
        assert_eq!((graph.num_atoms(), graph.num_bonds()), (6, 6), "{smarts}");
    }

    let raw_smarts = SmartsParseParams {
        allow_cxsmiles: false,
        parse_name: false,
        ..SmartsParseParams::default()
    };
    for (smarts, expected) in [
        (
            "C%",
            SmartsParseError::UnexpectedEnd("expected ring closure number".to_owned()),
        ),
        (
            "C%1",
            SmartsParseError::UnexpectedEnd("expected ring closure number".to_owned()),
        ),
        (
            "C%(1",
            SmartsParseError::UnexpectedEnd("expected ring closure number".to_owned()),
        ),
        (
            "C%0",
            SmartsParseError::UnexpectedCharacter {
                position: 3,
                character: '0',
                context: "invalid ring closure number".to_owned(),
            },
        ),
        (
            "C%00",
            SmartsParseError::UnexpectedCharacter {
                position: 3,
                character: '0',
                context: "invalid ring closure number".to_owned(),
            },
        ),
        (
            "C%()",
            SmartsParseError::UnexpectedCharacter {
                position: 4,
                character: ')',
                context: "invalid ring closure number".to_owned(),
            },
        ),
        (
            "C%1C",
            SmartsParseError::UnexpectedCharacter {
                position: 4,
                character: 'C',
                context: "invalid ring closure number".to_owned(),
            },
        ),
        (
            "C%(123456)CCCCC%(123456)",
            SmartsParseError::UnexpectedCharacter {
                position: 9,
                character: '6',
                context: "invalid ring closure number".to_owned(),
            },
        ),
        (
            "C%☃C",
            SmartsParseError::UnexpectedCharacter {
                position: 3,
                character: '☃',
                context: "unexpected character in SMARTS string".to_owned(),
            },
        ),
    ] {
        assert_eq!(
            parse_smarts(smarts, &raw_smarts).expect_err("invalid percent ring number"),
            expected,
            "{smarts:?}"
        );
    }
}

#[test]
fn q05_number_tokens_preserve_zero_width_and_signed_source_limit() {
    // A missing numeric token leaves the query's source default intact, while
    // an explicit zero is one token and any following digit remains available
    // to the enclosing atom-expression grammar.
    let ring = parse_source_case("[R]");
    let zero_rings = parse_source_case("[R0]");
    let zero_then_zero = parse_source_case("[R00]");
    let zero_then_one = parse_source_case("[R01]");
    let one_ring = parse_source_case("[R1]");
    assert_ne!(
        ring.atom(0).unwrap().predicate(),
        zero_rings.atom(0).unwrap().predicate()
    );
    assert_ne!(
        zero_rings.atom(0).unwrap().predicate(),
        zero_then_zero.atom(0).unwrap().predicate()
    );
    assert_ne!(
        zero_then_one.atom(0).unwrap().predicate(),
        one_ring.atom(0).unwrap().predicate()
    );

    // Each source [0] token reduces as a standalone number. The leading zero
    // therefore contributes a separate isotope predicate before 12C.
    let zero_prefixed_carbon = parse_source_case("[012C]");
    let twelve_carbon = parse_source_case("[12C]");
    assert_ne!(
        zero_prefixed_carbon.atom(0).unwrap().predicate(),
        twelve_carbon.atom(0).unwrap().predicate()
    );

    // The pinned int32 guard accepts 1,000,000,000 and 2,147,483,639, then
    // rejects 2,147,483,640 while folding its final digit.
    for (smarts, value) in [
        ("[R1000000000]", 1_000_000_000),
        ("[R2147483639]", 2_147_483_639),
    ] {
        let graph = parse_source_case(smarts);
        assert_eq!(
            graph.atom(0).unwrap().predicate(),
            &QueryNode::Predicate(AtomQueryPredicate::NumAtomRings(value)),
            "{smarts}"
        );
    }

    let raw_smarts = SmartsParseParams {
        allow_cxsmiles: false,
        parse_name: false,
        ..SmartsParseParams::default()
    };
    assert!(matches!(
        parse_smarts("[R2147483640]", &raw_smarts)
            .expect_err("source signed-int accumulation rejects this value"),
        SmartsParseError::InvalidAtomPrimitive { detail, .. }
            if detail == "number too large"
    ));
}

#[test]
fn q07e_ring_primitives_keep_source_identity_full_targets_and_defaults() {
    // Pinned smarts.ll creates AtomRingQuery(-1) for omitted R, smarts.yy
    // writes a numeric target into that same query and uses int factories for
    // r/k/x. SmartsWrite.cpp emits the corresponding token and target.
    let defaults = [
        ("[R]", AtomQueryPredicate::NumAtomRings(-1), "[R]"),
        ("[r]", AtomQueryPredicate::InRing, "[R]"),
        ("[k]", AtomQueryPredicate::InRing, "[R]"),
        ("[x]", AtomQueryPredicate::HasRingBond, "[x]"),
        ("[R0]", AtomQueryPredicate::NumAtomRings(0), "[R0]"),
        ("[r0]", AtomQueryPredicate::SmallestRingSize(0), "[r0]"),
        ("[k0]", AtomQueryPredicate::InRingOfSize(0), "[k0]"),
        ("[x0]", AtomQueryPredicate::RingBondCount(0), "[x0]"),
    ];
    let baseline_carrier = parse_source_case("[R]")
        .atom(0)
        .expect("default R query atom")
        .try_to_atom()
        .expect("default R carrier projection");
    for (smarts, predicate, expected_writer) in defaults {
        let graph = parse_source_case(smarts);
        let atom = graph.atom(0).expect("single ring query atom");
        assert_eq!(
            atom.predicate(),
            &QueryNode::Predicate(predicate),
            "{smarts}"
        );
        assert_eq!(atom.try_to_atom().unwrap(), baseline_carrier, "{smarts}");
        assert!(!atom.predicate_is_carrier_derived(), "{smarts}");
        assert_eq!(
            write_smarts(&graph, &SmartsWriteParams::default()).unwrap(),
            expected_writer,
            "{smarts}"
        );
    }

    for (token, target) in [
        ("R", 255),
        ("R", 256),
        ("R", 2_147_483_639),
        ("r", 255),
        ("r", 256),
        ("r", 2_147_483_639),
        ("k", 255),
        ("k", 256),
        ("k", 2_147_483_639),
        ("x", 255),
        ("x", 256),
        ("x", 2_147_483_639),
    ] {
        let smarts = format!("[{token}{target}]");
        let predicate = match token {
            "R" => AtomQueryPredicate::NumAtomRings(target),
            "r" => AtomQueryPredicate::SmallestRingSize(target),
            "k" => AtomQueryPredicate::InRingOfSize(target),
            "x" => AtomQueryPredicate::RingBondCount(target),
            _ => unreachable!("fixed ring token table"),
        };
        let graph = parse_source_case(&smarts);
        let atom = graph.atom(0).expect("single numeric ring query atom");
        assert_eq!(
            atom.predicate(),
            &QueryNode::Predicate(predicate),
            "{smarts}"
        );
        assert_eq!(atom.try_to_atom().unwrap(), baseline_carrier, "{smarts}");
        assert!(!atom.predicate_is_carrier_derived(), "{smarts}");
        assert_eq!(
            write_smarts(&graph, &SmartsWriteParams::default()).unwrap(),
            smarts,
            "{smarts}"
        );
    }

    let raw_smarts = SmartsParseParams {
        allow_cxsmiles: false,
        parse_name: false,
        ..SmartsParseParams::default()
    };
    for token in ["R", "r", "k", "x"] {
        let smarts = format!("[{token}2147483640]");
        assert_eq!(
            parse_smarts(&smarts, &raw_smarts).expect_err("source int guard rejects overflow"),
            SmartsParseError::InvalidAtomPrimitive {
                position: 1,
                detail: "number too large".to_owned(),
            },
            "{smarts}"
        );
    }

    let ring_atoms = (0..256)
        .map(|index| Atom::from_spec(AtomId::new(index), AtomSpec::new(Element::C)))
        .collect::<Vec<_>>();
    let ring_bonds = (0..256)
        .map(|index| {
            Bond::from_spec(
                BondId::new(index),
                BondSpec::new(
                    AtomId::new(index),
                    AtomId::new((index + 1) % 256),
                    BondOrder::Single,
                ),
            )
        })
        .collect::<Vec<_>>();
    let ring_256 = TopologyBlock::try_from_parts(ring_atoms, ring_bonds, Vec::new(), Vec::new())
        .expect("256-member ring topology");
    let acyclic = TopologyBlock::try_from_parts(
        vec![Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C))],
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .expect("acyclic control topology");

    let matches = |smarts: &str, target: &TopologyBlock| {
        !match_query(&parse_source_case(smarts), target)
            .expect("ring primitive query matching")
            .is_empty()
    };
    assert!(matches("[R]", &ring_256));
    assert!(!matches("[R]", &acyclic));
    assert!(!matches("[R256]", &acyclic));
    assert!(matches("[R0]", &acyclic));
    assert!(!matches("[R0]", &ring_256));
    assert!(matches("[r256]", &ring_256));
    assert!(!matches("[r256]", &acyclic));
    assert!(matches("[k256]", &ring_256));
    assert!(matches("[x0]", &acyclic));
    assert!(!matches("[x256]", &acyclic));
}

#[test]
fn q07e_ring_ranges_keep_source_query_classes_bounds_and_matching() {
    let baseline_carrier = parse_source_case("[R]")
        .atom(0)
        .expect("default R query atom")
        .try_to_atom()
        .expect("default R carrier projection");
    let range_data_function = |token: &str, lower: Option<i32>, upper: Option<i32>| match token {
        "R" => AtomRangeDataFunction::NumAtomRings,
        "r" => AtomRangeDataFunction::MinRingSize,
        "x" => AtomRangeDataFunction::RingBondCount,
        "k" => AtomRangeDataFunction::AtomRingSize {
            lower: lower.unwrap_or(-1),
            upper: upper.unwrap_or(-1),
            lower_open: false,
            upper_open: false,
        },
        _ => unreachable!("fixed ring token table"),
    };
    let assert_range = |smarts: &str,
                        expected_bounds: AtomRangeBounds,
                        expected_data_function: AtomRangeDataFunction| {
        let graph = parse_source_case(smarts);
        let atom = graph.atom(0).expect("single ring-range query atom");
        let QueryNode::Predicate(AtomQueryPredicate::Range(range)) = atom.predicate() else {
            panic!("{smarts} must lower to the generic range leaf");
        };
        assert_eq!(range.bounds(), expected_bounds, "{smarts}");
        assert_eq!(range.data_function(), expected_data_function, "{smarts}");
        assert_eq!(atom.try_to_atom().unwrap(), baseline_carrier, "{smarts}");
        assert!(!atom.predicate_is_carrier_derived(), "{smarts}");
        let written = write_smarts(&graph, &SmartsWriteParams::default()).unwrap();
        assert_eq!(written, smarts, "{smarts}");
        let reparsed = parse_source_case(&written);
        assert_eq!(
            reparsed.atom(0).unwrap().predicate(),
            atom.predicate(),
            "{smarts} source class and k data survive write/reparse"
        );
    };

    let maximum = 2_147_483_639;
    for token in ["R", "r", "k", "x"] {
        for value in [0, 255, 256, maximum] {
            let open_lower = format!("[{token}{{-{value}}}]");
            assert_range(
                &open_lower,
                AtomRangeBounds::GreaterEqual(value),
                range_data_function(token, None, Some(value)),
            );

            let open_upper = format!("[{token}{{{value}-}}]");
            assert_range(
                &open_upper,
                AtomRangeBounds::LessEqual(value),
                range_data_function(token, Some(value), None),
            );
        }

        for (lower, upper) in [(0, 0), (255, 256), (256, 255), (maximum, maximum)] {
            let smarts = format!("[{token}{{{lower}-{upper}}}]");
            assert_range(
                &smarts,
                AtomRangeBounds::Inclusive {
                    lower,
                    upper,
                    lower_open: false,
                    upper_open: false,
                },
                range_data_function(token, Some(lower), Some(upper)),
            );
        }
    }

    let raw_smarts = SmartsParseParams {
        allow_cxsmiles: false,
        parse_name: false,
        ..SmartsParseParams::default()
    };
    for token in ["R", "r", "k", "x"] {
        for (range, position) in [
            ("{2147483640-}", 2),
            ("{-2147483640}", 3),
            ("{1-2147483640}", 4),
        ] {
            let smarts = format!("[{token}{range}]");
            assert_eq!(
                parse_smarts(&smarts, &raw_smarts)
                    .expect_err("source signed-int accumulation rejects overflow"),
                SmartsParseError::InvalidAtomPrimitive {
                    position,
                    detail: "number too large".to_owned(),
                },
                "{smarts}"
            );
        }

        for range in ["{}", "{-}", "{1}", "{1--2}", "{1-2", "{1-2-3}"] {
            let smarts = format!("[{token}{range}]");
            let error = parse_smarts(&smarts, &raw_smarts)
                .expect_err("range form is absent from the pinned ring grammar");
            assert!(
                matches!(
                    error,
                    SmartsParseError::UnexpectedEnd(_)
                        | SmartsParseError::InvalidAtomPrimitive { .. }
                ),
                "{smarts}: {error:?}"
            );
        }
    }

    let ring_atoms = (0..256)
        .map(|index| Atom::from_spec(AtomId::new(index), AtomSpec::new(Element::C)))
        .collect::<Vec<_>>();
    let ring_bonds = (0..256)
        .map(|index| {
            Bond::from_spec(
                BondId::new(index),
                BondSpec::new(
                    AtomId::new(index),
                    AtomId::new((index + 1) % 256),
                    BondOrder::Single,
                ),
            )
        })
        .collect::<Vec<_>>();
    let ring_256 = TopologyBlock::try_from_parts(ring_atoms, ring_bonds, Vec::new(), Vec::new())
        .expect("256-member ring topology");
    let acyclic = TopologyBlock::try_from_parts(
        vec![Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C))],
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .expect("acyclic control topology");
    let matches = |smarts: &str, target: &TopologyBlock| {
        !match_query(&parse_source_case(smarts), target)
            .expect("ring range query matching")
            .is_empty()
    };

    for (smarts, ring_match, acyclic_match) in [
        // GreaterEqualQuery compares stored threshold >= observed; LessEqual
        // compares stored threshold <= observed (pinned Query headers).
        ("[R{-0}]", false, true),
        ("[R{-1}]", true, true),
        ("[R{-2}]", true, true),
        ("[R{0-}]", true, true),
        ("[R{1-}]", true, false),
        ("[R{2-}]", false, false),
        ("[R{1-1}]", true, false),
        ("[r{-255}]", false, true),
        ("[r{-256}]", true, true),
        ("[r{-257}]", true, true),
        ("[r{255-}]", true, false),
        ("[r{256-}]", true, false),
        ("[r{257-}]", false, false),
        ("[r{256-256}]", true, false),
        ("[k{-255}]", false, false),
        ("[k{-256}]", true, false),
        ("[k{-257}]", true, false),
        ("[k{255-}]", true, false),
        ("[k{256-}]", true, false),
        ("[k{257-}]", false, false),
        ("[k{256-256}]", true, false),
        ("[x{-1}]", false, true),
        ("[x{-2}]", true, true),
        ("[x{-3}]", true, true),
        ("[x{1-}]", true, false),
        ("[x{2-}]", true, false),
        ("[x{3-}]", false, false),
        ("[x{2-2}]", true, false),
        ("[!R{-1}]", false, false),
        ("[!r{256-}]", false, true),
        ("[!k{-256}]", false, true),
        ("[!x{2-}]", false, true),
        ("[D{-1}]", false, true),
        ("[D{-2}]", true, true),
        ("[D{-3}]", true, true),
        ("[D{1-}]", true, false),
        ("[D{2-}]", true, false),
        ("[D{3-}]", false, false),
        ("[D{2-2}]", true, false),
        ("[!D{-2}]", false, false),
    ] {
        assert_eq!(matches(smarts, &ring_256), ring_match, "{smarts} in ring");
        assert_eq!(
            matches(smarts, &acyclic),
            acyclic_match,
            "{smarts} in acyclic control"
        );
    }
}

#[test]
fn q07e_cx_ring_bond_lowering_keeps_exact_and_scan_targets_and_carrier() {
    fn contains_ring_bond_target(node: &QueryNode<AtomQueryPredicate>, target: i32) -> bool {
        match node {
            QueryNode::Predicate(AtomQueryPredicate::RingBondCount(value)) => *value == target,
            QueryNode::And(children) | QueryNode::Or(children) => children
                .iter()
                .any(|child| contains_ring_bond_target(child, target)),
            _ => false,
        }
    }

    let baseline_carrier = parse_source_case("C")
        .atom(0)
        .expect("ordinary carbon query atom")
        .try_to_atom()
        .expect("ordinary carbon carrier");
    for (smarts, target) in [("C |rb:0:3|", 3), ("C |rb:0:*|", 0xDEAD_BEEF_u32 as i32)] {
        let graph = parse_source_case(smarts);
        let atom = graph.atom(0).expect("CX ring-bond query atom");
        assert!(
            contains_ring_bond_target(atom.predicate(), target),
            "{smarts}"
        );
        assert_eq!(atom.try_to_atom().unwrap(), baseline_carrier, "{smarts}");
        assert!(!atom.predicate_is_carrier_derived(), "{smarts}");
    }
}

#[test]
fn q08_atom_boolean_precedence_and_not_carrier_reductions() {
    let c = QueryNode::predicate(AtomQueryPredicate::AtomType {
        atomic_number: 6,
        aromatic: false,
    });
    let n = QueryNode::predicate(AtomQueryPredicate::AtomType {
        atomic_number: 7,
        aromatic: false,
    });
    let o = QueryNode::predicate(AtomQueryPredicate::AtomType {
        atomic_number: 8,
        aromatic: false,
    });

    for (smarts, expected) in [
        (
            "[C,N&O]",
            QueryNode::or(vec![c.clone(), QueryNode::and(vec![n.clone(), o.clone()])]),
        ),
        (
            "[C&N,O]",
            QueryNode::or(vec![QueryNode::and(vec![c.clone(), n.clone()]), o.clone()]),
        ),
        (
            "[C;N,O]",
            QueryNode::and(vec![c.clone(), QueryNode::or(vec![n.clone(), o.clone()])]),
        ),
        (
            "[C,N;O]",
            QueryNode::and(vec![QueryNode::or(vec![c, n]), o]),
        ),
    ] {
        let graph = parse_source_case(smarts);
        let atom = graph.atom(0).expect("one Boolean-expression atom");
        assert_eq!(atom.predicate(), &expected, "{smarts}");
        assert!(!atom.predicate_is_carrier_derived(), "{smarts}");
    }

    let isotope_or = parse_source_case("[13C,N&O]");
    let isotope_or_atom = isotope_or.atom(0).expect("OR-reduced isotope query");
    assert_eq!(
        isotope_or_atom.predicate(),
        &QueryNode::or(vec![
            QueryNode::and(vec![
                QueryNode::predicate(AtomQueryPredicate::AtomType {
                    atomic_number: 6,
                    aromatic: false,
                }),
                QueryNode::predicate(AtomQueryPredicate::Isotope(13)),
            ]),
            QueryNode::and(vec![
                QueryNode::predicate(AtomQueryPredicate::AtomType {
                    atomic_number: 7,
                    aromatic: false,
                }),
                QueryNode::predicate(AtomQueryPredicate::AtomType {
                    atomic_number: 8,
                    aromatic: false,
                }),
            ]),
        ])
    );
    assert_eq!(
        isotope_or_atom.identity(),
        QueryAtomIdentity::Element(Element::DUMMY)
    );
    assert_eq!(isotope_or_atom.isotope(), None);
    assert_eq!(isotope_or_atom.formal_charge(), 0);
    assert_eq!(isotope_or_atom.explicit_hydrogens(), 0);
    assert!(!isotope_or_atom.predicate_is_carrier_derived());

    let isotope_semi = parse_source_case("[13C;N,O]");
    let isotope_semi_atom = isotope_semi.atom(0).expect("SEMI-reduced isotope query");
    assert_eq!(
        isotope_semi_atom.predicate(),
        &QueryNode::and(vec![
            QueryNode::and(vec![
                QueryNode::predicate(AtomQueryPredicate::AtomType {
                    atomic_number: 6,
                    aromatic: false,
                }),
                QueryNode::predicate(AtomQueryPredicate::Isotope(13)),
            ]),
            QueryNode::or(vec![
                QueryNode::predicate(AtomQueryPredicate::AtomType {
                    atomic_number: 7,
                    aromatic: false,
                }),
                QueryNode::predicate(AtomQueryPredicate::AtomType {
                    atomic_number: 8,
                    aromatic: false,
                }),
            ]),
        ])
    );
    assert_eq!(
        isotope_semi_atom.identity(),
        QueryAtomIdentity::Element(Element::C)
    );
    assert_eq!(isotope_semi_atom.isotope(), None);
    assert_eq!(isotope_semi_atom.formal_charge(), 0);
    assert_eq!(isotope_semi_atom.explicit_hydrogens(), 0);
    assert!(!isotope_semi_atom.predicate_is_carrier_derived());

    let not_cases = [
        (
            "13C",
            QueryNode::and(vec![
                QueryNode::predicate(AtomQueryPredicate::AtomType {
                    atomic_number: 6,
                    aromatic: false,
                }),
                QueryNode::predicate(AtomQueryPredicate::Isotope(13)),
            ]),
            false,
        ),
        (
            "+",
            QueryNode::predicate(AtomQueryPredicate::FormalCharge(1)),
            false,
        ),
        (
            "H",
            QueryNode::predicate(AtomQueryPredicate::HydrogenCount(1)),
            true,
        ),
    ];
    for (primitive, base_predicate, retains_no_implicit) in not_cases {
        for negation_count in 1..=3 {
            let bangs = "!".repeat(negation_count);
            let smarts = format!("[{bangs}{primitive}]");
            let expected_predicate = if negation_count % 2 == 1 {
                QueryNode::not(base_predicate.clone())
            } else {
                base_predicate.clone()
            };
            let graph = parse_source_case(&smarts);
            let atom = graph.atom(0).expect("one repeatedly negated query atom");
            assert_eq!(atom.predicate(), &expected_predicate, "{smarts}");
            assert_eq!(
                atom.identity(),
                QueryAtomIdentity::Element(Element::DUMMY),
                "{smarts}"
            );
            assert_eq!(atom.isotope(), None, "{smarts}");
            assert_eq!(atom.formal_charge(), 0, "{smarts}");
            assert_eq!(atom.explicit_hydrogens(), 0, "{smarts}");
            assert_eq!(atom.no_implicit(), retains_no_implicit, "{smarts}");
            assert!(!atom.predicate_is_carrier_derived(), "{smarts}");
        }
    }
}

#[test]
fn q06_range_bounds_cover_inclusive_open_ends_and_malformed_forms() {
    // Pinned smarts.yy has two-sided ranges plus missing-lower and
    // missing-upper reductions; two-sided bounds use false/false open flags.
    for (smarts, expected_bounds) in [
        (
            "[D{2-4}]",
            AtomRangeBounds::Inclusive {
                lower: 2,
                upper: 4,
                lower_open: false,
                upper_open: false,
            },
        ),
        ("[D{-2}]", AtomRangeBounds::GreaterEqual(2)),
        ("[D{2-}]", AtomRangeBounds::LessEqual(2)),
        (
            "[D{0-2}]",
            AtomRangeBounds::Inclusive {
                lower: 0,
                upper: 2,
                lower_open: false,
                upper_open: false,
            },
        ),
        ("[D{-0}]", AtomRangeBounds::GreaterEqual(0)),
    ] {
        let graph = parse_source_case(smarts);
        let predicate = graph.atom(0).expect("one range-query atom").predicate();
        let QueryNode::Predicate(AtomQueryPredicate::Range(range)) = predicate else {
            panic!("{smarts} must produce an atom range predicate, got {predicate:?}");
        };
        assert_eq!(range.bounds(), expected_bounds, "{smarts}");
        assert_eq!(
            range.data_function(),
            AtomRangeDataFunction::ExplicitDegree,
            "{smarts}"
        );
    }

    let raw_smarts = SmartsParseParams {
        allow_cxsmiles: false,
        parse_name: false,
        ..SmartsParseParams::default()
    };
    for smarts in [
        "[D{}]",
        "[D{2}]",
        "[D{2--4}]",
        "[D{-}]",
        "[D{2-4]",
        "[D{2-x}]",
    ] {
        let error = parse_smarts(smarts, &raw_smarts).expect_err("malformed atom range");
        assert!(
            matches!(
                &error,
                SmartsParseError::UnexpectedEnd(_) | SmartsParseError::InvalidAtomPrimitive { .. }
            ),
            "{smarts:?}: {error:?}"
        );
    }
}

#[test]
fn q07c_aromatic_flags_and_hybridization_tokens_keep_query_carriers_separate() {
    // Pinned smarts.ll gives generic a/A distinct query predicates, with only
    // `a` setting the carrier aromatic flag. simple_atom element tokens keep
    // their element identity and separately project the source aromatic bit.
    for (smarts, expected, atomic_number, aromatic) in [
        ("a", AtomQueryPredicate::IsAromatic(true), 0, true),
        ("[a]", AtomQueryPredicate::IsAromatic(true), 0, true),
        ("[A]", AtomQueryPredicate::IsAromatic(false), 0, false),
        (
            "[c]",
            AtomQueryPredicate::AtomType {
                atomic_number: 6,
                aromatic: true,
            },
            6,
            true,
        ),
        (
            "[C]",
            AtomQueryPredicate::AtomType {
                atomic_number: 6,
                aromatic: false,
            },
            6,
            false,
        ),
    ] {
        let graph = parse_source_case(smarts);
        let atom = graph.atom(0).expect("one aromatic/aliphatic query atom");
        assert_eq!(
            atom.predicate(),
            &QueryNode::predicate(expected),
            "{smarts}"
        );
        assert_eq!(atom.atomic_number(), atomic_number, "{smarts}");
        assert_eq!(atom.is_aromatic(), aromatic, "{smarts}");
        assert_eq!(atom.hybridization(), Hybridization::Unspecified, "{smarts}");
        assert!(!atom.predicate_is_carrier_derived(), "{smarts}");
    }

    // smarts.ll's six literal HYB_TOKEN actions map directly to RDKit's
    // HybridizationType values and leave QueryAtom's carrier unspecified.
    for (smarts, expected) in [
        ("[^0]", Hybridization::S),
        ("[^1]", Hybridization::Sp),
        ("[^2]", Hybridization::Sp2),
        ("[^3]", Hybridization::Sp3),
        ("[^4]", Hybridization::Sp3d),
        ("[^5]", Hybridization::Sp3d2),
    ] {
        let graph = parse_source_case(smarts);
        let atom = graph.atom(0).expect("one hybridization query atom");
        assert_eq!(
            atom.predicate(),
            &QueryNode::predicate(AtomQueryPredicate::HybridizationMatch(expected)),
            "{smarts}"
        );
        assert_eq!(atom.hybridization(), Hybridization::Unspecified, "{smarts}");
        assert_eq!(atom.atomic_number(), 0, "{smarts}");
        assert!(!atom.is_aromatic(), "{smarts}");
        assert!(!atom.predicate_is_carrier_derived(), "{smarts}");
    }
}

#[test]
fn q07d_degree_targets_keep_full_source_i32_and_default_carriers() {
    // smarts.ll creates D/X/v with query target 1; smarts.yy's shared
    // COMPLEX_ATOM_QUERY_TOKEN number reduction sets the original guarded
    // signed-int value. The pinned RDKit 2026.03.1 oracle preserves 256 and
    // 2147483639, and rejects 2147483640.
    let tokens: [(&str, fn(i32) -> AtomQueryPredicate); 3] = [
        ("D", AtomQueryPredicate::ExplicitDegree),
        ("X", AtomQueryPredicate::TotalDegree),
        ("v", AtomQueryPredicate::TotalValence),
    ];
    let targets = [
        ("", 1),
        ("0", 0),
        ("255", 255),
        ("256", 256),
        ("2147483639", 2_147_483_639),
    ];

    for (token, make_predicate) in tokens {
        for (suffix, target) in targets {
            let smarts = format!("[{token}{suffix}]");
            let graph = parse_source_case(&smarts);
            let atom = graph.atom(0).expect("one D/X/v query atom");
            assert_eq!(
                atom.predicate(),
                &QueryNode::predicate(make_predicate(target)),
                "{smarts}"
            );
            assert!(!atom.predicate_is_carrier_derived(), "{smarts}");
            assert_eq!(atom.atomic_number(), 0, "{smarts}");
            assert!(!atom.is_aromatic(), "{smarts}");
            assert_eq!(atom.formal_charge(), 0, "{smarts}");
            assert_eq!(atom.isotope(), None, "{smarts}");
            assert_eq!(atom.explicit_hydrogens(), 0, "{smarts}");
            assert!(!atom.no_implicit(), "{smarts}");
        }
    }

    let raw_smarts = SmartsParseParams {
        allow_cxsmiles: false,
        parse_name: false,
        ..SmartsParseParams::default()
    };
    for token in ["D", "X", "v"] {
        let smarts = format!("[{token}2147483640]");
        assert_eq!(
            parse_smarts(&smarts, &raw_smarts).expect_err("source signed-int number guard"),
            SmartsParseError::InvalidAtomPrimitive {
                position: 1,
                detail: "number too large".to_owned(),
            },
            "{smarts}"
        );
    }
}

#[test]
fn q07d_hydrogen_targets_keep_full_source_i32_and_narrow_carriers() {
    // smarts.ll emits `h` as IMPLICIT_H_ATOM_QUERY_TOKEN; smarts.yy maps it
    // without a number to HasImplicitHydrogen and with a number to the full
    // int ImplicitHydrogenCount target. H_TOKEN reductions build a full int
    // HydrogenCount target while Atom::setNumExplicitHs(unsigned int) stores
    // into its independent uint8_t carrier.
    let implicit_default = parse_source_case("[h]");
    let implicit_default_atom = implicit_default.atom(0).expect("one implicit-H query atom");
    assert_eq!(
        implicit_default_atom.predicate(),
        &QueryNode::predicate(AtomQueryPredicate::HasImplicitHydrogen)
    );
    assert_eq!(implicit_default_atom.atomic_number(), 0);
    assert_eq!(implicit_default_atom.isotope(), None);
    assert_eq!(implicit_default_atom.formal_charge(), 0);
    assert_eq!(implicit_default_atom.explicit_hydrogens(), 0);
    assert!(!implicit_default_atom.no_implicit());
    assert!(!implicit_default_atom.predicate_is_carrier_derived());

    for (suffix, target) in [
        ("0", 0),
        ("255", 255),
        ("256", 256),
        ("2147483639", 2_147_483_639),
    ] {
        let smarts = format!("[h{suffix}]");
        let graph = parse_source_case(&smarts);
        let atom = graph.atom(0).expect("one numbered implicit-H query atom");
        assert_eq!(
            atom.predicate(),
            &QueryNode::predicate(AtomQueryPredicate::ImplicitHydrogenCount(target)),
            "{smarts}"
        );
        assert_eq!(atom.atomic_number(), 0, "{smarts}");
        assert_eq!(atom.isotope(), None, "{smarts}");
        assert_eq!(atom.formal_charge(), 0, "{smarts}");
        assert_eq!(atom.explicit_hydrogens(), 0, "{smarts}");
        assert!(!atom.no_implicit(), "{smarts}");
        assert!(!atom.predicate_is_carrier_derived(), "{smarts}");
    }

    for (suffix, target, carrier) in [
        ("0", 0, 0),
        ("255", 255, 255),
        ("256", 256, 0),
        ("2147483639", 2_147_483_639, 247),
    ] {
        let smarts = format!("[H{suffix}]");
        let graph = parse_source_case(&smarts);
        let atom = graph.atom(0).expect("one numbered H-count query atom");
        assert_eq!(
            atom.predicate(),
            &QueryNode::predicate(AtomQueryPredicate::HydrogenCount(target)),
            "{smarts}"
        );
        assert_eq!(atom.atomic_number(), 0, "{smarts}");
        assert_eq!(atom.isotope(), None, "{smarts}");
        assert_eq!(atom.formal_charge(), 0, "{smarts}");
        assert_eq!(atom.explicit_hydrogens(), carrier, "{smarts}");
        assert!(atom.no_implicit(), "{smarts}");
        assert!(!atom.predicate_is_carrier_derived(), "{smarts}");
    }

    for (suffix, target, carrier) in [
        ("0", 0, 0),
        ("255", 255, 255),
        ("256", 256, 0),
        ("2147483639", 2_147_483_639, 247),
    ] {
        let smarts = format!("[13H{suffix}]");
        let graph = parse_source_case(&smarts);
        let atom = graph.atom(0).expect("one isotope-plus-H-count query atom");
        assert_eq!(
            atom.predicate(),
            &QueryNode::and(vec![
                QueryNode::predicate(AtomQueryPredicate::Isotope(13)),
                QueryNode::predicate(AtomQueryPredicate::HydrogenCount(target)),
            ]),
            "{smarts}"
        );
        assert_eq!(atom.atomic_number(), 0, "{smarts}");
        assert_eq!(atom.isotope(), Some(13), "{smarts}");
        assert_eq!(atom.formal_charge(), 0, "{smarts}");
        assert_eq!(atom.explicit_hydrogens(), carrier, "{smarts}");
        assert!(atom.no_implicit(), "{smarts}");
        assert!(!atom.predicate_is_carrier_derived(), "{smarts}");
    }

    // These source productions stay distinct from numbered H-count queries:
    // bare [H] and [13H] use hydrogen_atom; ordinary [CH] uses H_TOKEN in a
    // point query and therefore has an HCount(1) target plus carrier writes.
    let bare_hydrogen = parse_source_case("[H]");
    let bare_hydrogen_atom = bare_hydrogen.atom(0).expect("atomic hydrogen");
    assert_eq!(
        bare_hydrogen_atom.predicate(),
        &QueryNode::predicate(AtomQueryPredicate::AtomicNumber(1))
    );
    assert_eq!(bare_hydrogen_atom.explicit_hydrogens(), 0);
    assert!(!bare_hydrogen_atom.no_implicit());

    let isotopic_hydrogen = parse_source_case("[13H]");
    let isotopic_hydrogen_atom = isotopic_hydrogen.atom(0).expect("isotopic atomic hydrogen");
    assert_eq!(
        isotopic_hydrogen_atom.predicate(),
        &QueryNode::and(vec![
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(1)),
            QueryNode::predicate(AtomQueryPredicate::Isotope(13)),
        ])
    );
    assert_eq!(isotopic_hydrogen_atom.explicit_hydrogens(), 0);
    assert!(!isotopic_hydrogen_atom.no_implicit());

    let ordinary_hydrogen = parse_source_case("[CH]");
    let ordinary_hydrogen_atom = ordinary_hydrogen.atom(0).expect("ordinary carbon H query");
    assert_eq!(
        ordinary_hydrogen_atom.predicate(),
        &QueryNode::and(vec![
            QueryNode::predicate(AtomQueryPredicate::AtomType {
                atomic_number: 6,
                aromatic: false,
            }),
            QueryNode::predicate(AtomQueryPredicate::HydrogenCount(1)),
        ])
    );
    assert_eq!(ordinary_hydrogen_atom.atomic_number(), 6);
    assert_eq!(ordinary_hydrogen_atom.explicit_hydrogens(), 1);
    assert!(ordinary_hydrogen_atom.no_implicit());

    let raw_smarts = SmartsParseParams {
        allow_cxsmiles: false,
        parse_name: false,
        ..SmartsParseParams::default()
    };
    for (smarts, position) in [
        ("[h2147483640]", 1),
        ("[H2147483640]", 1),
        ("[13H2147483640]", 3),
    ] {
        assert_eq!(
            parse_smarts(smarts, &raw_smarts).expect_err("source signed-int number guard"),
            SmartsParseError::InvalidAtomPrimitive {
                position,
                detail: "number too large".to_owned(),
            },
            "{smarts}"
        );
    }
}

#[test]
fn q09_hydrogen_atom_and_h_count_forms_follow_distinct_source_reductions() {
    let atomic_hydrogen_cases = [
        (
            "[H]",
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(1)),
            None,
            0,
        ),
        (
            "[13H]",
            QueryNode::and(vec![
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(1)),
                QueryNode::predicate(AtomQueryPredicate::Isotope(13)),
            ]),
            Some(13),
            0,
        ),
        (
            "[H+]",
            QueryNode::and(vec![
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(1)),
                QueryNode::predicate(AtomQueryPredicate::FormalCharge(1)),
            ]),
            None,
            1,
        ),
        (
            "[13H+]",
            QueryNode::and(vec![
                QueryNode::and(vec![
                    QueryNode::predicate(AtomQueryPredicate::AtomicNumber(1)),
                    QueryNode::predicate(AtomQueryPredicate::Isotope(13)),
                ]),
                QueryNode::predicate(AtomQueryPredicate::FormalCharge(1)),
            ]),
            Some(13),
            1,
        ),
    ];
    for (smarts, expected_predicate, isotope, formal_charge) in atomic_hydrogen_cases {
        let graph = parse_source_case(smarts);
        let atom = graph.atom(0).expect("one atomic-hydrogen query");
        assert_eq!(
            atom.identity(),
            QueryAtomIdentity::from_atomic_number(1),
            "{smarts}"
        );
        assert_eq!(atom.predicate(), &expected_predicate, "{smarts}");
        assert_eq!(atom.isotope(), isotope, "{smarts}");
        assert_eq!(atom.formal_charge(), formal_charge, "{smarts}");
        assert_eq!(atom.explicit_hydrogens(), 0, "{smarts}");
        assert!(!atom.no_implicit(), "{smarts}");
        assert!(!atom.predicate_is_carrier_derived(), "{smarts}");
    }

    // A numeric suffix after H is the H_TOKEN count production, not an
    // atomic-hydrogen production. The source sets the count carrier and
    // noImplicit flag while leaving atomic identity at the query default.
    for (smarts, expected_predicate, isotope) in [
        (
            "[H1]",
            QueryNode::predicate(AtomQueryPredicate::HydrogenCount(1)),
            None,
        ),
        (
            "[13H1]",
            QueryNode::and(vec![
                QueryNode::predicate(AtomQueryPredicate::Isotope(13)),
                QueryNode::predicate(AtomQueryPredicate::HydrogenCount(1)),
            ]),
            Some(13),
        ),
    ] {
        let graph = parse_source_case(smarts);
        let atom = graph.atom(0).expect("one numbered H-count query");
        assert_eq!(
            atom.identity(),
            QueryAtomIdentity::Element(Element::DUMMY),
            "{smarts}"
        );
        assert_eq!(atom.predicate(), &expected_predicate, "{smarts}");
        assert_eq!(atom.isotope(), isotope, "{smarts}");
        assert_eq!(atom.formal_charge(), 0, "{smarts}");
        assert_eq!(atom.explicit_hydrogens(), 1, "{smarts}");
        assert!(atom.no_implicit(), "{smarts}");
        assert!(!atom.predicate_is_carrier_derived(), "{smarts}");
    }

    let carbon_hydrogen = parse_source_case("[CH]");
    let carbon_hydrogen_atom = carbon_hydrogen
        .atom(0)
        .expect("carbon with a hydrogen-count point query");
    assert_eq!(
        carbon_hydrogen_atom.predicate(),
        &QueryNode::and(vec![
            QueryNode::predicate(AtomQueryPredicate::AtomType {
                atomic_number: 6,
                aromatic: false,
            }),
            QueryNode::predicate(AtomQueryPredicate::HydrogenCount(1)),
        ])
    );
    assert_eq!(carbon_hydrogen_atom.atomic_number(), 6);
    assert_eq!(carbon_hydrogen_atom.explicit_hydrogens(), 1);
    assert!(carbon_hydrogen_atom.no_implicit());
    assert!(!carbon_hydrogen_atom.predicate_is_carrier_derived());
}

#[test]
fn q21_hydrogen_atom_and_count_queries_keep_source_classification() {
    let merge_hs = SmartsParseParams {
        merge_hs: true,
        ..SmartsParseParams::default()
    };

    for smarts in ["[C][H]", "[C][#1]"] {
        let graph = parse_smarts(smarts, &merge_hs)
            .unwrap_or_else(|error| panic!("pinned RDKit SMARTS {smarts:?}: {error}"));
        assert_eq!(graph.num_atoms(), 1, "{smarts}");
    }

    for smarts in ["[C][H1]", "[C][13H1]"] {
        let graph = parse_smarts(smarts, &merge_hs)
            .unwrap_or_else(|error| panic!("pinned RDKit SMARTS {smarts:?}: {error}"));
        assert_eq!(graph.num_atoms(), 2, "{smarts}");
    }

    let count_query = parse_smarts("[C][H1]", &merge_hs).expect("hydrogen-count query");
    assert_eq!(
        count_query
            .atom(1)
            .expect("hydrogen-count atom")
            .predicate(),
        &QueryNode::predicate(AtomQueryPredicate::HydrogenCount(1))
    );
    let isotope_count_query =
        parse_smarts("[C][13H1]", &merge_hs).expect("isotopic hydrogen-count query");
    assert_eq!(
        isotope_count_query
            .atom(1)
            .expect("isotopic hydrogen-count atom")
            .isotope(),
        Some(13)
    );
}

#[test]
fn q21_hydrogens_found_through_or_queries_are_unmergeable() {
    let merge_hs = SmartsParseParams {
        merge_hs: true,
        ..SmartsParseParams::default()
    };

    let disjunction = parse_smarts("[C][#1,#17]", &merge_hs).expect("hydrogen OR query");
    assert_eq!(disjunction.num_atoms(), 2);
    assert!(matches!(
        disjunction.atom(1).expect("hydrogen OR atom").predicate(),
        QueryNode::Or(_)
    ));

    let nested_disjunction =
        parse_smarts("[C][#6;#1,#17]", &merge_hs).expect("nested hydrogen OR query");
    assert_eq!(nested_disjunction.num_atoms(), 2);
    assert!(matches!(
        nested_disjunction
            .atom(1)
            .expect("nested hydrogen OR atom")
            .predicate(),
        QueryNode::And(children)
            if children.iter().any(|child| matches!(child, QueryNode::Or(_)))
    ));
}

#[test]
fn q21_recursive_smarts_hydrogen_predicates_are_not_atom_hydrogens() {
    let merge_hs = SmartsParseParams {
        merge_hs: true,
        ..SmartsParseParams::default()
    };
    let graph = parse_smarts("[C][$([#1])]", &merge_hs)
        .unwrap_or_else(|error| panic!("pinned recursive SMARTS: {error}"));
    assert_eq!(graph.num_atoms(), 2);
    assert!(matches!(
        graph.atom(1).expect("recursive query atom").predicate(),
        QueryNode::Predicate(AtomQueryPredicate::RecursiveSmarts(_))
    ));
}

#[test]
fn q22_nonrecursive_merge_counts_mapped_hydrogen_and_retains_isotope_by_default() {
    let merge_hs = SmartsParseParams {
        merge_hs: true,
        ..SmartsParseParams::default()
    };
    let graph = parse_smarts("[C]([H])([H:0])([2H])", &merge_hs)
        .expect("source-default query-H merge keeps isotope filtering");

    assert_eq!(graph.num_atoms(), 2);
    assert_eq!(
        graph.atom(0).expect("merged carbon").predicate(),
        &QueryNode::And(vec![
            QueryNode::predicate(AtomQueryPredicate::AtomType {
                atomic_number: 6,
                aromatic: false,
            }),
            QueryNode::Not(Box::new(QueryNode::Predicate(
                AtomQueryPredicate::HydrogenCount(0),
            ))),
            QueryNode::Not(Box::new(QueryNode::Predicate(
                AtomQueryPredicate::HydrogenCount(1),
            ))),
        ])
    );
    assert_eq!(
        graph
            .atom(1)
            .expect("isotopic hydrogen retained by default")
            .isotope(),
        Some(2)
    );
}

#[test]
fn q22_merge_hydrogen_neighbor_count_uses_source_unsigned_width() {
    let merge_hs = SmartsParseParams {
        merge_hs: true,
        ..SmartsParseParams::default()
    };
    let smarts = format!("[C]{}", "([H])".repeat(256));
    let graph = parse_smarts(&smarts, &merge_hs).expect("merge all 256 neighboring query H atoms");

    assert_eq!(graph.num_atoms(), 1);
    let QueryNode::And(children) = graph.atom(0).expect("merged carbon").predicate() else {
        panic!("source merge should add one H-count predicate per removed atom");
    };
    assert_eq!(children.len(), 257);
    for (hydrogen_count, child) in children.iter().skip(1).enumerate() {
        assert_eq!(
            child,
            &QueryNode::Not(Box::new(QueryNode::Predicate(
                AtomQueryPredicate::HydrogenCount(hydrogen_count as i32),
            )))
        );
    }
}

#[test]
fn q23_recursive_query_hydrogen_merge_descends_without_aliasing_or_losing_serials() {
    let smarts = "[$([C]([H])[$([N][H])_8])_7]";
    let source = parse_source_case(smarts);
    let source_before = source.clone();
    let QueryNode::Predicate(AtomQueryPredicate::RecursiveSmarts(source_outer)) = source
        .atom(0)
        .expect("outer recursive SMARTS atom")
        .predicate()
    else {
        panic!("outer atom should be a recursive SMARTS predicate");
    };
    assert_eq!(source_outer.serial_number(), 7);
    assert_eq!(
        source_outer.source_smarts(),
        Some("$([C]([H])[$([N][H])_8])")
    );
    assert_eq!(
        source_outer
            .query_graph()
            .expect("source outer recursive graph")
            .num_atoms(),
        3
    );

    // RecursiveStructureQuery::clone owns a detached nested graph; changing
    // the cloned query graph must not mutate the source query or its IDs.
    let mut cloned_outer = source_outer.clone();
    let cloned_graph = cloned_outer
        .query_graph_mut()
        .expect("cloned outer recursive graph");
    cloned_graph.set_prop("_q23_clone_probe", "detached");
    assert_eq!(cloned_graph.prop("_q23_clone_probe"), Some("detached"));
    assert_eq!(
        source_outer.query_graph().unwrap().prop("_q23_clone_probe"),
        None
    );
    assert_eq!(cloned_outer.serial_number(), source_outer.serial_number());
    assert_eq!(cloned_outer.source_smarts(), source_outer.source_smarts());
    assert_eq!(source, source_before);

    let merge_hs = SmartsParseParams {
        merge_hs: true,
        ..SmartsParseParams::default()
    };
    let merged = parse_smarts(smarts, &merge_hs).expect("merge query H through nested graphs");
    assert_eq!(source, source_before);
    assert_eq!(merged.num_atoms(), 1);
    let QueryNode::Predicate(AtomQueryPredicate::RecursiveSmarts(outer)) = merged
        .atom(0)
        .expect("merged outer recursive SMARTS atom")
        .predicate()
    else {
        panic!("merged outer atom should retain its recursive predicate");
    };
    assert_eq!(outer.serial_number(), 7);
    assert_eq!(outer.source_smarts(), source_outer.source_smarts());

    let inner = outer.query_graph().expect("merged outer recursive graph");
    assert_eq!(inner.num_atoms(), 2);
    assert_eq!(inner.num_bonds(), 1);
    assert_eq!(inner.atoms()[0].id(), AtomId::new(0));
    assert_eq!(inner.atoms()[1].id(), AtomId::new(1));
    assert_eq!(inner.bonds()[0].id(), BondId::new(0));
    assert_eq!(inner.bonds()[0].endpoints(), (0, 1));
    assert_eq!(
        inner.atoms()[0].predicate(),
        &QueryNode::And(vec![
            QueryNode::predicate(AtomQueryPredicate::AtomType {
                atomic_number: 6,
                aromatic: false,
            }),
            QueryNode::Not(Box::new(QueryNode::Predicate(
                AtomQueryPredicate::HydrogenCount(0),
            ))),
        ])
    );

    let QueryNode::Predicate(AtomQueryPredicate::RecursiveSmarts(nested)) =
        inner.atoms()[1].predicate()
    else {
        panic!("nested recursive atom should retain its recursive predicate");
    };
    assert_eq!(nested.serial_number(), 8);
    assert_eq!(nested.source_smarts(), Some("$([N][H])"));
    let deepest = nested
        .query_graph()
        .expect("merged deepest recursive graph");
    assert_eq!(deepest.num_atoms(), 1);
    assert_eq!(deepest.num_bonds(), 0);
    assert_eq!(deepest.atoms()[0].id(), AtomId::new(0));
    assert_eq!(
        deepest.atoms()[0].predicate(),
        &QueryNode::And(vec![
            QueryNode::predicate(AtomQueryPredicate::AtomType {
                atomic_number: 7,
                aromatic: false,
            }),
            QueryNode::Not(Box::new(QueryNode::Predicate(
                AtomQueryPredicate::HydrogenCount(0),
            ))),
        ])
    );
}

#[test]
fn q24_direction_pairs_missing_neighbors_and_query_carrier_provenance() {
    for (smarts, expected_stereo, expected_end_direction) in [
        ("C/C=C/C", BondStereo::Trans, BondDirection::EndUpRight),
        ("C/C=C\\C", BondStereo::Cis, BondDirection::EndDownRight),
    ] {
        let graph = parse_source_case(smarts);
        assert_eq!((graph.num_atoms(), graph.num_bonds()), (4, 3), "{smarts}");
        let double_bond = graph.bond(1).expect("central double bond");
        assert_eq!(double_bond.bond().order(), BondOrder::Double, "{smarts}");
        assert_eq!(double_bond.bond().stereo(), expected_stereo, "{smarts}");
        assert_eq!(
            double_bond.bond().stereo_atoms(),
            Some([AtomId::new(0), AtomId::new(3)]),
            "{smarts}"
        );
        assert_eq!(
            double_bond.predicate(),
            &QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Double)),
            "{smarts}"
        );
        assert!(!double_bond.predicate_is_carrier_derived(), "{smarts}");
        assert_eq!(
            graph.bond(0).unwrap().bond().direction(),
            BondDirection::EndUpRight,
            "{smarts}"
        );
        assert_eq!(
            graph.bond(2).unwrap().bond().direction(),
            expected_end_direction,
            "{smarts}"
        );
    }

    for (smarts, double_bond_index, directed_bond_index) in [("C=C/C", 0, 1), ("C/C=C", 1, 0)] {
        let graph = parse_source_case(smarts);
        let double_bond = graph.bond(double_bond_index).expect("double bond");
        assert_eq!(double_bond.bond().stereo(), BondStereo::None, "{smarts}");
        assert_eq!(double_bond.bond().stereo_atoms(), None, "{smarts}");
        assert_eq!(
            graph
                .bond(directed_bond_index)
                .expect("one directed neighbor")
                .bond()
                .direction(),
            BondDirection::EndUpRight,
            "{smarts}"
        );
    }

    // MolFromSmarts invokes setBondStereoFromDirections even when CX supplied
    // stereo and no slash direction exists; the finalizer clears its marker
    // while retaining the query bond predicate and CX stereo carrier.
    let cx_stereo = parse_source_case("FC=CF |c:1|");
    let cx_double_bond = cx_stereo.bond(1).expect("CX double bond");
    assert_eq!(cx_double_bond.bond().stereo(), BondStereo::Cis);
    assert_eq!(
        cx_double_bond.bond().stereo_atoms(),
        Some([AtomId::new(0), AtomId::new(3)])
    );
    assert_eq!(
        cx_double_bond.predicate(),
        &QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Double))
    );
    assert!(!cx_double_bond.predicate_is_carrier_derived());
    assert_eq!(cx_stereo.prop("_needsDetectBondStereo"), None);
}

#[test]
fn q25_chiral_permutation_boundaries_are_query_native_and_source_ordered() {
    let carbon_query = QueryNode::predicate(AtomQueryPredicate::AtomType {
        atomic_number: 6,
        aromatic: false,
    });
    let valid = [
        ("[C@TH]", ChiralTag::TetrahedralCcw, None),
        ("[C@TH1]", ChiralTag::TetrahedralCcw, None),
        ("[C@TH2]", ChiralTag::TetrahedralCw, None),
        ("[C@AL]", ChiralTag::Allene, Some(0)),
        ("[C@AL2]", ChiralTag::Allene, Some(2)),
        ("[C@SP]", ChiralTag::SquarePlanar, Some(0)),
        ("[C@SP3]", ChiralTag::SquarePlanar, Some(3)),
        ("[C@TB]", ChiralTag::TrigonalBipyramidal, Some(0)),
        ("[C@TB20]", ChiralTag::TrigonalBipyramidal, Some(20)),
        ("[C@OH]", ChiralTag::Octahedral, Some(0)),
        ("[C@OH30]", ChiralTag::Octahedral, Some(30)),
    ];
    for (smarts, expected_tag, expected_permutation) in valid {
        let graph = parse_source_case(smarts);
        assert_eq!(graph.num_atoms(), 1, "{smarts}");
        let atom = graph.atom(0).expect("query atom");
        assert_eq!(
            atom.identity(),
            QueryAtomIdentity::Element(Element::C),
            "{smarts}"
        );
        assert_eq!(atom.predicate(), &carbon_query, "{smarts}");
        assert!(!atom.predicate_is_carrier_derived(), "{smarts}");
        assert_eq!(atom.chiral_tag(), expected_tag, "{smarts}");
        assert_eq!(atom.chiral_permutation(), expected_permutation, "{smarts}");
    }

    for (smarts, expected_error) in [
        (
            "[C@TH3]",
            "invalid chiral permutation 3 for CHI_TETRAHEDRAL",
        ),
        ("[C@AL3]", "invalid chiral permutation 3 for CHI_ALLENE"),
        (
            "[C@SP4]",
            "invalid chiral permutation 4 for CHI_SQUAREPLANAR",
        ),
        (
            "[C@TB21]",
            "invalid chiral permutation 21 for CHI_TRIGONALBIPYRAMIDAL",
        ),
        (
            "[C@OH31]",
            "invalid chiral permutation 31 for CHI_OCTAHEDRAL",
        ),
    ] {
        let error = parse_smarts(smarts, &SmartsParseParams::default())
            .expect_err("permutation above the pinned limit must fail");
        assert!(
            error.to_string().contains(expected_error),
            "{smarts}: {error}"
        );
    }

    let explicit_zero = parse_smarts("[C@SP0]", &SmartsParseParams::default())
        .expect_err("an explicit zero is rejected by the SMARTS grammar");
    assert!(
        explicit_zero
            .to_string()
            .contains("chiral permutation cannot be zero"),
        "{explicit_zero}"
    );

    for (smarts, expected_error) in [
        (
            "[C@SP4][C@TB21]",
            "invalid chiral permutation 4 for CHI_SQUAREPLANAR",
        ),
        (
            "[C@TB21][C@SP4]",
            "invalid chiral permutation 21 for CHI_TRIGONALBIPYRAMIDAL",
        ),
    ] {
        let error = parse_smarts(smarts, &SmartsParseParams::default())
            .expect_err("the first invalid source atom determines the error");
        assert!(
            error.to_string().contains(expected_error),
            "{smarts}: {error}"
        );
    }
}

#[test]
fn q10_source_atom_carriers_compose_across_primitive_reductions() {
    // smarts.yy reduces isotope+simple_atom before the following charge
    // point-query. The lexer keeps generic `a` distinct from aromatic
    // element tokens, and each NOT reduction clears only source-named carrier
    // fields before later point-query carrier transfer.
    let cases = [
        (
            "[13c+2]",
            QueryNode::and(vec![
                QueryNode::and(vec![
                    QueryNode::predicate(AtomQueryPredicate::AtomType {
                        atomic_number: 6,
                        aromatic: true,
                    }),
                    QueryNode::predicate(AtomQueryPredicate::Isotope(13)),
                ]),
                QueryNode::predicate(AtomQueryPredicate::FormalCharge(2)),
            ]),
            QueryAtomIdentity::Element(Element::C),
            true,
            Some(13),
            2,
        ),
        (
            "[as+2]",
            QueryNode::and(vec![
                QueryNode::predicate(AtomQueryPredicate::AtomType {
                    atomic_number: 33,
                    aromatic: true,
                }),
                QueryNode::predicate(AtomQueryPredicate::FormalCharge(2)),
            ]),
            QueryAtomIdentity::Element(Element::AS),
            true,
            None,
            2,
        ),
        (
            "[a+]",
            QueryNode::and(vec![
                QueryNode::predicate(AtomQueryPredicate::IsAromatic(true)),
                QueryNode::predicate(AtomQueryPredicate::FormalCharge(1)),
            ]),
            QueryAtomIdentity::Element(Element::DUMMY),
            true,
            None,
            1,
        ),
        (
            "[!13C+2]",
            QueryNode::and(vec![
                QueryNode::not(QueryNode::and(vec![
                    QueryNode::predicate(AtomQueryPredicate::AtomType {
                        atomic_number: 6,
                        aromatic: false,
                    }),
                    QueryNode::predicate(AtomQueryPredicate::Isotope(13)),
                ])),
                QueryNode::predicate(AtomQueryPredicate::FormalCharge(2)),
            ]),
            QueryAtomIdentity::Element(Element::DUMMY),
            false,
            None,
            2,
        ),
    ];

    for (smarts, expected_predicate, identity, aromatic, isotope, formal_charge) in cases {
        let graph = parse_source_case(smarts);
        let atom = graph.atom(0).expect("one composed primitive query atom");
        assert_eq!(atom.predicate(), &expected_predicate, "{smarts}");
        assert_eq!(atom.identity(), identity, "{smarts}");
        assert_eq!(atom.is_aromatic(), aromatic, "{smarts}");
        assert_eq!(atom.isotope(), isotope, "{smarts}");
        assert_eq!(atom.formal_charge(), formal_charge, "{smarts}");
        assert_eq!(atom.explicit_hydrogens(), 0, "{smarts}");
        assert!(!atom.no_implicit(), "{smarts}");
        assert!(!atom.predicate_is_carrier_derived(), "{smarts}");
    }
}

#[test]
fn q07d_neighbor_targets_keep_full_source_i32_and_distinct_factories() {
    // smarts.ll gives z and Z separate has-neighbor factories; smarts.yy
    // replaces each with its corresponding QueryOps int-count factory only
    // when a number is present. Neither query reduction writes Atom carriers.
    let defaults = [
        ("z", AtomQueryPredicate::HasHeteroatomNeighbors),
        ("Z", AtomQueryPredicate::HasAliphaticHeteroatomNeighbors),
    ];
    for (token, predicate) in defaults {
        let smarts = format!("[{token}]");
        let graph = parse_source_case(&smarts);
        let atom = graph.atom(0).expect("one heteroatom-neighbor query atom");
        assert_eq!(
            atom.predicate(),
            &QueryNode::predicate(predicate),
            "{smarts}"
        );
        assert_eq!(atom.atomic_number(), 0, "{smarts}");
        assert!(!atom.is_aromatic(), "{smarts}");
        assert_eq!(atom.formal_charge(), 0, "{smarts}");
        assert_eq!(atom.isotope(), None, "{smarts}");
        assert_eq!(atom.explicit_hydrogens(), 0, "{smarts}");
        assert!(!atom.no_implicit(), "{smarts}");
        assert!(!atom.predicate_is_carrier_derived(), "{smarts}");
    }
    assert_ne!(
        parse_source_case("[z]").atom(0).unwrap().predicate(),
        parse_source_case("[Z]").atom(0).unwrap().predicate()
    );

    let tokens: [(&str, fn(i32) -> AtomQueryPredicate); 2] = [
        ("z", AtomQueryPredicate::NumHeteroatomNeighbors),
        ("Z", AtomQueryPredicate::NumAliphaticHeteroatomNeighbors),
    ];
    let targets = [
        ("0", 0),
        ("255", 255),
        ("256", 256),
        ("2147483639", 2_147_483_639),
    ];
    for (token, make_predicate) in tokens {
        for (suffix, target) in targets {
            let smarts = format!("[{token}{suffix}]");
            let graph = parse_source_case(&smarts);
            let atom = graph.atom(0).expect("one numbered neighbor query atom");
            assert_eq!(
                atom.predicate(),
                &QueryNode::predicate(make_predicate(target)),
                "{smarts}"
            );
            assert_eq!(atom.atomic_number(), 0, "{smarts}");
            assert!(!atom.is_aromatic(), "{smarts}");
            assert_eq!(atom.formal_charge(), 0, "{smarts}");
            assert_eq!(atom.isotope(), None, "{smarts}");
            assert_eq!(atom.explicit_hydrogens(), 0, "{smarts}");
            assert!(!atom.no_implicit(), "{smarts}");
            assert!(!atom.predicate_is_carrier_derived(), "{smarts}");
        }
    }

    let raw_smarts = SmartsParseParams {
        allow_cxsmiles: false,
        parse_name: false,
        ..SmartsParseParams::default()
    };
    for token in ["z", "Z"] {
        let smarts = format!("[{token}2147483640]");
        assert_eq!(
            parse_smarts(&smarts, &raw_smarts).expect_err("source signed-int number guard"),
            SmartsParseError::InvalidAtomPrimitive {
                position: 1,
                detail: "number too large".to_owned(),
            },
            "{smarts}"
        );
    }
}

#[test]
fn q07b_charge_targets_keep_full_i32_value_and_narrow_carrier() {
    // smarts.yy charge_spec keeps the signed number as the int query target;
    // Atom::setFormalCharge separately projects to its int8_t carrier.
    let cases = [
        ("[C+]", 6, 1, 1),
        ("[C-]", 6, -1, -1),
        ("[C+0]", 6, 0, 0),
        ("[C-0]", 6, 0, 0),
        ("[C+127]", 6, 127, 127),
        ("[C+128]", 6, 128, -128),
        ("[C-128]", 6, -128, -128),
        ("[C-129]", 6, -129, 127),
        ("[C+2147483639]", 6, 2_147_483_639, -9),
        ("[C-2147483639]", 6, -2_147_483_639, 9),
        ("[H+128]", 1, 128, -128),
        ("[H-129]", 1, -129, 127),
        ("[C++]", 6, 2, 2),
        ("[C--]", 6, -2, -2),
    ];

    for (smarts, atomic_number, target, carrier) in cases {
        let graph = parse_source_case(smarts);
        let atom = graph.atom(0).expect("one charge query atom");
        let atom_identity = if atomic_number == 1 {
            AtomQueryPredicate::AtomicNumber(1)
        } else {
            AtomQueryPredicate::AtomType {
                atomic_number: 6,
                aromatic: false,
            }
        };
        assert_eq!(
            atom.predicate(),
            &QueryNode::and(vec![
                QueryNode::predicate(atom_identity),
                QueryNode::predicate(AtomQueryPredicate::FormalCharge(target)),
            ]),
            "{smarts}"
        );
        assert_eq!(atom.atomic_number(), atomic_number, "{smarts}");
        assert_eq!(atom.formal_charge(), carrier, "{smarts}");
        assert!(!atom.predicate_is_carrier_derived(), "{smarts}");
    }

    let raw_smarts = SmartsParseParams {
        allow_cxsmiles: false,
        parse_name: false,
        ..SmartsParseParams::default()
    };
    for smarts in ["[C+2147483640]", "[C-2147483640]"] {
        assert_eq!(
            parse_smarts(smarts, &raw_smarts).expect_err("source int32 number guard"),
            SmartsParseError::InvalidAtomPrimitive {
                position: 2,
                detail: "number too large".to_owned(),
            },
            "{smarts}"
        );
    }
}

#[test]
fn q07b_isotope_targets_keep_i32_value_and_u16_carrier_projection() {
    // QueryOps accepts the source int target; Atom::setIsotope stores it in
    // uint16_t. The source projects that carrier independently, with zero
    // represented as absence by the detached QueryAtom.
    for (value, expected_carrier) in [
        (0, None),
        (13, Some(13)),
        (65_535, Some(65_535)),
        (65_536, None),
        (65_537, Some(1)),
        (2_147_483_639, Some(65_527)),
    ] {
        let carbon_smarts = format!("[{value}C]");
        let carbon = parse_source_case(&carbon_smarts);
        let carbon_atom = carbon.atom(0).expect("isotopic carbon query");
        assert_eq!(
            carbon_atom.predicate(),
            &QueryNode::and(vec![
                QueryNode::predicate(AtomQueryPredicate::AtomType {
                    atomic_number: 6,
                    aromatic: false,
                }),
                QueryNode::predicate(AtomQueryPredicate::Isotope(value)),
            ]),
            "{carbon_smarts}"
        );
        assert_eq!(carbon_atom.isotope(), expected_carrier, "{carbon_smarts}");
        assert!(
            !carbon_atom.predicate_is_carrier_derived(),
            "{carbon_smarts}"
        );

        let hydrogen_smarts = format!("[{value}H]");
        let hydrogen = parse_source_case(&hydrogen_smarts);
        let hydrogen_atom = hydrogen.atom(0).expect("isotopic hydrogen query");
        assert_eq!(
            hydrogen_atom.predicate(),
            &QueryNode::and(vec![
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(1)),
                QueryNode::predicate(AtomQueryPredicate::Isotope(value)),
            ]),
            "{hydrogen_smarts}"
        );
        assert_eq!(
            hydrogen_atom.isotope(),
            expected_carrier,
            "{hydrogen_smarts}"
        );
        assert!(
            !hydrogen_atom.predicate_is_carrier_derived(),
            "{hydrogen_smarts}"
        );
    }

    for value in [0, 65_535, 65_536, 65_537, 2_147_483_639] {
        let smarts = format!("[{value}]");
        let graph = parse_source_case(&smarts);
        let atom = graph.atom(0).expect("standalone isotope query");
        assert_eq!(
            atom.predicate(),
            &QueryNode::predicate(AtomQueryPredicate::Isotope(value)),
            "{smarts}"
        );
        // The pinned standalone `number` reduction creates only a query leaf;
        // unlike number+atom and number+H, it does not call setIsotope.
        assert_eq!(atom.isotope(), None, "{smarts}");
        assert!(!atom.predicate_is_carrier_derived(), "{smarts}");
    }

    // The pinned `number HASH_TOKEN number` action and the remaining
    // `number simple_atom` tokens preserve the source query child ordering and
    // apply the isotope carrier independently.
    let source_forms = [
        (
            "[13*]",
            QueryNode::predicate(AtomQueryPredicate::Isotope(13)),
            0,
            false,
        ),
        (
            "[13a]",
            QueryNode::and(vec![
                QueryNode::predicate(AtomQueryPredicate::IsAromatic(true)),
                QueryNode::predicate(AtomQueryPredicate::Isotope(13)),
            ]),
            0,
            true,
        ),
        (
            "[13A]",
            QueryNode::and(vec![
                QueryNode::predicate(AtomQueryPredicate::IsAromatic(false)),
                QueryNode::predicate(AtomQueryPredicate::Isotope(13)),
            ]),
            0,
            false,
        ),
        (
            "[13#6]",
            QueryNode::and(vec![
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                QueryNode::predicate(AtomQueryPredicate::Isotope(13)),
            ]),
            6,
            false,
        ),
    ];
    for (smarts, expected_predicate, atomic_number, aromatic) in source_forms {
        let graph = parse_source_case(smarts);
        let atom = graph.atom(0).expect("numeric isotope query atom");
        assert_eq!(atom.predicate(), &expected_predicate, "{smarts}");
        assert_eq!(atom.atomic_number(), atomic_number, "{smarts}");
        assert_eq!(atom.is_aromatic(), aromatic, "{smarts}");
        assert_eq!(atom.isotope(), Some(13), "{smarts}");
        assert!(!atom.predicate_is_carrier_derived(), "{smarts}");
    }

    let raw_smarts = SmartsParseParams {
        allow_cxsmiles: false,
        parse_name: false,
        ..SmartsParseParams::default()
    };
    for smarts in [
        "[2147483640C]",
        "[2147483640H]",
        "[2147483640]",
        "[2147483640#6]",
    ] {
        assert!(
            matches!(
                parse_smarts(smarts, &raw_smarts).expect_err("source int32 number guard"),
                SmartsParseError::InvalidAtomPrimitive { detail, .. }
                    if detail == "number too large"
            ),
            "{smarts}"
        );
    }
}
