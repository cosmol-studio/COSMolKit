//! Fixed construction vectors from RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8.
//! Sources: SmilesParse.{h,cpp}, smarts.{ll,yy}, SmilesParseOps.cpp,
//! AddHs.cpp, Canon.cpp and CXSmilesOps.cpp. These are noncorpus regressions;
//! expectations never invoke an oracle or derive values from the search owner.

use cosmolkit::{
    BINDING_CONTRACT, BindingItem, BindingOwner, BindingTypeRole, FunctionStatus, QueryGraph,
};

fn entry(id: &str) -> &'static cosmolkit::BindingContractEntry {
    let mut rows = BINDING_CONTRACT.iter().filter(|row| row.semantic_id == id);
    let row = rows.next().unwrap_or_else(|| panic!("missing {id}"));
    assert!(rows.next().is_none(), "duplicate {id}");
    row
}

#[test]
fn canonical_query_value_metadata_is_always_registered() {
    let _: Option<QueryGraph> = None;
    let row = entry("types.QueryGraph");
    assert_eq!(row.item, BindingItem::Type);
    assert_eq!(row.owner, BindingOwner::Type);
    assert_eq!(row.rust_path.replace(' ', ""), "crate::QueryGraph");
    assert_eq!(row.python_name, "QueryGraph");
    assert_eq!(row.javascript_name, "QueryGraph");
    assert_eq!(row.feature, "metadata");
    assert_eq!(row.status, FunctionStatus::Experimental);
    assert_eq!(row.type_role, Some(BindingTypeRole::Value));
    assert_eq!(row.callable, None);
}

#[test]
fn existing_element_prefix_retains_exact_order() {
    // Accepted b372 canonical Element methods occur between the two original
    // type entries. Retain all original prefix and QueryGraph placement checks.
    let mut expected = vec![
        "types.Element",
        "Element.from_atomic_number",
        "Element.from_symbol",
        "Element.atomic_number",
        "Element.symbol",
        "types.ElementInfo",
    ];
    if cfg!(feature = "cap-valence") {
        expected.push("module.element_info");
    }
    assert_eq!(
        BINDING_CONTRACT
            .iter()
            .take(expected.len())
            .map(|row| row.semantic_id)
            .collect::<Vec<_>>(),
        expected
    );
    assert_eq!(
        BINDING_CONTRACT[expected.len()].semantic_id,
        "types.QueryGraph"
    );
}

#[cfg(not(feature = "cap-search"))]
#[test]
fn search_constructor_is_absent_without_capability() {
    for id in [
        "types.SmartsParseParams",
        "types.SmartsParseError",
        "search.parse_smarts",
        "search.parse_smarts_with_params",
    ] {
        assert!(BINDING_CONTRACT.iter().all(|row| row.semantic_id != id));
    }
}

#[cfg(feature = "cap-search")]
mod enabled {
    use super::*;
    use cosmolkit::{
        AtomId, AtomQueryPredicate as A, AtomRangeBounds, AtomRangeDataFunction, AtomRangeQuery,
        BindingDefault, BindingKind, BondId, BondOrder, BondQueryPredicate as B, BondStereo,
        ChiralTag, Element, PropertyValue, QueryAtomIdentity, QueryNode as Q, SmartsParseError,
        SmartsParseParams, StateModel, search::parse_smarts_with_params as parse_smarts,
    };

    fn fixture_text(text: &cosmolkit::PropertyText) -> &str {
        std::str::from_utf8(text.as_bytes()).expect("original parser fixture is UTF-8")
    }

    fn parse(text: &str) -> QueryGraph {
        let graph = parse_smarts(text, &SmartsParseParams::default()).unwrap();
        graph.validate().unwrap();
        graph
    }

    fn atom_type(number: u8, aromatic: bool) -> Q<A> {
        Q::Predicate(A::AtomType {
            atomic_number: number,
            aromatic,
        })
    }

    fn parameters_equal(left: &SmartsParseParams, right: &SmartsParseParams) {
        assert_eq!(left.allow_cxsmiles, right.allow_cxsmiles);
        assert_eq!(left.strict_cxsmiles, right.strict_cxsmiles);
        assert_eq!(left.parse_name, right.parse_name);
        assert_eq!(left.merge_hs, right.merge_hs);
        assert_eq!(left.skip_cleanup, right.skip_cleanup);
        assert_eq!(left.debug_parse, right.debug_parse);
        assert_eq!(left.replacements, right.replacements);
    }

    #[test]
    fn construction_types_and_callable_have_exact_contract() {
        let _: fn(&str, &SmartsParseParams) -> Result<QueryGraph, SmartsParseError> = parse_smarts;
        let _: fn(&str) -> Result<QueryGraph, SmartsParseError> = cosmolkit::parse_smarts;
        assert_eq!(
            entry("search.parse_smarts")
                .callable
                .unwrap()
                .parameters
                .len(),
            1
        );
        assert_eq!(
            cosmolkit::search::parse_smarts("C").unwrap(),
            parse_smarts("C", &SmartsParseParams::default()).unwrap()
        );
        for (name, role) in [
            ("SmartsParseParams", BindingTypeRole::Parameter),
            ("SmartsParseError", BindingTypeRole::Error),
        ] {
            let row = entry(&format!("types.{name}"));
            assert_eq!(row.item, BindingItem::Type);
            assert_eq!(row.owner, BindingOwner::Type);
            assert_eq!(row.type_role, Some(role));
            assert_eq!(row.feature, "cap-search");
            assert_eq!(row.status, FunctionStatus::Experimental);
            assert_eq!(row.rust_path.replace(' ', ""), format!("crate::{name}"));
            assert_eq!(row.python_name, name);
            assert_eq!(row.javascript_name, name);
            assert_eq!(row.callable, None);
        }
        let row = entry("search.parse_smarts_with_params");
        assert_eq!(row.item, BindingItem::Callable);
        assert_eq!(row.owner, BindingOwner::Module);
        assert_eq!(
            row.rust_path.replace(' ', ""),
            "crate::parse_smarts_with_params"
        );
        assert_eq!(row.python_name, "parse_smarts_with_params");
        assert_eq!(row.javascript_name, "parseSmartsWithParams");
        assert_eq!(row.feature, "cap-search");
        assert_eq!(row.status, FunctionStatus::Experimental);
        let callable = row.callable.unwrap();
        assert_eq!(callable.kind, BindingKind::Module);
        assert_eq!(callable.receiver, None);
        assert_eq!(callable.state_model, StateModel::ValueReturning);
        assert_eq!(callable.operation_semantic_id, None);
        assert_eq!(callable.output_type.replace(' ', ""), "crate::QueryGraph");
        assert_eq!(
            callable.error_type.unwrap().replace(' ', ""),
            "crate::SmartsParseError"
        );
        assert_eq!(callable.parameters.len(), 2);
        for (parameter, (name, ty)) in callable
            .parameters
            .iter()
            .zip([("text", "&str"), ("params", "&crate::SmartsParseParams")])
        {
            assert_eq!(parameter.name, name);
            assert_eq!(parameter.type_name.replace(' ', ""), ty);
            assert_eq!(parameter.default, BindingDefault::Required);
        }
    }

    #[test]
    fn seven_defaults_and_parameters_are_preserved() {
        let params = SmartsParseParams::default();
        assert!(params.allow_cxsmiles && params.strict_cxsmiles && params.parse_name);
        assert!(!params.merge_hs && !params.skip_cleanup && !params.debug_parse);
        assert!(params.replacements.is_empty());
        let before = params.clone();
        let text = String::from("C(N)O");
        let text_before = text.clone();
        let graph = parse_smarts(&text, &params).unwrap();
        assert_eq!((graph.num_atoms(), graph.num_bonds()), (3, 2));
        assert!(parse_smarts("C?", &params).is_err());
        parameters_equal(&params, &before);
        assert_eq!(text, text_before);
        let empty = parse_smarts("", &params).unwrap();
        assert_eq!((empty.num_atoms(), empty.num_bonds()), (0, 0));
        empty.validate().unwrap();
    }

    #[test]
    fn primitive_predicates_preserve_source_order_and_carriers() {
        for (text, number, predicate) in [
            ("C", 6, atom_type(6, false)),
            ("n", 7, atom_type(7, true)),
            ("*", 0, Q::Predicate(A::Any)),
            ("[#119]", 119, Q::Predicate(A::AtomicNumber(119))),
        ] {
            let graph = parse(text);
            let atom = graph.atom(0).unwrap();
            assert_eq!(
                atom.identity(),
                QueryAtomIdentity::from_atomic_number(number)
            );
            assert_eq!(atom.predicate(), &predicate);
            assert!(!atom.predicate_is_carrier_derived());
        }
        let graph = parse("[13CH2+:7]");
        let atom = graph.atom(0).unwrap();
        assert_eq!(
            atom.predicate(),
            &Q::And(vec![
                Q::And(vec![
                    Q::And(vec![atom_type(6, false), Q::Predicate(A::Isotope(13))]),
                    Q::Predicate(A::HydrogenCount(2))
                ]),
                Q::Predicate(A::FormalCharge(1))
            ])
        );
        assert_eq!(atom.identity(), QueryAtomIdentity::Element(Element::C));
        assert_eq!(atom.isotope(), Some(13));
        assert_eq!(atom.formal_charge(), 1);
        assert_eq!(atom.explicit_hydrogens(), 2);
        assert_eq!(atom.atom_map(), Some(7));
        assert_eq!(
            parse("[#6,#7;!H0]").atom(0).unwrap().predicate(),
            &Q::And(vec![
                Q::Or(vec![
                    Q::Predicate(A::AtomicNumber(6)),
                    Q::Predicate(A::AtomicNumber(7))
                ]),
                Q::Not(Box::new(Q::Predicate(A::HydrogenCount(0))))
            ])
        );
        assert_eq!(
            parse("[D{2-4}]").atom(0).unwrap().predicate(),
            &Q::Predicate(A::Range(AtomRangeQuery::new(
                AtomRangeBounds::Inclusive {
                    lower: 2,
                    upper: 4,
                    lower_open: false,
                    upper_open: false
                },
                AtomRangeDataFunction::ExplicitDegree
            )))
        );
    }

    #[test]
    fn branches_rings_and_dative_bonds_preserve_source_order() {
        let branch = parse("C(N)O");
        assert_eq!(
            branch
                .bonds()
                .iter()
                .map(|bond| bond.endpoints())
                .collect::<Vec<_>>(),
            [(0, 1), (0, 2)]
        );
        let ring = parse("C1CC1");
        assert_eq!(
            ring.bonds()
                .iter()
                .map(|bond| bond.endpoints())
                .collect::<Vec<_>>(),
            [(0, 1), (1, 2), (2, 0)]
        );
        let dative = parse("N<-C");
        assert_eq!(dative.bond(0).unwrap().endpoints(), (1, 0));
        assert_eq!(dative.bond(0).unwrap().bond().order(), BondOrder::Dative);
        assert_eq!(
            dative.bond(0).unwrap().predicate(),
            &Q::Predicate(B::Order(BondOrder::Dative))
        );
        for (text, predicate) in [
            ("C~C", B::Any),
            ("C:C", B::Order(BondOrder::Aromatic)),
            ("C@C", B::IsInRing(true)),
        ] {
            assert_eq!(
                parse(text).bond(0).unwrap().predicate(),
                &Q::Predicate(predicate)
            );
        }
    }

    #[test]
    fn recursive_ids_and_graphs_remain_owned() {
        let graph = parse("[$(C),$(C),$(N)_7]");
        let Q::Or(children) = graph.atom(0).unwrap().predicate() else {
            panic!("ordered recursive OR")
        };
        assert_eq!(
            children.len(),
            2,
            "source OR reduction has two ordered children"
        );
        let Q::Or(left) = &children[0] else {
            panic!("left-associated recursive OR")
        };
        assert_eq!(
            left.len(),
            2,
            "source first OR reduction has two ordered children"
        );
        let recursive = [&left[0], &left[1], &children[1]]
            .into_iter()
            .map(|child| {
                let Q::Predicate(A::RecursiveSmarts(recursive)) = child else {
                    panic!("recursive leaf")
                };
                recursive
            })
            .collect::<Vec<_>>();
        assert_eq!(
            recursive
                .iter()
                .map(|leaf| leaf.serial_number())
                .collect::<Vec<_>>(),
            [100, 100, 7]
        );
        let mut copied = recursive[0].clone();
        copied
            .query_graph_mut()
            .unwrap()
            .set_prop("copy", "independent");
        assert_eq!(recursive[0].query_graph().unwrap().prop("copy"), None);
        assert_eq!(recursive[1].query_graph().unwrap().prop("copy"), None);
        assert_eq!(
            recursive[0]
                .query_graph()
                .unwrap()
                .atom(0)
                .unwrap()
                .predicate(),
            &atom_type(6, false)
        );
    }

    #[test]
    fn replacement_fixpoint_and_source_name_trim() {
        let mut params = SmartsParseParams {
            allow_cxsmiles: false,
            ..Default::default()
        };
        params.replacements.insert("N".into(), "C".into());
        params.replacements.insert("O".into(), "N".into());
        let before = params.clone();
        let graph = parse_smarts("O-O\tOriginal O\u{00a0}", &params).unwrap();
        assert_eq!(
            graph.name().unwrap().map(fixture_text),
            Some("Original O\u{00a0}")
        );
        assert_eq!((graph.num_atoms(), graph.num_bonds()), (2, 1));
        assert!(
            graph
                .atoms()
                .iter()
                .all(|atom| atom.predicate() == &atom_type(6, false))
        );
        assert!(parse_smarts("\tC name", &params).is_err());
        parameters_equal(&params, &before);
    }

    #[test]
    fn cx_options_names_and_error_priority() {
        let text = "C |$label$| note";
        let graph = parse(text);
        assert_eq!(
            graph.atom(0).unwrap().prop("atomLabel"),
            Some(&PropertyValue::String("label".into()))
        );
        assert_eq!(
            graph.prop("_CXSMILES_Data"),
            Some(&PropertyValue::String("|$label$|".into()))
        );
        assert_eq!(graph.name().unwrap().map(fixture_text), Some("note"));
        let no_name = SmartsParseParams {
            parse_name: false,
            ..Default::default()
        };
        let graph = parse_smarts(text, &no_name).unwrap();
        assert_eq!(graph.name().unwrap().map(fixture_text), None);
        assert!(matches!(
            parse_smarts("C note", &no_name),
            Err(SmartsParseError::CxSmiles(_))
        ));
        let lenient = SmartsParseParams {
            strict_cxsmiles: false,
            ..no_name.clone()
        };
        assert_eq!(
            parse_smarts("C note", &lenient)
                .unwrap()
                .name()
                .unwrap()
                .map(fixture_text),
            None
        );
        let no_cx = SmartsParseParams {
            allow_cxsmiles: false,
            ..Default::default()
        };
        let graph = parse_smarts(text, &no_cx).unwrap();
        assert_eq!(
            graph.name().unwrap().map(fixture_text),
            Some("|$label$| note")
        );
        assert_eq!(graph.atom(0).unwrap().prop("atomLabel"), None);
        assert!(matches!(
            parse_smarts("C |sense|", &Default::default()),
            Err(SmartsParseError::CxSmiles(_))
        ));
        assert!(matches!(
            parse_smarts("C- |sense|", &Default::default()),
            Err(SmartsParseError::UnexpectedEnd(_))
        ));
    }

    #[test]
    fn lenient_cx_keeps_commits_and_source_cursor() {
        let params = SmartsParseParams {
            strict_cxsmiles: false,
            ..Default::default()
        };
        for (text, prefix, label) in [
            ("C |sense| ignored", "|", None),
            ("C |$label$ s:0:x| ignored", "|$label$ s:0:", Some("label")),
        ] {
            let graph = parse_smarts(text, &params).unwrap();
            assert_eq!(
                graph.prop("_CXSMILES_Data"),
                Some(&PropertyValue::String(prefix.into()))
            );
            assert_eq!(graph.name().unwrap().map(fixture_text), None);
            assert_eq!(
                graph.atom(0).unwrap().prop("atomLabel"),
                label
                    .map(|value| PropertyValue::String(value.into()))
                    .as_ref()
            );
        }
        let graph = parse_smarts("C-C |$label$ wU:0.0,1.0| ignored", &params).unwrap();
        assert_eq!(
            graph.prop("_CXSMILES_Data"),
            Some(&PropertyValue::String("|$label$ wU:0.0,1.0".into()))
        );
        assert_eq!(
            graph.bond(0).unwrap().bond().direction(),
            cosmolkit::BondDirection::BeginWedge
        );
        assert_eq!(graph.name().unwrap().map(fixture_text), None);
    }

    #[test]
    fn merge_hs_preserves_query_counts_and_compacted_rows() {
        let params = SmartsParseParams {
            merge_hs: true,
            ..Default::default()
        };
        let graph = parse_smarts("[C]([H])([H:0])([2H])", &params).unwrap();
        assert_eq!((graph.num_atoms(), graph.num_bonds()), (2, 1));
        assert_eq!(
            graph.atom(0).unwrap().predicate(),
            &Q::And(vec![
                atom_type(6, false),
                Q::Not(Box::new(Q::Predicate(A::HydrogenCount(0)))),
                Q::Not(Box::new(Q::Predicate(A::HydrogenCount(1))))
            ])
        );
        assert_eq!(graph.atom(1).unwrap().isotope(), Some(2));
        assert_eq!(graph.atom(1).unwrap().id(), AtomId::new(1));
        assert_eq!(graph.bond(0).unwrap().id(), BondId::new(0));
        assert_eq!(graph.bond(0).unwrap().endpoints(), (0, 1));
        graph.validate().unwrap();
        assert_eq!(parse_smarts("[C][#1,#17]", &params).unwrap().num_atoms(), 2);
        let recursive = parse_smarts("[$([C][H])_7]", &params).unwrap();
        let Q::Predicate(A::RecursiveSmarts(leaf)) = recursive.atom(0).unwrap().predicate() else {
            panic!("recursive query")
        };
        assert_eq!(leaf.serial_number(), 7);
        assert_eq!(leaf.query_graph().unwrap().num_atoms(), 1);
    }

    #[test]
    fn bond_stereo_and_chiral_errors_are_typed() {
        for (text, stereo) in [
            ("C/C=C/C", BondStereo::Trans),
            ("C/C=C\\C", BondStereo::Cis),
        ] {
            let graph = parse(text);
            let bond = graph.bond(1).unwrap().bond();
            assert_eq!(bond.stereo(), stereo);
            assert_eq!(bond.stereo_atoms(), Some([AtomId::new(0), AtomId::new(3)]));
        }
        assert!(
            matches!(parse_smarts("[C@SP4]", &Default::default()), Err(SmartsParseError::Parse(message)) if message.contains("invalid chiral permutation 4"))
        );
        assert!(
            matches!(parse_smarts("[C@SP0]", &Default::default()), Err(SmartsParseError::InvalidAtomPrimitive { detail, .. }) if detail == "chiral permutation cannot be zero")
        );
    }

    #[test]
    fn cleanup_and_skip_cleanup_preserve_source_parser_state() {
        // smarts.yy initializes source component starts; getUnspecifiedQueryBond
        // stores the source marker. CleanupAfterParsing clears both by default.
        let params = SmartsParseParams {
            skip_cleanup: true,
            ..Default::default()
        };
        let graph = parse_smarts("CC.O", &params).unwrap();
        assert_eq!(
            graph.atom(0).unwrap().prop("_SmilesStart"),
            Some(&PropertyValue::Int(1))
        );
        assert_eq!(
            graph.atom(2).unwrap().prop("_SmilesStart"),
            Some(&PropertyValue::Int(1))
        );
        assert_eq!(
            graph.bond(0).unwrap().bond().prop("_unspecifiedOrder"),
            Some(&PropertyValue::Int(1))
        );
        assert!(
            graph
                .bond(0)
                .unwrap()
                .bond()
                .prop("_cxsmilesBondIdx")
                .is_some()
        );
        let cleaned = parse("CC.O");
        assert!(
            cleaned
                .atoms()
                .iter()
                .all(|atom| atom.prop("_SmilesStart").is_none())
        );
        assert!(
            cleaned
                .bonds()
                .iter()
                .all(|bond| bond.bond().prop("_unspecifiedOrder").is_none()
                    && bond.bond().prop("_cxsmilesBondIdx").is_none())
        );
    }

    #[test]
    fn attachment_cleanup_follows_source_ap_branches() {
        // SmilesParseOps.cpp CleanupAfterParsing: only atomic number 0 and
        // exactly _AP1/_AP2 receive the typed attachment-point integer.
        for (text, value) in [
            ("* |$_AP1$|", Some(1)),
            ("* |$_AP2$|", Some(2)),
            ("C |$_AP1$|", None),
            ("* |$_AP3$|", None),
        ] {
            assert_eq!(
                parse(text).atom(0).unwrap().prop("_fromAttchpt"),
                value.map(PropertyValue::Int).as_ref(),
                "{text}"
            );
        }
        // types.h spells the common property value _fromAttchpt.
        // Retain the old literal as an explicit absence check.
        for text in ["* |$_AP1$|", "* |$_AP2$|"] {
            assert_eq!(parse(text).atom(0).unwrap().prop("_fromAttachPoint"), None);
        }
        let params = SmartsParseParams {
            skip_cleanup: true,
            ..Default::default()
        };
        assert_eq!(
            parse_smarts("* |$_AP1$|", &params)
                .unwrap()
                .atom(0)
                .unwrap()
                .prop("_fromAttchpt"),
            None
        );
    }

    #[test]
    fn cx_ring_bond_indices_follow_source_parse_order() {
        // smarts.yy allocates closing-ring parse index 2 before the oxygen
        // bond. CloseMolRings appends the ring at final bond row 3.
        let graph = parse("C1CC1O |Z:2|");
        assert_eq!(graph.bond(2).unwrap().endpoints(), (2, 3));
        assert_eq!(graph.bond(2).unwrap().bond().order(), BondOrder::Single);
        assert_eq!(graph.bond(3).unwrap().endpoints(), (2, 0));
        assert_eq!(graph.bond(3).unwrap().bond().order(), BondOrder::Zero);
    }

    #[test]
    fn atom_chirality_is_finalized_in_source_neighbor_order() {
        // AdjustAtomChiralityFlags / Canon::chiralAtomNeedsTagInversion: initial
        // explicit H and ring-order permutations contribute distinct swaps.
        assert_eq!(
            parse("[C@](Cl)(F)C").atom(0).unwrap().chiral_tag(),
            ChiralTag::TetrahedralCcw
        );
        assert_eq!(
            parse("[C@H](Cl)(F)C").atom(0).unwrap().chiral_tag(),
            ChiralTag::TetrahedralCw
        );
        assert_eq!(
            parse("F[C@]1(Br)I.Cl1").atom(1).unwrap().chiral_tag(),
            // GetBondOrdering = [0,3,1,2], storage = [0,1,2,3]: two swaps.
            // Degree four prevents the additional inversion. Keep original
            // input/assertion and correct its source-inconsistent expected tag.
            ChiralTag::TetrahedralCcw
        );
        assert_eq!(
            parse("F[C@H]1CCCC1").atom(1).unwrap().chiral_tag(),
            ChiralTag::TetrahedralCw
        );
    }

    #[test]
    fn diagnostic_gap_and_syntax_errors_keep_typed_categories() {
        let params = SmartsParseParams {
            debug_parse: true,
            ..Default::default()
        };
        let before = params.clone();
        for text in ["C", "C?"] {
            assert_eq!(
                parse_smarts(text, &params),
                Err(SmartsParseError::UnsupportedFeature(
                    "Bison debug_parse diagnostic output"
                ))
            );
        }
        parameters_equal(&params, &before);
        assert_eq!(
            parse_smarts("[C", &Default::default()),
            Err(SmartsParseError::UnclosedBracket(0))
        );
        assert_eq!(
            parse_smarts("C?", &Default::default()),
            Err(SmartsParseError::UnexpectedCharacter {
                position: 2,
                character: b'?',
                context: "unexpected character in SMARTS string".to_owned()
            })
        );
        assert_eq!(
            parse_smarts("C1CC", &Default::default()),
            Err(SmartsParseError::Parse("unclosed ring".to_owned()))
        );
        assert!(matches!(
            parse_smarts("[;C]", &Default::default()),
            Err(SmartsParseError::InvalidAtomPrimitive { .. })
        ));
        assert!(matches!(
            parse_smarts("[$()]", &Default::default()),
            Err(SmartsParseError::InvalidAtomPrimitive { .. })
        ));
        assert!(matches!(
            parse_smarts("C(", &Default::default()),
            Err(SmartsParseError::UnexpectedCharacter { .. })
        ));
    }
}
