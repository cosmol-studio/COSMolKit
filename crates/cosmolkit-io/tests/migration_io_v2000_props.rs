use cosmolkit_io::{
    MolBlockReadParams, MolBlockRecord, QueryMolBlockRecord, SdfReadError, read_mol_block_detached,
    read_mol_block_detached_with_params,
};
use cosmolkit_model::{
    AtomId, AtomQueryPredicate, BondId, QueryNode, SGroupBracketStyle, SGroupConnection,
    SubstanceGroupKind,
};
use cosmolkit_types::{BondOrder, Element};

fn atom_line(symbol: &str, charge_code: i32) -> String {
    format!(
        "{:>10.4}{:>10.4}{:>10.4} {symbol:<3}{:>2}{charge_code:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}{:>3}",
        0.0, 0.0, 0.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    )
}

fn bond_line(begin: usize, end: usize, order: usize) -> String {
    format!(
        "{begin:>3}{end:>3}{order:>3}{:>3}{:>3}{:>3}{:>3}",
        0, 0, 0, 0
    )
}

fn molblock(atoms: &[String], bonds: &[String], properties: &str) -> String {
    format!(
        "props\n  pinned          2D\ncomment\n{:>3}{:>3}  0  0  0  0  0  0  0  0999 V2000\n{}{}{}",
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
        properties,
    )
}

fn concrete(
    block: &str,
) -> (
    cosmolkit_model::TopologyBlock,
    cosmolkit_model::MoleculeProperties,
) {
    match read_mol_block_detached(block).expect("concrete V2000 property record") {
        MolBlockRecord::Concrete {
            topology,
            properties,
            ..
        } => (topology, properties),
        MolBlockRecord::Query(_) => panic!("expected concrete record"),
    }
}

fn query(block: &str) -> QueryMolBlockRecord {
    match read_mol_block_detached(block).expect("query V2000 property record") {
        MolBlockRecord::Query(record) => record,
        MolBlockRecord::Concrete { .. } => panic!("expected query record"),
    }
}

fn query_contains(
    node: &QueryNode<AtomQueryPredicate>,
    predicate: &impl Fn(&AtomQueryPredicate) -> bool,
) -> bool {
    match node {
        QueryNode::Predicate(value) => predicate(value),
        QueryNode::And(children) | QueryNode::Or(children) | QueryNode::Xor(children) => children
            .iter()
            .any(|child| query_contains(child, predicate)),
        QueryNode::Not(child) => query_contains(child, predicate),
    }
}

#[test]
fn charge_radical_and_isotope_records_preserve_reset_multiline_and_signed_rules() {
    let input = molblock(
        &[atom_line("C", 1), atom_line("N", 2), atom_line("O", 3)],
        &[],
        concat!(
            "M  CHG  1   1  -1\n",
            "M  CHG  1   2   2\n",
            "M  RAD  3   1   0   2   1   3   2\n",
            "M  ISO  3   1  13   2  -1   3    \n",
            "M  END\n",
        ),
    );
    let (topology, _) = concrete(&input);
    assert_eq!(
        topology
            .atoms
            .iter()
            .map(|atom| atom.formal_charge())
            .collect::<Vec<_>>(),
        [-1, 2, 0]
    );
    assert_eq!(
        topology
            .atoms
            .iter()
            .map(|atom| atom.radical_electrons())
            .collect::<Vec<_>>(),
        [0, 2, 1]
    );
    assert_eq!(topology.atoms[0].isotope(), Some(13));
    assert_eq!(topology.atoms[1].isotope(), None);
    assert_eq!(topology.atoms[2].isotope(), None);
}

#[test]
fn property_numeric_boundaries_and_bookmark_error_order_are_structured() {
    let atom = atom_line("C", 0);
    let iso_bad_atom = molblock(&[atom.clone()], &[], "M  ISO  1   2    \nM  END\n");
    assert!(
        read_mol_block_detached(&iso_bad_atom)
            .unwrap_err()
            .to_string()
            .contains("Atom index 2 out of range")
    );

    let apo_bad_atom_and_value = molblock(&[atom.clone()], &[], "M  APO  1   2   9\nM  END\n");
    assert!(
        read_mol_block_detached(&apo_bad_atom_and_value)
            .unwrap_err()
            .to_string()
            .contains("Atom index 2 out of range")
    );

    assert_eq!(
        read_mol_block_detached(&molblock(
            &[atom.clone()],
            &[],
            "M  CHG  1   1 128\nM  END\n"
        )),
        Err(SdfReadError::Unsupported(
            "V2000 CHG charge outside the detached i8 charge model"
        ))
    );
    let (topology, _) = concrete(&molblock(&[atom], &[], "M  ISO  1   1 999\nM  END\n"));
    assert_eq!(topology.atoms[0].isotope(), Some(999));
}

#[test]
fn new_and_old_atom_lists_preserve_source_element_negation_and_constraints() {
    let mut constrained_atom = atom_line("C", 0);
    constrained_atom.replace_range(42..45, "  2");
    let one_member = query(&molblock(
        &[constrained_atom],
        &[],
        "M  ALS   1  1 T N   \nM  END\n",
    ));
    let atom = &one_member.query.atoms()[0];
    assert_eq!(atom.atom().element(), Element::N);
    assert_eq!(atom.atom().prop("_MolFileAtomQuery"), Some("1"));
    assert!(query_contains(atom.predicate(), &|predicate| {
        predicate == &AtomQueryPredicate::AtomicNumberNotIn(vec![7])
    }));
    assert!(query_contains(atom.predicate(), &|predicate| {
        predicate == &AtomQueryPredicate::ImplicitHydrogenCountLessEqual(1)
    }));

    let multi_member = query(&molblock(
        &[atom_line("C", 0)],
        &[],
        "M  ALS   1  2 F N   O   \nM  END\n",
    ));
    assert_eq!(
        multi_member.query.atoms()[0].atom().element(),
        Element::DUMMY
    );
    assert!(query_contains(
        multi_member.query.atoms()[0].predicate(),
        &|predicate| predicate == &AtomQueryPredicate::AtomicNumberIn(vec![7, 8])
    ));

    let old = query(&molblock(
        &[atom_line("C", 0)],
        &[],
        "  1 T    2   6   8\nM  END\n",
    ));
    assert_eq!(old.query.atoms()[0].atom().element(), Element::C);
    assert!(query_contains(
        old.query.atoms()[0].predicate(),
        &|predicate| { predicate == &AtomQueryPredicate::AtomicNumberNotIn(vec![6, 8]) }
    ));
}

#[test]
fn substitution_unsaturation_and_ring_count_cover_every_source_branch() {
    let atoms = (0..8).map(|_| atom_line("C", 0)).collect::<Vec<_>>();
    let bonds = vec![bond_line(1, 2, 1), bond_line(2, 3, 1)];
    let record = query(&molblock(
        &atoms,
        &bonds,
        concat!(
            "M  SUB  4   1  -1   2  -2   3   6   4   0\n",
            "M  UNS  2   4   0   5   1\n",
            "M  RBC  5   4  -1   5  -2   6   1   7   3   8   4\n",
            "M  END\n",
        ),
    ));
    let predicates = record
        .query
        .atoms()
        .iter()
        .map(|atom| atom.predicate())
        .collect::<Vec<_>>();
    for (index, expected) in [(0, 0), (1, 2), (2, 6)] {
        assert!(query_contains(predicates[index], &|predicate| {
            predicate == &AtomQueryPredicate::ExplicitDegree(expected)
        }));
    }
    assert!(query_contains(predicates[4], &|predicate| {
        predicate == &AtomQueryPredicate::IsUnsaturated
    }));
    for (index, expected) in [(3, 0), (4, 0xDEAD_BEEF), (5, 1), (6, 3)] {
        assert!(query_contains(predicates[index], &|predicate| {
            predicate == &AtomQueryPredicate::RingBondCount(expected)
        }));
    }
    assert!(query_contains(predicates[7], &|predicate| {
        predicate == &AtomQueryPredicate::RingBondCountLessEqual(4)
    }));

    for line in [
        "M  SUB  1   1   7",
        "M  UNS  1   1   2",
        "M  RBC  1   1   5",
        "M  SUB  1   9   0",
    ] {
        let bad = molblock(&[atom_line("C", 0)], &[], &format!("{line}\nM  END\n"));
        assert!(read_mol_block_detached(&bad).is_err(), "{line}");
    }
}

#[test]
fn rgroup_and_marvin_smarts_records_create_canonical_typed_queries() {
    let rgroup = query(&molblock(
        &[atom_line("*", 0)],
        &[],
        "M  RGP  1   1   7\nM  END\n",
    ));
    let atom = rgroup.query.atoms()[0].atom();
    assert_eq!(atom.prop("_MolFileRLabel"), Some("7"));
    assert_eq!(atom.prop("dummyLabel"), Some("R7"));
    assert_eq!(atom.isotope(), Some(7));

    let marvin = query(&molblock(
        &[atom_line("C", 0)],
        &[],
        "M  MRV SMA   1 [#6,#7]\nM  END\n",
    ));
    let atom = &marvin.query.atoms()[0];
    assert_eq!(atom.atom().prop("MRV SMA"), Some("[#6,#7]"));
    assert_eq!(atom.atom().prop("_MolFileAtomQuery"), Some("1"));
    assert!(query_contains(atom.predicate(), &|predicate| {
        matches!(
            predicate,
            AtomQueryPredicate::RecursiveSmarts(recursive)
                if recursive.source_smarts() == Some("[#6,#7]")
                    && recursive.query_graph().is_some()
        )
    }));

    let ignored = concrete(&molblock(
        &[atom_line("C", 0)],
        &[],
        "M  MRV FOO   1 impossible\nM  END\n",
    ));
    assert_eq!(ignored.0.atoms[0].prop("MRV SMA"), None);

    let invalid = molblock(&[atom_line("C", 0)], &[], "M  MRV SMA   1 [\nM  END\n");
    assert!(
        read_mol_block_detached(&invalid)
            .unwrap_err()
            .to_string()
            .contains("Cannot parse smarts")
    );
}

#[test]
fn aliases_extended_records_skip_and_termination_keep_exact_state() {
    let input = molblock(
        &[atom_line("C", 0), atom_line("N", 0)],
        &[bond_line(1, 2, 1)],
        concat!(
            "A    1\ncarbon alias\n",
            "V    2 atom value\n",
            "G  1\ndeprecated\n",
            "M  PXA   1 pxa payload\n",
            "M  APO  2   1   3   2   1\n",
            "M  LIN  1   1   3   2   0\n",
            "S  SKP  1\n",
            "M  CHG  1   1  -1\n",
            "M  END\n",
        ),
    );
    let (topology, properties) = concrete(&input);
    assert_eq!(topology.atoms[0].prop("molFileAlias"), Some("carbon alias"));
    assert_eq!(topology.atoms[1].prop("molFileValue"), Some("atom value"));
    assert_eq!(topology.atoms[0].prop("_MolFile_PXA"), Some(" pxa payload"));
    assert_eq!(topology.atoms[0].prop("molAttachPoint"), Some("-1"));
    assert_eq!(topology.atoms[1].prop("molAttachPoint"), Some("1"));
    assert_eq!(topology.atoms[0].formal_charge(), 0);
    assert_eq!(properties.prop("_MolFileLinkNodes"), Some("1 3 1 1 2"));

    let missing_end = input.replace("M  END\n", "");
    assert!(
        read_mol_block_detached(&missing_end)
            .unwrap_err()
            .to_string()
            .contains("M  END missing")
    );
    let negative_skip = molblock(&[], &[], "S  SKP -1\nM  END\n");
    assert!(
        read_mol_block_detached(&negative_skip)
            .unwrap_err()
            .to_string()
            .contains("negative skip value")
    );
}

#[test]
fn zero_bond_charge_hydrogen_and_attachment_records_preserve_typed_effects() {
    let input = molblock(
        &[atom_line("B", 0), atom_line("N", 0)],
        &[bond_line(1, 2, 1)],
        concat!(
            "M  ZBO  1   1   0\n",
            "M  ZCH  2   1   1   2  -1\n",
            "M  HYD  2   1   3   2  -1\n",
            "M  APO  1   2   0\n",
            "M  END\n",
        ),
    );
    let (topology, _) = concrete(&input);
    assert_eq!(topology.bonds[0].order(), BondOrder::Zero);
    assert_eq!(topology.atoms[0].formal_charge(), 1);
    assert_eq!(topology.atoms[1].formal_charge(), -1);
    assert_eq!(topology.atoms[0].explicit_hydrogens(), 3);
    assert_eq!(topology.atoms[0].prop("_ZBO_H"), Some("1"));
    assert_eq!(topology.atoms[1].explicit_hydrogens(), 0);
    assert_eq!(topology.atoms[1].prop("molAttachPoint"), None);
}

#[test]
fn every_typed_v2000_sgroup_record_resolves_ids_and_ordered_state() {
    let sdi = format!(
        "M  SDI   1  4{:>10.4}{:>10.4}{:>10.4}{:>10.4}",
        0.0, 1.0, 2.0, 3.0
    );
    let sbv = format!("M  SBV   1   1{:>10.4}{:>10.4}", 0.5, 0.25);
    let sdt = format!(
        "M  SDT   2 {:<30}{:<2}{:<20}{:<2}{}",
        "FIELD", "T", "INFO", "Q", "OP"
    );
    let input = molblock(
        &[atom_line("C", 0), atom_line("O", 0)],
        &[bond_line(1, 2, 1)],
        &format!(
            concat!(
                "M  STY  2   2 DAT   1 SUP\n",
                "M  SST  1   1 ALT\n",
                "M  SLB  1   1   7\n",
                "M  SCN  1   1 HT\n",
                "M  SDS EXP  1   1\n",
                "M  SAL   1  1   1\n",
                "M  SBL   1  1   1\n",
                "M  SPA   1  1   1\n",
                "M  SMT   1 Me\n",
                "{}\n{}\n{}\n",
                "M  SDD   2 display spec\n",
                "M  SCD   2 first value\n",
                "M  SED   2 second value\n",
                "M  SPL  1   2   1\n",
                "M  SNC  1   2   5\n",
                "M  SAP   1  1   1   2 AP\n",
                "M  SCL   2 CLASS\n",
                "M  SBT  1   2   1\n",
                "M  END\n",
            ),
            sdi, sbv, sdt
        ),
    );
    let (topology, _) = concrete(&input);
    assert_eq!(topology.substance_groups.len(), 2);
    let sup = &topology.substance_groups[0];
    assert_eq!(sup.id().index(), 0);
    assert_eq!(sup.rdkit_sequence_id(), Some(1));
    assert_eq!(sup.external_id(), Some(7));
    assert_eq!(sup.kind(), &SubstanceGroupKind::Superatom);
    assert_eq!(sup.atoms(), &[AtomId::new(0)]);
    assert_eq!(sup.parent_atoms(), &[AtomId::new(0)]);
    assert_eq!(sup.bonds(), &[BondId::new(0)]);
    assert_eq!(sup.label(), Some("Me"));
    assert_eq!(sup.subtype(), Some("ALT"));
    assert_eq!(sup.connection(), Some(&SGroupConnection::HeadToTail));
    assert_eq!(sup.expansion_state(), Some("E"));
    assert_eq!(sup.display().unwrap().brackets[0].p1, [0.0, 1.0]);
    assert_eq!(sup.cstates()[0].vector, [0.5, 0.25]);
    assert_eq!(sup.attach_points()[0].atom, AtomId::new(0));
    assert_eq!(sup.attach_points()[0].leaving_atom, Some(AtomId::new(1)));

    let dat = &topology.substance_groups[1];
    assert_eq!(dat.id().index(), 1);
    assert_eq!(dat.rdkit_sequence_id(), Some(2));
    assert_eq!(dat.parent(), Some(sup.id()));
    assert_eq!(dat.component_number(), Some(5));
    assert_eq!(dat.class(), Some("CLASS"));
    assert_eq!(dat.bracket_style(), Some(&SGroupBracketStyle::Parenthesis));
    let data = dat.data().unwrap();
    assert_eq!(data.field_name.as_deref(), Some("FIELD"));
    assert_eq!(data.field_type.as_deref(), Some("T"));
    assert_eq!(data.field_info.as_deref(), Some("INFO"));
    assert_eq!(data.query_type.as_deref(), Some("Q"));
    assert_eq!(data.query_op.as_deref(), Some("OP"));
    assert_eq!(data.field_display.as_deref(), Some("display spec"));
    assert_eq!(data.values, ["first valuesecond value"]);
}

#[test]
fn sgroup_strict_errors_and_non_strict_cleanup_never_keep_partial_groups() {
    let base = molblock(
        &[atom_line("C", 0), atom_line("C", 0)],
        &[bond_line(1, 2, 1)],
        concat!(
            "M  STY  2   1 SUP   2 DAT\n",
            "M  SAL   1  1   1\n",
            "M  SAL   2  1   2\n",
            "M  SBT  1   1   2\n",
            "M  END\n",
        ),
    );
    assert!(read_mol_block_detached(&base).is_err());
    let relaxed = read_mol_block_detached_with_params(
        &base,
        MolBlockReadParams {
            strict_parsing: false,
        },
    )
    .expect("non-strict invalid SGroup cleanup");
    let MolBlockRecord::Concrete { topology, .. } = relaxed else {
        panic!("expected concrete SGroup record")
    };
    assert_eq!(topology.substance_groups.len(), 1);
    assert_eq!(topology.substance_groups[0].rdkit_sequence_id(), Some(2));

    let crs = molblock(&[], &[], "M  CRS  1\nM  END\n");
    assert!(
        read_mol_block_detached(&crs)
            .unwrap_err()
            .to_string()
            .contains("Unsupported SGroup subtype 'M  CRS'")
    );

    let short_sap = molblock(
        &[atom_line("C", 0), atom_line("C", 0)],
        &[bond_line(1, 2, 1)],
        "M  STY  1   1 SUP\nM  SAL   1  1   1\nM  SBL   1  1   1\nM  SAP   1  1   1\nM  END\n",
    );
    assert!(read_mol_block_detached(&short_sap).is_err());
    let relaxed = read_mol_block_detached_with_params(
        &short_sap,
        MolBlockReadParams {
            strict_parsing: false,
        },
    )
    .expect("non-strict SAP inference");
    let MolBlockRecord::Concrete { topology, .. } = relaxed else {
        panic!("expected concrete SAP record")
    };
    assert_eq!(
        topology.substance_groups[0].attach_points()[0].leaving_atom,
        Some(AtomId::new(1))
    );
}

#[test]
fn data_sgroup_line_limit_and_property_blank_recovery_match_pinned_policy() {
    let too_many_scd = molblock(
        &[atom_line("C", 0)],
        &[],
        concat!(
            "M  STY  1   1 DAT\n",
            "M  SCD   1 one\n",
            "M  SCD   1 two\n",
            "M  SCD   1 three\n",
            "M  SCD   1 four\n",
            "M  SED   1 end\n",
            "M  END\n",
        ),
    );
    assert!(
        read_mol_block_detached(&too_many_scd)
            .unwrap_err()
            .to_string()
            .contains("too many consecutive SCD")
    );
    let relaxed = read_mol_block_detached_with_params(
        &too_many_scd,
        MolBlockReadParams {
            strict_parsing: false,
        },
    )
    .expect("non-strict SCD continuation");
    let MolBlockRecord::Concrete { topology, .. } = relaxed else {
        panic!("expected concrete DAT SGroup")
    };
    assert_eq!(
        topology.substance_groups[0].data().unwrap().values,
        ["onetwothreefourend"]
    );

    let leading_blank = molblock(&[], &[], "\nM  END\n");
    assert!(read_mol_block_detached(&leading_blank).is_err());
    assert!(
        read_mol_block_detached_with_params(
            &leading_blank,
            MolBlockReadParams {
                strict_parsing: false,
            }
        )
        .is_ok()
    );

    let later_blank = molblock(&[], &[], "M  CHG  0\n\nM  END\n");
    assert!(
        read_mol_block_detached_with_params(
            &later_blank,
            MolBlockReadParams {
                strict_parsing: false,
            }
        )
        .is_err()
    );
}
