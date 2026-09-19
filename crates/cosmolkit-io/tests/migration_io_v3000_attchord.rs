use cosmolkit_io::{MolBlockRecord, SdfReadError, read_mol_block_detached};

fn v3000_with_atoms(atoms: &[&str]) -> String {
    let mut block = format!(
        concat!(
            "ATTCHORD\n",
            "  COSMolKit\n",
            "\n",
            "  0  0  0  0  0  0  0  0  0  0999 V3000\n",
            "M  V30 BEGIN CTAB\n",
            "M  V30 COUNTS {} 0 0 0 0\n",
            "M  V30 BEGIN ATOM\n",
        ),
        atoms.len()
    );
    for atom in atoms {
        block.push_str("M  V30 ");
        block.push_str(atom);
        block.push('\n');
    }
    block.push_str("M  V30 END ATOM\nM  V30 END CTAB\nM  END\n");
    block
}

#[test]
fn v3000_template_attchord_preserves_pair_order_labels_and_forward_row_targets() {
    let block = v3000_with_atoms(&[
        "42 C 0 0 0 0 ATTCHORD=(4 3 Br 2 Al)",
        "7 N 1 0 0 0",
        "99 O 2 0 0 0",
    ]);

    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("template ATTCHORD")
    else {
        panic!("ordinary atoms must produce concrete topology");
    };
    let order = topology.atoms[0]
        .template_attachment_order()
        .expect("typed template attachment order");
    assert_eq!(order.entries().len(), 2);
    assert_eq!(order.entries()[0].target().index(), 2);
    assert_eq!(order.entries()[0].label(), "Br");
    assert_eq!(order.entries()[1].target().index(), 1);
    assert_eq!(order.entries()[1].label(), "Al");
    assert!(topology.atoms[1].template_attachment_order().is_none());
    assert!(topology.atoms[2].template_attachment_order().is_none());
}

#[test]
fn v3000_integer_attchord_remains_distinct_from_template_state() {
    let block = v3000_with_atoms(&["17 C 0 0 0 0 ATTCHORD=3"]);
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("ordinary integer ATTCHORD")
    else {
        panic!("ordinary atom must produce concrete topology");
    };
    let atom = &topology.atoms[0];
    assert_eq!(atom.prop("molAttachOrder"), Some("3"));
    assert!(atom.template_attachment_order().is_none());
}

#[test]
fn v3000_query_atom_carries_canonical_template_attchord_state() {
    let block = v3000_with_atoms(&[
        "10 * 0 0 0 0 ATTCHORD=(4 3 Cx 2 br)",
        "80 N 1 0 0 0",
        "20 O 2 0 0 0",
    ]);
    let MolBlockRecord::Query(record) =
        read_mol_block_detached(&block).expect("query template ATTCHORD")
    else {
        panic!("wildcard atom must remain query topology");
    };
    let order = record.query.atoms()[0]
        .atom()
        .template_attachment_order()
        .expect("query carrier preserves canonical atom state");
    let observed = order
        .entries()
        .iter()
        .map(|entry| (entry.target().index(), entry.label()))
        .collect::<Vec<_>>();
    assert_eq!(observed, vec![(2, "Cx"), (1, "br")]);
}

#[test]
fn v3000_template_attchord_rejects_invalid_shape_and_duplicate_pairs() {
    for (value, reason) in [
        ("(3 2 Al 3)", "odd item count"),
        ("(4 2 Al)", "declared count mismatch"),
        ("(4 2 Al 2 Br)", "duplicate target"),
        ("(4 2 Al 3 Al)", "duplicate label"),
        ("(2 0 Al)", "zero row"),
        ("(2 nope Al)", "nonnumeric row"),
    ] {
        let block = v3000_with_atoms(&[
            &format!("1 C 0 0 0 0 ATTCHORD={value}"),
            "2 N 1 0 0 0",
            "3 O 2 0 0 0",
        ]);
        assert!(
            matches!(
                read_mol_block_detached(&block),
                Err(SdfReadError::Parse(message)) if message.contains("Invalid ATTCHORD value")
            ),
            "{reason} must be a structured parse failure"
        );
    }
}

#[test]
fn v3000_template_attchord_rejects_out_of_range_after_complete_atom_table() {
    let concrete = v3000_with_atoms(&["5 C 0 0 0 0 ATTCHORD=(2 3 Al)", "90 N 1 0 0 0"]);
    assert!(matches!(
        read_mol_block_detached(&concrete),
        Err(SdfReadError::Topology(_))
    ));

    let query = v3000_with_atoms(&["5 * 0 0 0 0 ATTCHORD=(2 3 Al)", "90 N 1 0 0 0"]);
    assert!(matches!(
        read_mol_block_detached(&query),
        Err(SdfReadError::QueryGraph(_))
    ));
}

#[test]
fn v3000_template_attchord_rejects_repeated_separators_like_the_source() {
    // `boost::split` keeps the empty field between repeated separators, so the
    // token count no longer matches `itemCount + 1` and the record fails.
    for value in ["(2  2 Al)", "(2\t\t2 Al)", "(2 \t 2 Al)"] {
        let block = v3000_with_atoms(&[&format!("1 C 0 0 0 0 ATTCHORD={value}"), "2 N 1 0 0 0"]);
        assert!(
            matches!(
                read_mol_block_detached(&block),
                Err(SdfReadError::Parse(message)) if message.contains("Invalid ATTCHORD value")
            ),
            "repeated separators in {value} must be rejected"
        );
    }
}

#[test]
fn v3000_template_attchord_accepts_empty_labels_and_mixed_separators() {
    // A trailing separator supplies the final (empty) label exactly like the
    // pinned branch; tabs are delimiters just like spaces.
    let block = v3000_with_atoms(&[
        "1 C 0 0 0 0 ATTCHORD=(2 2 )",
        "2 N 1 0 0 0",
        "3 O 2 0 0 0 ATTCHORD=(6 3 Br\t2 Al 1 )",
    ]);
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("empty labels are supported source behavior")
    else {
        panic!("ordinary atoms must produce concrete topology");
    };
    let first = topology.atoms[0]
        .template_attachment_order()
        .expect("first carrier");
    assert_eq!(first.entries().len(), 1);
    assert_eq!(first.entries()[0].target().index(), 1);
    assert_eq!(first.entries()[0].label(), "");

    let second = topology.atoms[2]
        .template_attachment_order()
        .expect("second carrier");
    let observed = second
        .entries()
        .iter()
        .map(|entry| (entry.target().index(), entry.label()))
        .collect::<Vec<_>>();
    assert_eq!(observed, vec![(2, "Br"), (1, "Al"), (0, "")]);
}

#[test]
fn v3000_template_attchord_rejects_zero_negative_and_overflow_indices() {
    for (value, reason) in [
        ("(2 0 Al)", "zero row is rejected before subtraction"),
        (
            "(2 -1 Al)",
            "negative rows wrap in the source and fail closed here",
        ),
        (
            "(2 -2 Al)",
            "negative rows wrap in the source and fail closed here",
        ),
        ("(4294967296 1 Al)", "item-count overflow resolves to zero"),
        ("(2 4294967296 Al)", "index overflow resolves to zero"),
        ("(2 +0 Al)", "explicit positive zero is still zero"),
    ] {
        let block = v3000_with_atoms(&[&format!("1 C 0 0 0 0 ATTCHORD={value}"), "2 N 1 0 0 0"]);
        assert!(
            matches!(
                read_mol_block_detached(&block),
                Err(SdfReadError::Parse(message)) if message.contains("Invalid ATTCHORD value")
            ),
            "{reason}: {value} must be a structured parse failure"
        );
    }
}

#[test]
fn v3000_template_attchord_rejects_empty_and_all_space_inner_values() {
    for value in ["()", "(   )"] {
        let block = v3000_with_atoms(&[&format!("1 C 0 0 0 0 ATTCHORD={value}"), "2 N 1 0 0 0"]);
        assert!(
            matches!(
                read_mol_block_detached(&block),
                Err(SdfReadError::Parse(message)) if message.contains("Invalid ATTCHORD value")
            ),
            "empty inner value {value} must be rejected"
        );
    }
}

#[test]
fn v3000_template_attchord_keeps_source_literal_last_character_strip() {
    // The source strips the first and last characters without checking for
    // the closing parenthesis, so `(2 2 Al` parses with label `A`; malformed
    // input must not be normalized into anything else.
    let block = v3000_with_atoms(&["1 C 0 0 0 0 ATTCHORD=(2 2 Al", "2 N 1 0 0 0"]);
    let MolBlockRecord::Concrete { topology, .. } =
        read_mol_block_detached(&block).expect("source accepts the literal last-char strip")
    else {
        panic!("ordinary atoms must produce concrete topology");
    };
    let order = topology.atoms[0]
        .template_attachment_order()
        .expect("carrier from the source-shaped strip");
    assert_eq!(order.entries().len(), 1);
    assert_eq!(order.entries()[0].target().index(), 1);
    assert_eq!(order.entries()[0].label(), "A");
}
