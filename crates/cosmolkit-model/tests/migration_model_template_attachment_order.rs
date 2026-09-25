use std::collections::BTreeMap;

use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Element, QueryAtom, QueryGraph, QueryGraphError, TemplateAttachment,
    TemplateAttachmentOrder, TemplateAttachmentOrderError, TopologyBlock, TopologyEditError,
    TopologyValidationError,
};

fn attachment_order() -> TemplateAttachmentOrder {
    TemplateAttachmentOrder::new(vec![
        TemplateAttachment::new(AtomId::new(2), "Al"),
        TemplateAttachment::new(AtomId::new(1), "Br"),
    ])
    .expect("fixture attachment order must be valid")
}

fn atom(id: usize, spec: AtomSpec) -> Atom {
    Atom::from_spec(AtomId::new(id), spec)
}

#[test]
fn template_attachment_order_preserves_order_labels_and_clone_state() {
    let order = attachment_order();
    let spec = AtomSpec::new(Element::C).with_template_attachment_order(order.clone());
    let concrete = atom(0, spec.clone());
    let cloned = concrete.clone();

    assert_eq!(
        order
            .entries()
            .iter()
            .map(|entry| (entry.target(), entry.label()))
            .collect::<Vec<_>>(),
        vec![(AtomId::new(2), "Al"), (AtomId::new(1), "Br")]
    );
    assert_eq!(spec.template_attachment_order(), Some(&order));
    assert_eq!(concrete.template_attachment_order(), Some(&order));
    assert_eq!(cloned.template_attachment_order(), Some(&order));
}

#[test]
fn template_attachment_order_rejects_empty_and_duplicate_targets() {
    assert_eq!(
        TemplateAttachmentOrder::new(Vec::new()),
        Err(TemplateAttachmentOrderError::Empty)
    );
    assert_eq!(
        TemplateAttachmentOrder::new(vec![
            TemplateAttachment::new(AtomId::new(1), "Al"),
            TemplateAttachment::new(AtomId::new(1), "Br"),
        ]),
        Err(TemplateAttachmentOrderError::DuplicateTarget {
            first_position: 0,
            duplicate_position: 1,
            target: AtomId::new(1),
        })
    );
}

#[test]
fn template_attachment_order_rejects_duplicate_labels_without_rejecting_empty_labels() {
    assert_eq!(
        TemplateAttachmentOrder::new(vec![
            TemplateAttachment::new(AtomId::new(1), "Al"),
            TemplateAttachment::new(AtomId::new(2), "Al"),
        ]),
        Err(TemplateAttachmentOrderError::DuplicateLabel {
            first_position: 0,
            duplicate_position: 1,
            label: "Al".to_owned(),
        })
    );

    let empty_label =
        TemplateAttachmentOrder::new(vec![TemplateAttachment::new(AtomId::new(0), "")])
            .expect("the pinned parser does not reject an otherwise unique empty label");
    assert_eq!(empty_label.entries()[0].label(), "");
}

#[test]
fn concrete_topology_rejects_template_attachment_targets_outside_atom_table() {
    let invalid = TemplateAttachmentOrder::new(vec![TemplateAttachment::new(AtomId::new(2), "Al")])
        .expect("local order is valid before graph-range validation");
    let result = TopologyBlock::try_from_parts(
        vec![
            atom(
                0,
                AtomSpec::new(Element::C).with_template_attachment_order(invalid),
            ),
            atom(1, AtomSpec::new(Element::C)),
        ],
        Vec::new(),
        Vec::new(),
        Vec::new(),
    );

    assert!(matches!(
        result,
        Err(TopologyValidationError::TemplateAttachmentOrder {
            atom,
            source: TemplateAttachmentOrderError::TargetOutOfRange {
                position: 0,
                target,
                atom_count: 2,
            },
        }) if atom == AtomId::new(0) && target == AtomId::new(2)
    ));
}

#[test]
fn query_atoms_share_the_concrete_carrier_and_range_validation() {
    let order = attachment_order();
    let graph = QueryGraph::from_parts(
        vec![
            QueryAtom::new(
                AtomId::new(0),
                AtomSpec::new(Element::C).with_template_attachment_order(order.clone()),
            ),
            QueryAtom::new(AtomId::new(1), AtomSpec::new(Element::C)),
            QueryAtom::new(AtomId::new(2), AtomSpec::new(Element::C)),
        ],
        Vec::new(),
        BTreeMap::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .expect("query graph uses the same valid atom-table references");
    assert_eq!(graph.atoms()[0].template_attachment_order(), Some(&order));
    assert_eq!(graph.atoms()[0].template_attachment_order(), Some(&order));

    let invalid = TemplateAttachmentOrder::new(vec![TemplateAttachment::new(AtomId::new(3), "Al")])
        .expect("local order is valid before graph-range validation");
    let result = QueryGraph::from_parts(
        vec![QueryAtom::new(
            AtomId::new(0),
            AtomSpec::new(Element::C).with_template_attachment_order(invalid),
        )],
        Vec::new(),
        BTreeMap::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    );
    assert!(matches!(
        result,
        Err(QueryGraphError::TemplateAttachmentOrder {
            atom,
            source: TemplateAttachmentOrderError::TargetOutOfRange {
                position: 0,
                target,
                atom_count: 1,
            },
        }) if atom == AtomId::new(0) && target == AtomId::new(3)
    ));
}

#[test]
fn topology_reorder_remaps_carrier_targets_and_preserves_entry_order() {
    let source = TopologyBlock::try_from_parts(
        vec![
            atom(
                0,
                AtomSpec::new(Element::C).with_template_attachment_order(attachment_order()),
            ),
            atom(1, AtomSpec::new(Element::N)),
            atom(2, AtomSpec::new(Element::O)),
        ],
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .expect("source topology is valid");

    let (reordered, mapping) = source
        .reordered_atoms(&[AtomId::new(2), AtomId::new(0), AtomId::new(1)])
        .expect("complete permutation must preserve attachment references");
    assert_eq!(
        mapping.atoms.old_to_new,
        vec![
            Some(AtomId::new(1)),
            Some(AtomId::new(2)),
            Some(AtomId::new(0))
        ]
    );
    assert_eq!(
        reordered.atoms[1]
            .template_attachment_order()
            .expect("carrier state is preserved")
            .entries()
            .iter()
            .map(|entry| (entry.target(), entry.label()))
            .collect::<Vec<_>>(),
        vec![(AtomId::new(0), "Al"), (AtomId::new(2), "Br")]
    );
}

#[test]
fn batch_delete_rejects_a_surviving_carrier_with_a_removed_target_atomically() {
    let source = TopologyBlock::try_from_parts(
        vec![
            atom(
                0,
                AtomSpec::new(Element::C).with_template_attachment_order(
                    TemplateAttachmentOrder::new(vec![TemplateAttachment::new(
                        AtomId::new(2),
                        "Al",
                    )])
                    .unwrap(),
                ),
            ),
            atom(1, AtomSpec::new(Element::N)),
            atom(2, AtomSpec::new(Element::O)),
        ],
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .unwrap();
    let snapshot = source.clone();
    let mut edit = source.begin_batch_edit().unwrap();
    edit.remove_atom(AtomId::new(2)).unwrap();

    assert!(matches!(
        edit.finish(),
        Err(TopologyEditError::TemplateAttachmentRemap {
            carrier,
            source: TemplateAttachmentOrderError::TargetRemoved {
                position: 0,
                target,
            },
        }) if carrier == AtomId::new(0) && target == AtomId::new(2)
    ));
    assert_eq!(
        source, snapshot,
        "detached edit failure cannot alter the source"
    );
    source.validate().unwrap();
}

#[test]
fn batch_delete_of_the_carrier_drops_its_attachment_state_with_the_row() {
    let source = TopologyBlock::try_from_parts(
        vec![
            atom(
                0,
                AtomSpec::new(Element::C).with_template_attachment_order(attachment_order()),
            ),
            atom(1, AtomSpec::new(Element::N)),
            atom(2, AtomSpec::new(Element::O)),
        ],
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .unwrap();
    let mut edit = source.begin_batch_edit().unwrap();
    edit.remove_atom(AtomId::new(0)).unwrap();

    let (result, mapping) = edit
        .finish()
        .expect("removing the carrier removes its state");
    assert_eq!(mapping.atoms.old_to_new[0], None);
    assert_eq!(result.atoms.len(), 2);
    assert!(
        result
            .atoms
            .iter()
            .all(|atom| atom.template_attachment_order().is_none())
    );
    result.validate().unwrap();
}

#[test]
fn shared_remap_supports_offset_combination_maps_without_a_parallel_combine_api() {
    let order = TemplateAttachmentOrder::new(vec![TemplateAttachment::new(AtomId::new(1), "tail")])
        .unwrap();
    let remapped = order
        .remapped(&[Some(AtomId::new(3)), Some(AtomId::new(4))])
        .expect("an append/combination offset map is a complete reference map");
    assert_eq!(remapped.entries()[0].target(), AtomId::new(4));
    assert_eq!(remapped.entries()[0].label(), "tail");
}

#[test]
fn topology_edits_without_template_attachment_state_keep_existing_behavior() {
    let source = TopologyBlock::try_from_parts(
        vec![
            atom(0, AtomSpec::new(Element::C)),
            atom(1, AtomSpec::new(Element::N)),
        ],
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .unwrap();
    let (reordered, _) = source
        .reordered_atoms(&[AtomId::new(1), AtomId::new(0)])
        .unwrap();
    assert_eq!(reordered.atoms[0].element(), Element::N);
    assert_eq!(reordered.atoms[1].element(), Element::C);

    let mut edit = source.begin_batch_edit().unwrap();
    edit.remove_atom(AtomId::new(1)).unwrap();
    let (trimmed, _) = edit.finish().unwrap();
    assert_eq!(trimmed.atoms.len(), 1);
    assert_eq!(trimmed.atoms[0].element(), Element::C);
}
