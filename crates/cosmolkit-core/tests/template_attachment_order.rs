use cosmolkit_core::{DetachedPathSubgraph, PathError, SubtopologyParams, subtopology_from_path};
use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, TemplateAttachment, TemplateAttachmentOrder,
    TemplateAttachmentOrderError, TopologyBlock,
};
use cosmolkit_types::{BondOrder, Element};

fn source_topology() -> TopologyBlock {
    let order = TemplateAttachmentOrder::new(vec![TemplateAttachment::new(AtomId::new(2), "port")])
        .unwrap();
    let atoms = vec![
        Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C).with_template_attachment_order(order),
        ),
        Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::N)),
        Atom::from_spec(AtomId::new(2), AtomSpec::new(Element::O)),
        Atom::from_spec(AtomId::new(3), AtomSpec::new(Element::S)),
    ];
    let bonds = [(0, 1), (0, 2), (2, 3)]
        .into_iter()
        .enumerate()
        .map(|(index, (begin, end))| {
            Bond::from_spec(
                BondId::new(index),
                BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
            )
        })
        .collect();
    TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap()
}

fn assert_remapped_carrier(subgraph: &DetachedPathSubgraph) {
    let order = match subgraph {
        DetachedPathSubgraph::Concrete(topology) => topology.atoms[0]
            .template_attachment_order()
            .expect("concrete carrier preserves typed state"),
        DetachedPathSubgraph::Query { graph, .. } => graph.atoms()[0]
            .template_attachment_order()
            .expect("query carrier preserves typed state"),
    };
    assert_eq!(order.entries()[0].target(), AtomId::new(1));
    assert_eq!(order.entries()[0].label(), "port");
}

#[test]
fn extraction_remaps_template_attachment_targets_for_concrete_and_query_paths() {
    let source = source_topology();
    let snapshot = source.clone();

    for copy_as_query in [false, true] {
        let result = subtopology_from_path(
            &source,
            &[BondId::new(1)],
            &SubtopologyParams { copy_as_query },
        )
        .expect("selected carrier and target form a valid detached subgraph");
        assert_eq!(
            result.mapping.atoms.old_to_new,
            vec![Some(AtomId::new(0)), None, Some(AtomId::new(1)), None]
        );
        assert_remapped_carrier(&result.subgraph);
    }

    assert_eq!(
        source, snapshot,
        "extraction cannot mutate the source topology"
    );
}

#[test]
fn extraction_fails_atomically_when_a_surviving_carrier_loses_its_target() {
    let source = source_topology();
    let snapshot = source.clone();

    for copy_as_query in [false, true] {
        assert!(matches!(
            subtopology_from_path(
                &source,
                &[BondId::new(0)],
                &SubtopologyParams { copy_as_query },
            ),
            Err(PathError::TemplateAttachmentRemap {
                carrier,
                source: TemplateAttachmentOrderError::TargetRemoved {
                    position: 0,
                    target,
                },
            }) if carrier == AtomId::new(0) && target == AtomId::new(2)
        ));
    }

    assert_eq!(source, snapshot);
    source.validate().unwrap();
}

#[test]
fn extraction_without_template_attachment_state_keeps_existing_behavior() {
    let atoms = vec![
        Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
        Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::O)),
    ];
    let bonds = vec![Bond::from_spec(
        BondId::new(0),
        BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Double),
    )];
    let source = TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap();

    let result =
        subtopology_from_path(&source, &[BondId::new(0)], &SubtopologyParams::default()).unwrap();
    let DetachedPathSubgraph::Concrete(result) = result.subgraph else {
        panic!("default extraction must remain concrete");
    };
    assert_eq!(result.atoms.len(), 2);
    assert_eq!(result.bonds.len(), 1);
    result.validate().unwrap();
}
