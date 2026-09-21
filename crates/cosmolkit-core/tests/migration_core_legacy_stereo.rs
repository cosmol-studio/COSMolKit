#![cfg(feature = "op-contracts-strict")]

use cosmolkit_core::{
    RingSearchParams, ValenceModel, assign_legacy_stereochemistry, assign_valence_for_topology,
    symmetrized_sssr,
};
use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, StereoGroup, StereoGroupKind, TopologyBlock,
};
use cosmolkit_types::{BondDirection, BondOrder, BondStereo, ChiralTag, Element};

fn topology_from_specs(
    atoms: Vec<AtomSpec>,
    specs: Vec<BondSpec>,
    groups: Vec<StereoGroup>,
) -> TopologyBlock {
    TopologyBlock::try_from_parts(
        atoms
            .into_iter()
            .enumerate()
            .map(|(index, spec)| Atom::from_spec(AtomId::new(index), spec))
            .collect(),
        specs
            .into_iter()
            .enumerate()
            .map(|(index, spec)| Bond::from_spec(BondId::new(index), spec))
            .collect(),
        Vec::new(),
        groups,
    )
    .unwrap()
}

fn topology(atom_count: usize, specs: Vec<BondSpec>, groups: Vec<StereoGroup>) -> TopologyBlock {
    topology_from_specs(
        (0..atom_count).map(|_| AtomSpec::new(Element::C)).collect(),
        specs,
        groups,
    )
}

fn finish(input: TopologyBlock) -> TopologyBlock {
    let valence = assign_valence_for_topology(&input, ValenceModel::RdkitLike).unwrap();
    let rings = symmetrized_sssr(&input, &RingSearchParams::default()).unwrap();
    assign_legacy_stereochemistry(input, &valence, &rings).unwrap()
}

fn three_coordinate_center_with_graph_hydrogen(isotope: u16) -> TopologyBlock {
    topology_from_specs(
        vec![
            AtomSpec::new(Element::C)
                .with_explicit_hydrogens(1)
                .with_no_implicit(true)
                .with_chiral_tag(ChiralTag::TetrahedralCw),
            AtomSpec::new(Element::H).with_isotope(isotope),
            AtomSpec::new(Element::F),
            AtomSpec::new(Element::CL),
        ],
        vec![
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            BondSpec::new(AtomId::new(0), AtomId::new(2), BondOrder::Single),
            BondSpec::new(AtomId::new(0), AtomId::new(3), BondOrder::Single),
        ],
        vec![],
    )
}

#[test]
fn mol_post_legacy_audit_isotope_zero_graph_hydrogen_is_protium() {
    let protium = finish(three_coordinate_center_with_graph_hydrogen(0));
    assert_eq!(protium.atoms[0].chiral_tag(), ChiralTag::Unspecified);
    assert_eq!(protium.atoms[0].prop("_CIPCode"), None);

    let deuterium = finish(three_coordinate_center_with_graph_hydrogen(2));
    assert_eq!(deuterium.atoms[0].chiral_tag(), ChiralTag::TetrahedralCw);
    assert!(matches!(
        deuterium.atoms[0].prop("_CIPCode"),
        Some("R" | "S")
    ));
}

#[test]
fn mol_post_legacy_audit_inert_graph_does_not_materialize_cip_ranks() {
    let output = finish(topology(
        2,
        vec![BondSpec::new(
            AtomId::new(0),
            AtomId::new(1),
            BondOrder::Single,
        )],
        vec![],
    ));
    assert!(
        output
            .atoms
            .iter()
            .all(|atom| atom.prop("_CIPRank").is_none())
    );
    assert!(
        output
            .atoms
            .iter()
            .all(|atom| atom.prop("_ChiralityPossible").is_none())
    );
}

#[test]
fn mol_post_legacy_audit_potential_center_materializes_lazy_rank_state() {
    let output = finish(topology_from_specs(
        vec![
            AtomSpec::new(Element::C)
                .with_explicit_hydrogens(1)
                .with_no_implicit(true),
            AtomSpec::new(Element::F),
            AtomSpec::new(Element::CL),
            AtomSpec::new(Element::BR),
        ],
        vec![
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            BondSpec::new(AtomId::new(0), AtomId::new(2), BondOrder::Single),
            BondSpec::new(AtomId::new(0), AtomId::new(3), BondOrder::Single),
        ],
        vec![],
    ));
    assert_eq!(output.atoms[0].prop("_ChiralityPossible"), Some("1"));
    assert!(
        output
            .atoms
            .iter()
            .all(|atom| atom.prop("_CIPRank").is_some())
    );
}

#[test]
fn mol_post_legacy_closure_cleans_stale_wedge_but_preserves_atrop_control() {
    let stale = topology(
        2,
        vec![
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single)
                .with_direction(BondDirection::BeginWedge),
        ],
        vec![],
    );
    assert_eq!(finish(stale).bonds[0].direction(), BondDirection::None);

    let atrop_control = topology(
        3,
        vec![
            BondSpec::new(AtomId::new(1), AtomId::new(0), BondOrder::Single)
                .with_direction(BondDirection::BeginDash),
            BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Single)
                .with_stereo(BondStereo::AtropCw),
        ],
        vec![],
    );
    assert_eq!(
        finish(atrop_control).bonds[0].direction(),
        BondDirection::BeginDash
    );
}

#[test]
fn mol_post_legacy_closure_cleans_unused_direction_around_unresolved_double_bond() {
    let input = topology(
        4,
        vec![
            BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Double),
            BondSpec::new(AtomId::new(1), AtomId::new(0), BondOrder::Single)
                .with_direction(BondDirection::EndUpRight),
            BondSpec::new(AtomId::new(2), AtomId::new(3), BondOrder::Single),
        ],
        vec![],
    );
    let output = finish(input);
    assert_eq!(output.bonds[0].stereo(), BondStereo::None);
    assert_eq!(output.bonds[1].direction(), BondDirection::None);
}

#[test]
fn mol_post_legacy_closure_cleans_atrop_groups_before_general_groups() {
    let input = topology(
        2,
        vec![
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single)
                .with_stereo(BondStereo::AtropCcw),
        ],
        vec![
            StereoGroup::new(
                StereoGroupKind::Absolute,
                vec![AtomId::new(0), AtomId::new(1)],
                vec![],
            )
            .with_id(17),
        ],
    );
    let output = finish(input);
    assert_eq!(output.stereo_groups.len(), 1);
    assert_eq!(output.stereo_groups[0].atoms(), &[]);
    assert_eq!(output.stereo_groups[0].bonds(), &[BondId::new(0)]);
    assert_eq!(output.stereo_groups[0].id(), None);
}
