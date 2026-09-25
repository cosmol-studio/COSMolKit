#![cfg(feature = "op-contracts-strict")]

use cosmolkit_core::{
    LegacyStereoError, RingSearchParams, ValenceModel, assign_legacy_stereochemistry,
    assign_legacy_stereochemistry_for_depiction, assign_legacy_stereochemistry_with_flags,
    assign_valence_for_topology, cleanup_stereo_groups, symmetrized_sssr,
};
use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, StereoGroup, StereoGroupKind, TopologyBlock,
};
use cosmolkit_types::{BondOrder, BondStereo, ChiralTag, Element};

fn topology_from_specs(atoms: Vec<AtomSpec>, bonds: Vec<BondSpec>) -> TopologyBlock {
    topology_with_stereo_groups(atoms, bonds, Vec::new())
}

fn topology_with_stereo_groups(
    atoms: Vec<AtomSpec>,
    bonds: Vec<BondSpec>,
    stereo_groups: Vec<StereoGroup>,
) -> TopologyBlock {
    TopologyBlock::try_from_parts(
        atoms
            .into_iter()
            .enumerate()
            .map(|(index, spec)| Atom::from_spec(AtomId::new(index), spec))
            .collect(),
        bonds
            .into_iter()
            .enumerate()
            .map(|(index, spec)| Bond::from_spec(BondId::new(index), spec))
            .collect(),
        Vec::new(),
        stereo_groups,
    )
    .unwrap()
}

fn raw_large_ring_cis(with_stale_properties: bool) -> TopologyBlock {
    // Matched raw topology from pinned RDKit 2026.03.1 input
    // `C1CCCCC=CCCC1 |c:5|` parsed with sanitize=false/removeHs=false.
    // RDKit has a ten-membered ring, Cis on bond 5, stereo atoms (4, 7), and
    // no directional single bonds at rows 4 or 6. Under the fixed legacy
    // profile, cleanStereo=true clears that Cis state before CX extension
    // output; cleanStereo=false retains it. The upstream oracle profile and
    // exact strings remain frozen in the owning SMILES regression.
    let atoms = (0..10)
        .map(|index| {
            let mut atom = AtomSpec::new(Element::C);
            if with_stale_properties && index == 0 {
                atom = atom
                    .with_prop("_CIPCode", "S")
                    .unwrap()
                    .with_prop("_ChiralityPossible", "0")
                    .unwrap()
                    .with_prop("_ringStereochemCand", "1")
                    .unwrap()
                    .with_prop("_ringStereoAtoms", "4,7")
                    .unwrap()
                    .with_prop("user_marker", "keep")
                    .unwrap();
            }
            atom
        })
        .collect();
    let bonds = (0..10)
        .map(|index| {
            let next = (index + 1) % 10;
            let order = if index == 5 {
                BondOrder::Double
            } else {
                BondOrder::Single
            };
            let mut bond = BondSpec::new(AtomId::new(index), AtomId::new(next), order);
            if index == 5 {
                bond = bond
                    .with_stereo(BondStereo::Cis)
                    .with_stereo_atoms(AtomId::new(4), AtomId::new(7));
                if with_stale_properties {
                    bond = bond
                        .with_prop("_CIPCode", "Z")
                        .unwrap()
                        .with_prop("user_marker", "keep")
                        .unwrap();
                }
            }
            bond
        })
        .collect();
    topology_from_specs(atoms, bonds)
}

fn potential_tetrahedral_center() -> TopologyBlock {
    topology_from_specs(
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
    )
}

fn assign_with_flags(
    topology: TopologyBlock,
    clean_it: bool,
    flag_possible_stereo_centers: bool,
) -> Result<TopologyBlock, LegacyStereoError> {
    let valence = assign_valence_for_topology(&topology, ValenceModel::RdkitLike).unwrap();
    let rings = symmetrized_sssr(&topology, &RingSearchParams::default()).unwrap();
    assign_legacy_stereochemistry_with_flags(
        topology,
        &valence,
        &rings,
        clean_it,
        flag_possible_stereo_centers,
    )
}

fn assign_fixed_clean(topology: TopologyBlock) -> TopologyBlock {
    let valence = assign_valence_for_topology(&topology, ValenceModel::RdkitLike).unwrap();
    let rings = symmetrized_sssr(&topology, &RingSearchParams::default()).unwrap();
    assign_legacy_stereochemistry(topology, &valence, &rings).unwrap()
}

fn assign_fixed_depiction(topology: TopologyBlock) -> TopologyBlock {
    let valence = assign_valence_for_topology(&topology, ValenceModel::RdkitLike).unwrap();
    let rings = symmetrized_sssr(&topology, &RingSearchParams::default()).unwrap();
    assign_legacy_stereochemistry_for_depiction(topology, &valence, &rings).unwrap()
}

#[test]
fn cx_cleanup_stereo_groups_keeps_valid_members_empty_groups_ids_and_order() {
    // Pinned RDKit Chirality.cpp::cleanupStereoGroups keeps a group unchanged
    // when every listed atom and bond is valid, including empty groups.
    let groups = vec![
        StereoGroup::new(
            StereoGroupKind::Or,
            vec![AtomId::new(0), AtomId::new(1)],
            vec![BondId::new(0)],
        )
        .with_id(7),
        StereoGroup::new(StereoGroupKind::And, vec![AtomId::new(2)], Vec::new()).with_id(19),
        StereoGroup::new(StereoGroupKind::Absolute, Vec::new(), Vec::new()).with_id(23),
    ];
    let mut topology = topology_with_stereo_groups(
        vec![
            AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::TetrahedralCw),
            AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::TetrahedralCcw),
            AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::TetrahedralCw),
        ],
        vec![
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single)
                .with_stereo(BondStereo::AtropCw),
        ],
        groups,
    );
    let before = topology.clone();

    cleanup_stereo_groups(&mut topology);

    assert_eq!(
        topology, before,
        "fully valid groups remain in source order"
    );
    assert_eq!(topology.stereo_groups[0].id(), Some(7));
    assert_eq!(topology.stereo_groups[1].id(), Some(19));
    assert_eq!(topology.stereo_groups[2].id(), Some(23));
}

#[test]
fn cx_cleanup_stereo_groups_filters_mixed_members_and_preserves_read_id_order() {
    let mut topology = topology_with_stereo_groups(
        vec![
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::TetrahedralCw),
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::TetrahedralCcw),
            AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::TetrahedralCw),
        ],
        vec![
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Single)
                .with_stereo(BondStereo::AtropCw),
            BondSpec::new(AtomId::new(2), AtomId::new(3), BondOrder::Single),
            BondSpec::new(AtomId::new(3), AtomId::new(4), BondOrder::Single)
                .with_stereo(BondStereo::AtropCcw),
        ],
        vec![
            StereoGroup::new(
                StereoGroupKind::Or,
                vec![AtomId::new(0), AtomId::new(3), AtomId::new(1)],
                vec![
                    BondId::new(2),
                    BondId::new(3),
                    BondId::new(1),
                    BondId::new(0),
                ],
            )
            .with_id(41),
        ],
    );
    let before = topology.clone();

    cleanup_stereo_groups(&mut topology);

    let expected_group = StereoGroup::new(
        StereoGroupKind::Or,
        vec![AtomId::new(3), AtomId::new(1)],
        vec![BondId::new(3), BondId::new(1)],
    )
    .with_id(41);
    assert_eq!(topology.stereo_groups, vec![expected_group]);
    assert_eq!(topology.atoms, before.atoms);
    assert_eq!(topology.bonds, before.bonds);
    assert_eq!(topology.adjacency, before.adjacency);
}

#[test]
fn cx_cleanup_stereo_groups_drops_invalid_groups_but_keeps_empty_group() {
    let empty_group =
        StereoGroup::new(StereoGroupKind::Absolute, Vec::new(), Vec::new()).with_id(8);
    let mut topology = topology_with_stereo_groups(
        vec![
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
        ],
        vec![
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            BondSpec::new(AtomId::new(2), AtomId::new(3), BondOrder::Single)
                .with_stereo(BondStereo::AtropCw),
        ],
        vec![
            StereoGroup::new(StereoGroupKind::Or, vec![AtomId::new(0)], Vec::new()).with_id(5),
            StereoGroup::new(StereoGroupKind::And, Vec::new(), vec![BondId::new(0)]).with_id(6),
            StereoGroup::new(
                StereoGroupKind::Or,
                vec![AtomId::new(1)],
                vec![BondId::new(1)],
            )
            .with_id(7),
            empty_group.clone(),
        ],
    );

    cleanup_stereo_groups(&mut topology);

    assert_eq!(topology.stereo_groups, vec![empty_group]);
}

#[test]
fn cx_cleanup_stereo_groups_accepts_only_atrop_bond_only_groups() {
    let valid_group =
        StereoGroup::new(StereoGroupKind::And, Vec::new(), vec![BondId::new(0)]).with_id(77);
    let mut topology = topology_with_stereo_groups(
        vec![
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
        ],
        vec![
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single)
                .with_stereo(BondStereo::AtropCcw),
            BondSpec::new(AtomId::new(2), AtomId::new(3), BondOrder::Single),
        ],
        vec![
            valid_group.clone(),
            StereoGroup::new(StereoGroupKind::Or, Vec::new(), vec![BondId::new(1)]).with_id(78),
        ],
    );

    cleanup_stereo_groups(&mut topology);

    assert_eq!(topology.stereo_groups, vec![valid_group]);
}

#[test]
fn cx_legacy_profile_cleans_the_exact_raw_large_ring_cis_state() {
    let input = raw_large_ring_cis(false);
    let before = input.clone();
    assert_eq!(input.atoms.len(), 10);
    assert_eq!(input.bonds.len(), 10);
    assert_eq!(input.bonds[5].stereo(), BondStereo::Cis);
    assert_eq!(
        input.bonds[5].stereo_atoms(),
        Some([AtomId::new(4), AtomId::new(7)])
    );

    let clean_false = assign_with_flags(input.clone(), false, false).unwrap();
    let old_false_false = assign_fixed_depiction(input.clone());
    assert_eq!(clean_false, old_false_false);
    assert_eq!(clean_false.bonds[5].stereo(), BondStereo::Cis);
    assert_eq!(
        clean_false.bonds[5].stereo_atoms(),
        Some([AtomId::new(4), AtomId::new(7)])
    );

    let cx_profile = assign_with_flags(input.clone(), true, false).unwrap();
    assert_eq!(cx_profile.bonds[5].stereo(), BondStereo::None);
    assert_eq!(cx_profile.bonds[5].stereo_atoms(), None);

    let old_true_true = assign_fixed_clean(input.clone());
    let exact_flags_true_true = assign_with_flags(input.clone(), true, true).unwrap();
    assert_eq!(old_true_true, exact_flags_true_true);
    assert_eq!(old_true_true.bonds[5].stereo(), BondStereo::None);
    assert_eq!(
        input, before,
        "detached assignment must leave its input intact"
    );
}

#[test]
fn cx_legacy_possible_center_flag_is_independent_of_cleaning() {
    let input = potential_tetrahedral_center();
    let before = input.clone();

    let no_possible_centers = assign_with_flags(input.clone(), false, false).unwrap();
    let possible_centers = assign_with_flags(input.clone(), false, true).unwrap();
    assert_eq!(
        no_possible_centers.atoms[0].prop("_ChiralityPossible"),
        None
    );
    assert!(
        no_possible_centers
            .atoms
            .iter()
            .all(|atom| atom.prop("_CIPRank").is_none())
    );
    assert_eq!(
        possible_centers.atoms[0].prop("_ChiralityPossible"),
        Some("1")
    );
    assert!(
        possible_centers
            .atoms
            .iter()
            .all(|atom| atom.prop("_CIPRank").is_some())
    );

    let clean_without_possible_centers = assign_with_flags(input.clone(), true, false).unwrap();
    assert_eq!(
        clean_without_possible_centers.atoms[0].prop("_ChiralityPossible"),
        None
    );
    let clean_with_possible_centers = assign_fixed_clean(input.clone());
    assert_eq!(
        clean_with_possible_centers,
        assign_with_flags(input.clone(), true, true).unwrap()
    );
    assert_eq!(
        clean_with_possible_centers.atoms[0].prop("_ChiralityPossible"),
        Some("1")
    );
    assert_eq!(
        input, before,
        "detached assignment must leave its input intact"
    );
}

#[test]
fn cx_legacy_clean_flag_controls_stereo_property_cleanup() {
    let input = raw_large_ring_cis(true);
    let before = input.clone();

    let clean_false = assign_with_flags(input.clone(), false, false).unwrap();
    assert_eq!(clean_false.atoms[0].prop("_CIPCode"), Some("S"));
    assert_eq!(clean_false.atoms[0].prop("_ChiralityPossible"), Some("0"));
    assert_eq!(clean_false.atoms[0].prop("_ringStereochemCand"), Some("1"));
    assert_eq!(clean_false.atoms[0].prop("_ringStereoAtoms"), Some("4,7"));
    assert_eq!(clean_false.bonds[5].prop("_CIPCode"), Some("Z"));

    let clean_true = assign_with_flags(input.clone(), true, false).unwrap();
    assert_eq!(clean_true.atoms[0].prop("_CIPCode"), None);
    assert_eq!(clean_true.atoms[0].prop("_ChiralityPossible"), None);
    assert_eq!(clean_true.atoms[0].prop("_ringStereochemCand"), None);
    assert_eq!(clean_true.atoms[0].prop("_ringStereoAtoms"), None);
    assert_eq!(clean_true.bonds[5].prop("_CIPCode"), None);
    assert_eq!(clean_true.atoms[0].prop("user_marker"), Some("keep"));
    assert_eq!(clean_true.bonds[5].prop("user_marker"), Some("keep"));
    assert_eq!(
        input, before,
        "detached assignment must leave its input intact"
    );
}

#[test]
fn cx_legacy_profile_propagates_topology_validation_errors() {
    let valid = topology_from_specs(
        vec![AtomSpec::new(Element::C), AtomSpec::new(Element::C)],
        Vec::new(),
    );
    let valence = assign_valence_for_topology(&valid, ValenceModel::RdkitLike).unwrap();
    let rings = symmetrized_sssr(&valid, &RingSearchParams::default()).unwrap();
    let mut invalid = valid;
    invalid.atoms.pop();
    let before = invalid.clone();

    let error =
        assign_legacy_stereochemistry_with_flags(invalid.clone(), &valence, &rings, true, false)
            .unwrap_err();
    assert!(matches!(error, LegacyStereoError::InvalidTopology(_)));
    assert_eq!(
        invalid, before,
        "failed assignment must leave its input intact"
    );
}
