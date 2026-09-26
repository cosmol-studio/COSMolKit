use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, CoordinateBlock, TopologyBlock,
};
use cosmolkit_search::{
    SearchTarget, SmartsParseParams, SubstructMatchParams, get_substruct_matches_with_params,
    parse_smarts,
};
use cosmolkit_types::{BondOrder, Element};

fn topology(elements: &[Element], edges: &[(usize, usize)]) -> TopologyBlock {
    let atoms = elements
        .iter()
        .copied()
        .enumerate()
        .map(|(index, element)| Atom::from_spec(AtomId::new(index), AtomSpec::new(element)))
        .collect();
    let bonds = edges
        .iter()
        .enumerate()
        .map(|(index, &(begin, end))| {
            Bond::from_spec(
                BondId::new(index),
                BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
            )
        })
        .collect();
    TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
        .expect("fixed generic-group target topology is valid")
}

fn topology_with_bond_orders(
    elements: &[Element],
    edges: &[(usize, usize)],
    bond_orders: &[BondOrder],
) -> TopologyBlock {
    assert_eq!(edges.len(), bond_orders.len());
    let atoms = elements
        .iter()
        .copied()
        .enumerate()
        .map(|(index, element)| Atom::from_spec(AtomId::new(index), AtomSpec::new(element)))
        .collect();
    let bonds = edges
        .iter()
        .zip(bond_orders)
        .enumerate()
        .map(|(index, (&(begin, end), &order))| {
            Bond::from_spec(
                BondId::new(index),
                BondSpec::new(AtomId::new(begin), AtomId::new(end), order),
            )
        })
        .collect();
    TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
        .expect("fixed ordered generic-group target topology is valid")
}

fn aromatic_carbon_pair(aromatic_atom: bool, aromatic_bond: bool) -> TopologyBlock {
    let atoms = [0, 1]
        .into_iter()
        .map(|index| {
            Atom::from_spec(
                AtomId::new(index),
                AtomSpec::new(Element::C).with_aromatic(aromatic_atom && index == 0),
            )
        })
        .collect();
    let bonds = [Bond::from_spec(
        BondId::new(0),
        BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single)
            .with_aromatic(aromatic_bond),
    )]
    .into_iter()
    .collect();
    TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
        .expect("fixed aromatic generic-group target topology is valid")
}

fn alkenyl_carbon_pair(aromatic_atom: bool, aromatic_bond: bool) -> TopologyBlock {
    let atoms = [0, 1]
        .into_iter()
        .map(|index| {
            Atom::from_spec(
                AtomId::new(index),
                AtomSpec::new(Element::C).with_aromatic(aromatic_atom && index == 0),
            )
        })
        .collect();
    let bonds = [Bond::from_spec(
        BondId::new(0),
        BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Double)
            .with_aromatic(aromatic_bond),
    )]
    .into_iter()
    .collect();
    TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
        .expect("fixed alkenyl generic-group target topology is valid")
}

fn generic_label_matches(
    smarts: &str,
    generic_atom_index: usize,
    label: &str,
    topology: &TopologyBlock,
    enabled: bool,
) -> bool {
    let mut query = parse_smarts(smarts, &SmartsParseParams::default()).expect("parse query");
    query
        .atoms_mut()
        .get_mut(generic_atom_index)
        .expect("generic label query atom exists")
        .set_prop("_QueryAtomGenericLabel", label)
        .expect("set generic group label");
    let coordinates = CoordinateBlock::default();
    let target = SearchTarget::new(topology, &coordinates, &topology.stereo_groups, None, None);
    let params = SubstructMatchParams {
        use_generic_matchers: enabled,
        ..SubstructMatchParams::default()
    };
    !get_substruct_matches_with_params(&target, &query, &params).is_empty()
}

#[test]
fn q42_generic_group_dispatch_keeps_all_source_alias_pairs() {
    // Each long name and short label below selects the same function in the
    // pinned GenericGroups.h::genericMatchers table.
    let aliases = [
        ("Group", "G"),
        ("GroupH", "GH"),
        ("Group*", "G*"),
        ("GroupH*", "GH*"),
        ("Alkyl", "ALK"),
        ("AlkylH", "ALH"),
        ("Alkenyl", "AEL"),
        ("AlkenylH", "AEH"),
        ("Alkynyl", "AYL"),
        ("AlkynylH", "AYH"),
        ("Carbocyclic", "CBC"),
        ("CarbocyclicH", "CBH"),
        ("Carbocycloalkyl", "CAL"),
        ("CarbocycloalkylH", "CAH"),
        ("Carbocycloalkenyl", "CEL"),
        ("CarbocycloalkenylH", "CEH"),
        ("Carboaryl", "ARY"),
        ("CarboarylH", "ARH"),
        ("Cyclic", "CYC"),
        ("CyclicH", "CYH"),
        ("Acyclic", "ACY"),
        ("AcyclicH", "ACH"),
        ("Carboacyclic", "ABC"),
        ("CarboacyclicH", "ABH"),
        ("Heteroacyclic", "AHC"),
        ("HeteroacyclicH", "AHH"),
        ("Alkoxy", "AOX"),
        ("AlkoxyH", "AOH"),
        ("Heterocyclic", "CHC"),
        ("HeterocyclicH", "CHH"),
        ("Heteroaryl", "HAR"),
        ("HeteroarylH", "HAH"),
        ("NoCarbonRing", "CXX"),
        ("NoCarbonRingH", "CXH"),
    ];
    let carbon_pair = topology(&[Element::C, Element::C], &[(0, 1)]);
    let carbon_hydrogen = topology(&[Element::C, Element::H], &[(0, 1)]);
    let carbon_oxygen = topology(&[Element::C, Element::O], &[(0, 1)]);
    let oxygen_carbon = topology(&[Element::O, Element::C], &[(0, 1)]);
    let carbon_chain = topology(&[Element::C, Element::C, Element::C], &[(0, 1), (1, 2)]);
    let cyclopropane = topology(
        &[Element::C, Element::C, Element::C],
        &[(0, 1), (1, 2), (2, 0)],
    );
    let cases = [
        ("C*", 1, &carbon_pair),
        ("C*", 1, &carbon_hydrogen),
        ("C*", 1, &carbon_oxygen),
        ("[O]C", 1, &oxygen_carbon),
        ("C*", 1, &carbon_chain),
        ("C*", 1, &cyclopropane),
    ];

    for (long_name, short_label) in aliases {
        for (smarts, atom_index, target) in cases {
            assert_eq!(
                generic_label_matches(smarts, atom_index, long_name, target, true),
                generic_label_matches(smarts, atom_index, short_label, target, true),
                "source alias {short_label} must select {long_name} for {smarts}"
            );
        }
    }
}

#[test]
fn q42_generic_group_dispatch_preserves_source_h_ignore_ring_and_label_rules() {
    let carbon_hydrogen = topology(&[Element::C, Element::H], &[(0, 1)]);
    let carbon_oxygen = topology(&[Element::C, Element::O], &[(0, 1)]);
    let oxygen_carbon = topology(&[Element::O, Element::C], &[(0, 1)]);
    let carbon_chain = topology(&[Element::C, Element::C, Element::C], &[(0, 1), (1, 2)]);
    let cyclopropane = topology(
        &[Element::C, Element::C, Element::C],
        &[(0, 1), (1, 2), (2, 0)],
    );

    // GenericGroups.cpp::IsHydrogen permits the H-only alternatives when the
    // mapped H has degree one; the non-H Group matcher rejects that sidechain.
    assert!(!generic_label_matches(
        "C*",
        1,
        "Group",
        &carbon_hydrogen,
        true
    ));
    assert!(!generic_label_matches("C*", 1, "G", &carbon_hydrogen, true));
    assert!(generic_label_matches(
        "C*",
        1,
        "GroupH",
        &carbon_hydrogen,
        true
    ));
    assert!(generic_label_matches("C*", 1, "GH", &carbon_hydrogen, true));
    assert!(generic_label_matches(
        "C*",
        1,
        "GroupH*",
        &carbon_hydrogen,
        true
    ));
    assert!(generic_label_matches(
        "C*",
        1,
        "GH*",
        &carbon_hydrogen,
        true
    ));

    // The star family requires a ring bond in the unmatched sidechain.
    assert!(generic_label_matches(
        "C*",
        1,
        "Group*",
        &cyclopropane,
        true
    ));
    assert!(generic_label_matches("C*", 1, "G*", &cyclopropane, true));
    assert!(!generic_label_matches(
        "C*",
        1,
        "Group*",
        &carbon_chain,
        true
    ));

    // The mapped oxygen is ignored by the alkyl walk. Without the source ignore
    // bitmap it would be visited and rejected by AlkylAtomMatcher.
    assert!(generic_label_matches(
        "[O]C",
        1,
        "Alkyl",
        &oxygen_carbon,
        true
    ));

    // A known label rejects O, but unknown labels are ignored by the source
    // dispatcher. Disabling generic matchers also preserves the plain query.
    assert!(!generic_label_matches(
        "C*",
        1,
        "Alkyl",
        &carbon_oxygen,
        true
    ));
    assert!(generic_label_matches(
        "C*",
        1,
        "Alkyl",
        &carbon_oxygen,
        false
    ));
    assert!(generic_label_matches(
        "C*",
        1,
        "not-a-generic-label",
        &carbon_oxygen,
        true
    ));

    // GenericGroups.cpp skips query atoms of degree > 1 before reading labels.
    // If this center were dispatched to GroupStar, its matched neighbors would
    // leave no ring bond and the otherwise valid chain match would be rejected.
    assert!(generic_label_matches(
        "CCC",
        1,
        "Group*",
        &carbon_chain,
        true
    ));
}

#[test]
fn q64_generic_dispatch_checks_each_single_atom_label_and_ignores_unknown_labels() {
    let target_topology = topology(&[Element::C, Element::O], &[]);
    let mut query = parse_smarts("C.O", &SmartsParseParams::default())
        .expect("fixed disconnected Q64 query parses");
    query
        .atoms_mut()
        .get_mut(0)
        .expect("first Q64 query atom exists")
        .set_prop("_QueryAtomGenericLabel", "not-a-generic-label")
        .expect("set first Q64 generic label");
    query
        .atoms_mut()
        .get_mut(1)
        .expect("second Q64 query atom exists")
        .set_prop("_QueryAtomGenericLabel", "also-not-a-generic-label")
        .expect("set second Q64 generic label");

    let coordinates = CoordinateBlock::default();
    let target = SearchTarget::new(
        &target_topology,
        &coordinates,
        &target_topology.stereo_groups,
        None,
        None,
    );
    let generic_params = SubstructMatchParams {
        use_generic_matchers: true,
        ..SubstructMatchParams::default()
    };
    assert!(
        !get_substruct_matches_with_params(&target, &query, &generic_params).is_empty(),
        "unknown source labels are ignored for every single-atom query"
    );

    // The first unknown label is ignored, but dispatch continues to the later
    // known Alkyl label, whose source matcher rejects the mapped oxygen.
    query
        .atoms_mut()
        .get_mut(1)
        .expect("second Q64 query atom exists")
        .set_prop("_QueryAtomGenericLabel", "Alkyl")
        .expect("set later known Q64 generic label");
    assert!(get_substruct_matches_with_params(&target, &query, &generic_params).is_empty());
}

#[test]
fn q44_all_atoms_match_traverses_the_unignored_component_and_combines_witnesses() {
    let carbon_pair = topology(&[Element::C, Element::C], &[(0, 1)]);
    let carbon_hydrogen = topology(&[Element::C, Element::H], &[(0, 1)]);
    let oxygen_carbon = topology(&[Element::O, Element::C], &[(0, 1)]);
    let extended_chain = topology(
        &[Element::O, Element::C, Element::C, Element::C],
        &[(0, 1), (1, 2), (2, 3)],
    );
    let disconnected_oxygen = topology(
        &[Element::O, Element::C, Element::C, Element::O],
        &[(0, 1), (1, 2)],
    );
    let oxygen_branch = topology(
        &[Element::O, Element::C, Element::C, Element::O],
        &[(0, 1), (1, 2), (2, 3)],
    );
    let ring_branch = topology(
        &[Element::O, Element::C, Element::C, Element::C, Element::C],
        &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 2)],
    );
    let cyclopropane = topology(
        &[Element::C, Element::C, Element::C],
        &[(0, 1), (1, 2), (2, 0)],
    );
    let carbon_chain = topology(&[Element::C, Element::C, Element::C], &[(0, 1), (1, 2)]);
    let alkene = topology_with_bond_orders(
        &[Element::O, Element::C, Element::C, Element::C],
        &[(0, 1), (1, 2), (2, 3)],
        &[BondOrder::Single, BondOrder::Single, BondOrder::Double],
    );
    let saturated_chain = topology_with_bond_orders(
        &[Element::O, Element::C, Element::C, Element::C],
        &[(0, 1), (1, 2), (2, 3)],
        &[BondOrder::Single, BondOrder::Single, BondOrder::Single],
    );

    // AllAtomsMatch tests the root for a required atom witness, while its H
    // wrapper takes the source IsHydrogen shortcut before the component walk.
    assert!(generic_label_matches("C*", 1, "Group", &carbon_pair, true));
    assert!(!generic_label_matches(
        "C*",
        1,
        "Group",
        &carbon_hydrogen,
        true
    ));
    assert!(generic_label_matches(
        "C*",
        1,
        "GroupH",
        &carbon_hydrogen,
        true
    ));

    // The mapped O is preset in AllAtomsMatch's by-value ignore mask. The
    // traversal follows only the connected unignored sidechain and rejects an
    // unignored O reached from the generic root.
    assert!(generic_label_matches(
        "[O]C",
        1,
        "Alkyl",
        &oxygen_carbon,
        true
    ));
    assert!(generic_label_matches(
        "[O]C*",
        2,
        "Alkyl",
        &extended_chain,
        true
    ));
    assert!(generic_label_matches(
        "[O]C*",
        2,
        "Alkyl",
        &disconnected_oxygen,
        true
    ));
    assert!(!generic_label_matches(
        "[O]C*",
        2,
        "Alkyl",
        &oxygen_branch,
        true
    ));

    // Alkyl's bond callback rejects ring bonds; GroupStar's separate required
    // ring-bond witness makes the final atom-and-bond conjunction true only in
    // the ring target.
    assert!(!generic_label_matches(
        "[O]C*",
        2,
        "Alkyl",
        &ring_branch,
        true
    ));
    assert!(generic_label_matches(
        "C*",
        1,
        "Group*",
        &cyclopropane,
        true
    ));
    assert!(!generic_label_matches(
        "C*",
        1,
        "Group*",
        &carbon_chain,
        true
    ));

    // Unsaturated helpers require an allowed bond traversal and at least one
    // source-defined double-bond witness.
    assert!(generic_label_matches("[O]C*", 2, "Alkenyl", &alkene, true));
    assert!(!generic_label_matches(
        "[O]C*",
        2,
        "Alkenyl",
        &saturated_chain,
        true
    ));

    // Known labels reject failing components; unknown labels and disabled
    // generic matching preserve the ordinary query result.
    assert!(!generic_label_matches(
        "[O]C*",
        2,
        "Alkyl",
        &oxygen_branch,
        true
    ));
    assert!(generic_label_matches(
        "[O]C*",
        2,
        "unknown-label",
        &oxygen_branch,
        true
    ));
    assert!(generic_label_matches(
        "[O]C*",
        2,
        "Alkyl",
        &oxygen_branch,
        false
    ));
}

#[test]
fn q45_group_variants_keep_source_hydrogen_and_ring_fallbacks() {
    let isolated_hydrogen = topology(&[Element::H], &[]);
    let isolated_carbon = topology(&[Element::C], &[]);
    let carbon_hydrogen = topology(&[Element::C, Element::H], &[(0, 1)]);
    let carbon_pair = topology(&[Element::C, Element::C], &[(0, 1)]);
    let carbon_chain = topology(&[Element::C, Element::C, Element::C], &[(0, 1), (1, 2)]);
    let cyclopropane = topology(
        &[Element::C, Element::C, Element::C],
        &[(0, 1), (1, 2), (2, 0)],
    );

    // Group requires a non-H witness. Its H variant accepts only a degree-one
    // H shortcut; an isolated H falls through and still fails.
    assert!(generic_label_matches(
        "*",
        0,
        "Group",
        &isolated_carbon,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "Group",
        &isolated_hydrogen,
        true
    ));
    assert!(generic_label_matches(
        "*",
        0,
        "GroupH",
        &isolated_carbon,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "GroupH",
        &isolated_hydrogen,
        true
    ));
    assert!(!generic_label_matches(
        "C*",
        1,
        "Group",
        &carbon_hydrogen,
        true
    ));
    assert!(generic_label_matches(
        "C*",
        1,
        "GroupH",
        &carbon_hydrogen,
        true
    ));

    // The star variants additionally require a ring bond. Their H variant
    // bypasses that requirement only for source IsHydrogen's degree-one case.
    assert!(!generic_label_matches(
        "*",
        0,
        "Group*",
        &isolated_carbon,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "GroupH*",
        &isolated_hydrogen,
        true
    ));
    assert!(generic_label_matches(
        "C*",
        1,
        "Group*",
        &cyclopropane,
        true
    ));
    assert!(generic_label_matches(
        "C*",
        1,
        "GroupH*",
        &cyclopropane,
        true
    ));
    assert!(!generic_label_matches(
        "C*",
        1,
        "Group*",
        &carbon_chain,
        true
    ));
    assert!(!generic_label_matches(
        "C*",
        1,
        "GroupH*",
        &carbon_chain,
        true
    ));
    assert!(generic_label_matches("C*", 1, "GroupH", &carbon_pair, true));

    // Unrecognized source labels are ignored, and disabling generic matching
    // leaves the ordinary wildcard query intact.
    assert!(generic_label_matches(
        "*",
        0,
        "unknown-label",
        &isolated_hydrogen,
        true
    ));
    assert!(generic_label_matches(
        "*",
        0,
        "Group*",
        &isolated_carbon,
        false
    ));
}

#[test]
fn q46_alkyl_variants_keep_source_component_carbon_and_bond_rules() {
    let carbon_chain = topology(&[Element::C, Element::C, Element::C], &[(0, 1), (1, 2)]);
    let carbon_hydrogen = topology(&[Element::C, Element::H], &[(0, 1)]);
    let isolated_hydrogen = topology(&[Element::H], &[]);
    let carbon_oxygen = topology(&[Element::C, Element::O], &[(0, 1)]);
    let disconnected_oxygen = topology(&[Element::C, Element::O], &[]);
    let double_bond =
        topology_with_bond_orders(&[Element::C, Element::C], &[(0, 1)], &[BondOrder::Double]);
    let cyclopropane = topology(
        &[Element::C, Element::C, Element::C],
        &[(0, 1), (1, 2), (2, 0)],
    );
    let aromatic_atom = aromatic_carbon_pair(true, false);
    let aromatic_bond = aromatic_carbon_pair(false, true);
    let oxygen_carbon_chain = topology(
        &[Element::O, Element::C, Element::C, Element::C],
        &[(0, 1), (1, 2), (2, 3)],
    );
    let oxygen_branch = topology(
        &[Element::O, Element::C, Element::C, Element::O],
        &[(0, 1), (1, 2), (2, 3)],
    );

    // Both variants accept nonaromatic C/H single-bond components. Alkyl
    // additionally requires a carbon witness; AlkylH permits an H-only root.
    assert!(generic_label_matches("*", 0, "Alkyl", &carbon_chain, true));
    assert!(generic_label_matches("*", 0, "ALK", &carbon_chain, true));
    assert!(generic_label_matches("*", 0, "AlkylH", &carbon_chain, true));
    assert!(generic_label_matches(
        "*",
        0,
        "Alkyl",
        &carbon_hydrogen,
        true
    ));
    assert!(generic_label_matches(
        "*",
        0,
        "AlkylH",
        &carbon_hydrogen,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "Alkyl",
        &isolated_hydrogen,
        true
    ));
    assert!(generic_label_matches(
        "*",
        0,
        "AlkylH",
        &isolated_hydrogen,
        true
    ));

    // Traversal covers only the root's connected component and starts with
    // query-mapped atoms ignored; an unignored connected O fails the atom rule.
    assert!(generic_label_matches(
        "*",
        0,
        "Alkyl",
        &disconnected_oxygen,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "Alkyl",
        &carbon_oxygen,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "AlkylH",
        &carbon_oxygen,
        true
    ));
    assert!(generic_label_matches(
        "[O]C*",
        2,
        "Alkyl",
        &oxygen_carbon_chain,
        true
    ));
    assert!(!generic_label_matches(
        "[O]C*",
        2,
        "Alkyl",
        &oxygen_branch,
        true
    ));

    // Every traversed bond must be single, nonaromatic, and outside rings;
    // atom aromaticity is independently rejected.
    for invalid in [&double_bond, &cyclopropane, &aromatic_atom, &aromatic_bond] {
        assert!(!generic_label_matches("*", 0, "Alkyl", invalid, true));
        assert!(!generic_label_matches("*", 0, "AlkylH", invalid, true));
    }
}

#[test]
fn q47_alkenyl_variants_require_source_double_bond_witness() {
    let alkene =
        topology_with_bond_orders(&[Element::C, Element::C], &[(0, 1)], &[BondOrder::Double]);
    let alkane = topology(&[Element::C, Element::C], &[(0, 1)]);
    let alkyne =
        topology_with_bond_orders(&[Element::C, Element::C], &[(0, 1)], &[BondOrder::Triple]);
    let ring_alkene = topology_with_bond_orders(
        &[Element::C, Element::C, Element::C],
        &[(0, 1), (1, 2), (2, 0)],
        &[BondOrder::Double, BondOrder::Single, BondOrder::Single],
    );
    let disconnected_oxygen = topology_with_bond_orders(
        &[Element::C, Element::C, Element::O],
        &[(0, 1)],
        &[BondOrder::Double],
    );
    let oxygen_branch = topology_with_bond_orders(
        &[Element::C, Element::C, Element::C, Element::O],
        &[(0, 1), (1, 2), (1, 3)],
        &[BondOrder::Single, BondOrder::Double, BondOrder::Single],
    );
    let mapped_alkene = topology_with_bond_orders(
        &[Element::C, Element::C, Element::C, Element::C],
        &[(0, 1), (1, 2), (2, 3)],
        &[BondOrder::Single, BondOrder::Single, BondOrder::Double],
    );
    let carbon_hydrogen = topology(&[Element::C, Element::H], &[(0, 1)]);
    let isolated_hydrogen = topology(&[Element::H], &[]);
    let aromatic_atom = alkenyl_carbon_pair(true, false);
    let aromatic_bond = alkenyl_carbon_pair(false, true);

    // Acyclic nonaromatic double bonds are the required witness. Single or
    // triple bonds cannot substitute for it.
    assert!(generic_label_matches("*", 0, "Alkenyl", &alkene, true));
    assert!(generic_label_matches("*", 0, "AEL", &alkene, true));
    assert!(generic_label_matches("*", 0, "AlkenylH", &alkene, true));
    for invalid in [&alkane, &alkyne] {
        assert!(!generic_label_matches("*", 0, "Alkenyl", invalid, true));
        assert!(!generic_label_matches("*", 0, "AlkenylH", invalid, true));
    }

    // Every component bond must remain outside rings and every visited atom
    // must be nonaromatic C/H, even when a double-bond witness is present.
    for invalid in [&ring_alkene, &oxygen_branch, &aromatic_atom, &aromatic_bond] {
        assert!(!generic_label_matches("*", 0, "Alkenyl", invalid, true));
        assert!(!generic_label_matches("*", 0, "AlkenylH", invalid, true));
    }
    assert!(generic_label_matches(
        "*",
        0,
        "Alkenyl",
        &disconnected_oxygen,
        true
    ));

    // The mapped C atoms are ignored before the generic walk, leaving only
    // the unmatched C=C component to satisfy the source bond witness.
    assert!(generic_label_matches(
        "[C]C*",
        2,
        "Alkenyl",
        &mapped_alkene,
        true
    ));

    // AlkenylH adds the degree-one IsHydrogen shortcut; an isolated H has
    // degree zero and does not take that branch.
    assert!(!generic_label_matches(
        "C*",
        1,
        "Alkenyl",
        &carbon_hydrogen,
        true
    ));
    assert!(generic_label_matches(
        "C*",
        1,
        "AlkenylH",
        &carbon_hydrogen,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "AlkenylH",
        &isolated_hydrogen,
        true
    ));
}

#[test]
fn q48_alkynyl_variants_require_source_triple_bond_witness() {
    let alkyne =
        topology_with_bond_orders(&[Element::C, Element::C], &[(0, 1)], &[BondOrder::Triple]);
    let alkene =
        topology_with_bond_orders(&[Element::C, Element::C], &[(0, 1)], &[BondOrder::Double]);
    let alkane = topology(&[Element::C, Element::C], &[(0, 1)]);
    let ring_alkyne = topology_with_bond_orders(
        &[Element::C, Element::C, Element::C],
        &[(0, 1), (1, 2), (2, 0)],
        &[BondOrder::Triple, BondOrder::Single, BondOrder::Single],
    );
    let disconnected_oxygen = topology_with_bond_orders(
        &[Element::C, Element::C, Element::O],
        &[(0, 1)],
        &[BondOrder::Triple],
    );
    let oxygen_branch = topology_with_bond_orders(
        &[Element::C, Element::C, Element::C, Element::O],
        &[(0, 1), (1, 2), (1, 3)],
        &[BondOrder::Single, BondOrder::Triple, BondOrder::Single],
    );
    let mapped_alkyne = topology_with_bond_orders(
        &[Element::C, Element::C, Element::C, Element::C],
        &[(0, 1), (1, 2), (2, 3)],
        &[BondOrder::Single, BondOrder::Single, BondOrder::Triple],
    );
    let carbon_hydrogen = topology(&[Element::C, Element::H], &[(0, 1)]);
    let isolated_hydrogen = topology(&[Element::H], &[]);

    // Only an acyclic triple bond witnesses Alkynyl. Double and single bonds
    // remain permitted only when a triple bond is also present.
    assert!(generic_label_matches("*", 0, "Alkynyl", &alkyne, true));
    assert!(generic_label_matches("*", 0, "AYL", &alkyne, true));
    assert!(generic_label_matches("*", 0, "AlkynylH", &alkyne, true));
    for invalid in [&alkene, &alkane] {
        assert!(!generic_label_matches("*", 0, "Alkynyl", invalid, true));
        assert!(!generic_label_matches("*", 0, "AlkynylH", invalid, true));
    }

    // Ring bonds and connected non-C/H atoms reject the component, while a
    // disconnected O stays outside the traversal. Mapped atoms are ignored.
    for invalid in [&ring_alkyne, &oxygen_branch] {
        assert!(!generic_label_matches("*", 0, "Alkynyl", invalid, true));
        assert!(!generic_label_matches("*", 0, "AlkynylH", invalid, true));
    }
    assert!(generic_label_matches(
        "*",
        0,
        "Alkynyl",
        &disconnected_oxygen,
        true
    ));
    assert!(generic_label_matches(
        "[C]C*",
        2,
        "Alkynyl",
        &mapped_alkyne,
        true
    ));

    // The H variant uses IsHydrogen's degree-one shortcut, but an isolated H
    // falls through and still lacks the required triple-bond witness.
    assert!(!generic_label_matches(
        "C*",
        1,
        "Alkynyl",
        &carbon_hydrogen,
        true
    ));
    assert!(generic_label_matches(
        "C*",
        1,
        "AlkynylH",
        &carbon_hydrogen,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "AlkynylH",
        &isolated_hydrogen,
        true
    ));
}

#[test]
fn q49_acyclic_variants_use_source_ring_and_non_h_witnesses() {
    let carbon_chain = topology(&[Element::C, Element::C, Element::C], &[(0, 1), (1, 2)]);
    let cyclopropane = topology(
        &[Element::C, Element::C, Element::C],
        &[(0, 1), (1, 2), (2, 0)],
    );
    let ring_tail = topology(
        &[Element::C, Element::C, Element::C, Element::C, Element::C],
        &[(0, 1), (1, 2), (2, 0), (2, 3), (3, 4)],
    );
    let disconnected_ring = topology(
        &[Element::C, Element::C, Element::C, Element::C],
        &[(1, 2), (2, 3), (3, 1)],
    );
    let carbon_hydrogen = topology(&[Element::C, Element::H], &[(0, 1)]);
    let isolated_hydrogen = topology(&[Element::H], &[]);

    // Both variants accept acyclic connected components. Acyclic additionally
    // requires at least one non-H atom, while AcyclicH has no such witness.
    assert!(generic_label_matches(
        "*",
        0,
        "Acyclic",
        &carbon_chain,
        true
    ));
    assert!(generic_label_matches("*", 0, "ACY", &carbon_chain, true));
    assert!(generic_label_matches(
        "*",
        0,
        "AcyclicH",
        &carbon_chain,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "Acyclic",
        &isolated_hydrogen,
        true
    ));
    assert!(generic_label_matches(
        "*",
        0,
        "AcyclicH",
        &isolated_hydrogen,
        true
    ));
    assert!(generic_label_matches(
        "C*",
        1,
        "ACH",
        &carbon_hydrogen,
        true
    ));

    // Any ring atom in the traversed component rejects both variants, but a
    // disconnected ring does not alter the isolated carbon component.
    assert!(!generic_label_matches(
        "*",
        0,
        "Acyclic",
        &cyclopropane,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "AcyclicH",
        &cyclopropane,
        true
    ));
    assert!(generic_label_matches(
        "*",
        0,
        "Acyclic",
        &disconnected_ring,
        true
    ));
    assert!(generic_label_matches(
        "*",
        0,
        "AcyclicH",
        &disconnected_ring,
        true
    ));

    // Query-mapped ring atoms are preset in the source ignore bitmap, so the
    // acyclic matcher sees only the terminal generic atom's remaining component.
    assert!(!generic_label_matches("*", 0, "Acyclic", &ring_tail, true));
    assert!(generic_label_matches(
        "C1CC1C*", 4, "Acyclic", &ring_tail, true
    ));
    assert!(generic_label_matches(
        "C1CC1C*", 4, "AcyclicH", &ring_tail, true
    ));
}

#[test]
fn q50_carboacyclic_variants_keep_source_c_h_ring_rules() {
    let carbon_chain = topology(&[Element::C, Element::C, Element::C], &[(0, 1), (1, 2)]);
    let carbon_alkene =
        topology_with_bond_orders(&[Element::C, Element::C], &[(0, 1)], &[BondOrder::Double]);
    let carbon_hydrogen = topology(&[Element::C, Element::H], &[(0, 1)]);
    let isolated_hydrogen = topology(&[Element::H], &[]);
    let cyclopropane = topology(
        &[Element::C, Element::C, Element::C],
        &[(0, 1), (1, 2), (2, 0)],
    );
    let disconnected_ring = topology(
        &[Element::C, Element::C, Element::C, Element::C],
        &[(1, 2), (2, 3), (3, 1)],
    );
    let oxygen_branch = topology(&[Element::C, Element::C, Element::O], &[(0, 1), (1, 2)]);
    let oxygen_mapped_chain = topology(
        &[Element::O, Element::C, Element::C, Element::C],
        &[(0, 1), (1, 2), (2, 3)],
    );
    let aromatic_carbon = aromatic_carbon_pair(true, false);

    // The source accepts acyclic C/H components. Only Carboacyclic requires
    // a carbon witness; neither source helper constrains bond order or atom
    // aromaticity beyond atomic number and ring membership.
    assert!(generic_label_matches(
        "*",
        0,
        "Carboacyclic",
        &carbon_chain,
        true
    ));
    assert!(generic_label_matches("*", 0, "ABC", &carbon_chain, true));
    assert!(generic_label_matches(
        "*",
        0,
        "CarboacyclicH",
        &carbon_alkene,
        true
    ));
    assert!(generic_label_matches(
        "*",
        0,
        "Carboacyclic",
        &aromatic_carbon,
        true
    ));
    assert!(generic_label_matches(
        "*",
        0,
        "CarboacyclicH",
        &aromatic_carbon,
        true
    ));
    assert!(!generic_label_matches(
        "C*",
        1,
        "Carboacyclic",
        &carbon_hydrogen,
        true
    ));
    assert!(generic_label_matches(
        "C*",
        1,
        "CarboacyclicH",
        &carbon_hydrogen,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "Carboacyclic",
        &isolated_hydrogen,
        true
    ));
    assert!(generic_label_matches(
        "*",
        0,
        "ABH",
        &isolated_hydrogen,
        true
    ));

    // O and ring atoms are rejected when the traversal reaches them; a
    // disconnected ring is outside the root component, and query-mapped O is
    // preset in the ignore bitmap.
    for invalid in [&cyclopropane, &oxygen_branch] {
        assert!(!generic_label_matches(
            "*",
            0,
            "Carboacyclic",
            invalid,
            true
        ));
        assert!(!generic_label_matches(
            "*",
            0,
            "CarboacyclicH",
            invalid,
            true
        ));
    }
    assert!(generic_label_matches(
        "*",
        0,
        "Carboacyclic",
        &disconnected_ring,
        true
    ));
    assert!(generic_label_matches(
        "[O]C*",
        2,
        "Carboacyclic",
        &oxygen_mapped_chain,
        true
    ));
}

#[test]
fn q51_heteroacyclic_variants_use_source_acyclic_hetero_witnesses() {
    let carbon_oxygen_chain = topology(&[Element::C, Element::O, Element::C], &[(0, 1), (1, 2)]);
    let carbon_oxygen_double =
        topology_with_bond_orders(&[Element::C, Element::O], &[(0, 1)], &[BondOrder::Double]);
    let carbon_chain = topology(&[Element::C, Element::C, Element::C], &[(0, 1), (1, 2)]);
    let oxygen_ring = topology(
        &[Element::C, Element::O, Element::C],
        &[(0, 1), (1, 2), (2, 0)],
    );
    let disconnected_ring = topology(
        &[Element::O, Element::C, Element::C, Element::C],
        &[(1, 2), (2, 3), (3, 1)],
    );
    let oxygen_mapped_chain = topology(&[Element::O, Element::C, Element::C], &[(0, 1), (1, 2)]);
    let carbon_hydrogen = topology(&[Element::C, Element::H], &[(0, 1)]);
    let isolated_hydrogen = topology(&[Element::H], &[]);

    // The source requires an acyclic component with at least one atom whose
    // atomic number is neither carbon nor hydrogen; it has no bond callback.
    assert!(generic_label_matches(
        "*",
        0,
        "Heteroacyclic",
        &carbon_oxygen_chain,
        true
    ));
    assert!(generic_label_matches(
        "C~*",
        1,
        "Heteroacyclic",
        &carbon_oxygen_double,
        true
    ));
    assert!(!generic_label_matches(
        "C*",
        1,
        "Heteroacyclic",
        &carbon_chain,
        true
    ));

    // Ring membership rejects every atom reached in a cyclic component, while
    // a disconnected ring is outside the isolated oxygen root's traversal.
    assert!(!generic_label_matches(
        "*",
        0,
        "Heteroacyclic",
        &oxygen_ring,
        true
    ));
    assert!(generic_label_matches(
        "*",
        0,
        "Heteroacyclic",
        &disconnected_ring,
        true
    ));

    // Query-mapped atoms are preignored by AllAtomsMatch and cannot provide
    // the hetero witness, even though oxygen is connected beyond mapped C.
    assert!(!generic_label_matches(
        "[O]C*",
        2,
        "Heteroacyclic",
        &oxygen_mapped_chain,
        true
    ));
    assert!(!generic_label_matches(
        "[O]C*",
        2,
        "HeteroacyclicH",
        &oxygen_mapped_chain,
        true
    ));

    // HeteroacyclicH adds only the degree-one IsHydrogen shortcut; isolated H
    // falls through to the ordinary hetero-witness requirement and fails.
    assert!(!generic_label_matches(
        "C*",
        1,
        "Heteroacyclic",
        &carbon_hydrogen,
        true
    ));
    assert!(generic_label_matches(
        "C*",
        1,
        "HeteroacyclicH",
        &carbon_hydrogen,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "HeteroacyclicH",
        &isolated_hydrogen,
        true
    ));
}

#[test]
fn q52_alkoxyacyclic_variants_keep_source_attachment_and_component_rules() {
    let sidechain_target =
        |aromatic_start: bool, internal_order: BondOrder, internal_aromatic: bool| {
            let atoms = [
                Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::O)),
                Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::O)),
                Atom::from_spec(
                    AtomId::new(2),
                    AtomSpec::new(Element::C).with_aromatic(aromatic_start),
                ),
                Atom::from_spec(AtomId::new(3), AtomSpec::new(Element::C)),
            ]
            .into_iter()
            .collect();
            let bonds = [
                Bond::from_spec(
                    BondId::new(0),
                    BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
                ),
                Bond::from_spec(
                    BondId::new(1),
                    BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Single),
                ),
                Bond::from_spec(
                    BondId::new(2),
                    BondSpec::new(AtomId::new(2), AtomId::new(3), internal_order)
                        .with_aromatic(internal_aromatic),
                ),
            ]
            .into_iter()
            .collect();
            TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
                .expect("fixed alkoxy sidechain topology is valid")
        };

    let single_sidechain = sidechain_target(false, BondOrder::Single, false);
    let aromatic_start = sidechain_target(true, BondOrder::Single, false);
    let double_sidechain = sidechain_target(false, BondOrder::Double, false);
    let aromatic_sidechain_bond = sidechain_target(false, BondOrder::Single, true);
    let oxygen_carbon_hydrogen = topology(&[Element::C, Element::O, Element::H], &[(0, 1), (1, 2)]);
    let carbon_oxygen_carbon = topology(&[Element::C, Element::O, Element::C], &[(0, 1), (1, 2)]);
    let oxygen_oxygen_carbon_ring = topology(
        &[Element::O, Element::O, Element::C, Element::C, Element::C],
        &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 2)],
    );
    let disconnected_ring = topology(
        &[
            Element::O,
            Element::O,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
        ],
        &[(0, 1), (1, 2), (2, 3), (4, 5), (5, 6), (6, 4)],
    );
    let oxygen_carbon = topology(&[Element::O, Element::C], &[(0, 1)]);
    let isolated_oxygen = topology(&[Element::O], &[]);
    let degree_two_carbon = topology(&[Element::C, Element::C, Element::C], &[(0, 1), (1, 2)]);
    let carbon_hydrogen = topology(&[Element::C, Element::H], &[(0, 1)]);
    let isolated_hydrogen = topology(&[Element::H], &[]);

    // Query-mapped O/O are skipped when selecting the first sidechain neighbor;
    // the remaining C/H component needs a carbon witness and allowed bonds.
    assert!(generic_label_matches(
        "OO",
        1,
        "Alkoxy",
        &single_sidechain,
        true
    ));
    assert!(generic_label_matches(
        "OO",
        1,
        "AOX",
        &single_sidechain,
        true
    ));
    assert!(!generic_label_matches(
        "CO",
        1,
        "Alkoxy",
        &oxygen_carbon_hydrogen,
        true
    ));
    assert!(generic_label_matches(
        "CO",
        1,
        "Alkoxy",
        &carbon_oxygen_carbon,
        true
    ));

    // The sidechain rejects aromatic atoms, non-single/aromatic bonds, and
    // ring bonds; a disconnected ring is outside the unignored component.
    for invalid in [&aromatic_start, &double_sidechain, &aromatic_sidechain_bond] {
        assert!(!generic_label_matches("OO", 1, "Alkoxy", invalid, true));
    }
    assert!(!generic_label_matches(
        "OO",
        1,
        "Alkoxy",
        &oxygen_oxygen_carbon_ring,
        true
    ));
    assert!(generic_label_matches(
        "OO",
        1,
        "Alkoxy",
        &disconnected_ring,
        true
    ));

    // The source dispatcher skips this COC label because its query degree is
    // greater than one, so this does not claim the helper's no-neighbor path.
    assert!(generic_label_matches(
        "COC",
        1,
        "Alkoxy",
        &carbon_oxygen_carbon,
        true
    ));
    assert!(!generic_label_matches(
        "O",
        0,
        "Alkoxy",
        &oxygen_carbon,
        true
    ));
    assert!(!generic_label_matches(
        "O",
        0,
        "Alkoxy",
        &isolated_oxygen,
        true
    ));
    assert!(!generic_label_matches(
        "C",
        0,
        "Alkoxy",
        &degree_two_carbon,
        true
    ));

    // AlkoxyH adds only the source degree-one hydrogen shortcut, then falls
    // back to the same alkoxy matcher for other atoms.
    assert!(!generic_label_matches(
        "C*",
        1,
        "Alkoxy",
        &carbon_hydrogen,
        true
    ));
    assert!(generic_label_matches(
        "C*",
        1,
        "AlkoxyH",
        &carbon_hydrogen,
        true
    ));
    assert!(generic_label_matches(
        "C*",
        1,
        "AOH",
        &carbon_hydrogen,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "AlkoxyH",
        &isolated_hydrogen,
        true
    ));
    assert!(generic_label_matches(
        "OO",
        1,
        "AlkoxyH",
        &single_sidechain,
        true
    ));
}

#[test]
fn q53_ring_atom_checks_preserve_root_ignore_and_fused_ring_membership() {
    let cyclopropane = topology(
        &[Element::C, Element::C, Element::C],
        &[(0, 1), (1, 2), (2, 0)],
    );
    let heterocyclopropane = topology(
        &[Element::C, Element::O, Element::C],
        &[(0, 1), (1, 2), (2, 0)],
    );
    let oxygen_ring = topology(
        &[Element::O, Element::O, Element::O],
        &[(0, 1), (1, 2), (2, 0)],
    );
    let carbon_chain = topology(&[Element::C, Element::C, Element::C], &[(0, 1), (1, 2)]);
    let fused_carbon_rings = topology(
        &[
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
        ],
        &[
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 0),
            (3, 6),
            (6, 7),
            (7, 8),
            (8, 9),
            (9, 2),
        ],
    );
    let fused_hetero_ring = topology(
        &[
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::O,
            Element::C,
            Element::C,
            Element::C,
        ],
        &[
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 0),
            (3, 6),
            (6, 7),
            (7, 8),
            (8, 9),
            (9, 2),
        ],
    );
    let carbon_hydrogen = topology(&[Element::C, Element::H], &[(0, 1)]);

    // checkAtomRing exempts its fused-ring root from the ignore and atom
    // predicates. Other query-mapped ring atoms remain ignored and reject the
    // generic match, even when both atoms are carbon.
    assert!(generic_label_matches(
        "*",
        0,
        "Carbocyclic",
        &cyclopropane,
        true
    ));
    assert!(generic_label_matches("*", 0, "CBC", &cyclopropane, true));
    assert!(!generic_label_matches(
        "C*",
        1,
        "Carbocyclic",
        &cyclopropane,
        true
    ));

    // The atom predicate applies to every atom in each ring; a heteroatom
    // rejects Carbocyclic, and a carbon rejects NoCarbonRing.
    assert!(!generic_label_matches(
        "*",
        0,
        "Carbocyclic",
        &heterocyclopropane,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "NoCarbonRing",
        &heterocyclopropane,
        true
    ));
    assert!(generic_label_matches(
        "*",
        0,
        "NoCarbonRing",
        &oxygen_ring,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "Carbocyclic",
        &carbon_chain,
        true
    ));

    // FusedRingMatch checks each SSSR ring sharing at least two atoms with the
    // accumulated ring system, so the second ring's O also rejects the match.
    assert!(generic_label_matches(
        "*",
        0,
        "Carbocyclic",
        &fused_carbon_rings,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "Carbocyclic",
        &fused_hetero_ring,
        true
    ));

    // The source H variant keeps IsHydrogen's degree-one shortcut before its
    // ring matcher; the ordinary variant has no such shortcut.
    assert!(!generic_label_matches(
        "C*",
        1,
        "Carbocyclic",
        &carbon_hydrogen,
        true
    ));
    assert!(generic_label_matches(
        "C*",
        1,
        "CarbocyclicH",
        &carbon_hydrogen,
        true
    ));
}

#[test]
fn q54_ring_bond_checks_preserve_matcher_and_per_ring_witness_rules() {
    let ring_edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)];
    let carbon_ring = [Element::C; 6];
    let saturated_orders = [BondOrder::Single; 6];
    let cyclohexane = topology_with_bond_orders(&carbon_ring, &ring_edges, &saturated_orders);
    let cyclohexene = topology_with_bond_orders(
        &carbon_ring,
        &ring_edges,
        &[
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
        ],
    );
    let aromatic_order_ring = topology_with_bond_orders(
        &carbon_ring,
        &ring_edges,
        &[
            BondOrder::Aromatic,
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
        ],
    );
    let carbon_chain_with_double = topology_with_bond_orders(
        &[Element::C, Element::C, Element::C],
        &[(0, 1), (1, 2)],
        &[BondOrder::Double, BondOrder::Single],
    );
    let fused_edges = [
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 5),
        (5, 0),
        (3, 6),
        (6, 7),
        (7, 8),
        (8, 9),
        (9, 2),
    ];
    let fused_carbons = [Element::C; 10];
    let fused_saturated_orders = [BondOrder::Single; 11];
    let fused_saturated =
        topology_with_bond_orders(&fused_carbons, &fused_edges, &fused_saturated_orders);
    let fused_one_unsaturated_ring = topology_with_bond_orders(
        &fused_carbons,
        &fused_edges,
        &[
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
        ],
    );
    let fused_two_unsaturated_rings = topology_with_bond_orders(
        &fused_carbons,
        &fused_edges,
        &[
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
        ],
    );
    let carbon_hydrogen = topology(&[Element::C, Element::H], &[(0, 1)]);

    // checkBondRing evaluates its all-bond matcher before the optional
    // per-ring witness. A degree-zero wildcard root still carries its source
    // ignore bit, which does not suppress ring-bond checks.
    assert!(generic_label_matches(
        "*",
        0,
        "Carbocycloalkyl",
        &cyclohexane,
        true
    ));
    assert!(generic_label_matches(
        "*",
        0,
        "Carbocycloalkenyl",
        &cyclohexene,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "Carbocycloalkyl",
        &cyclohexene,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "Carbocycloalkenyl",
        &cyclohexane,
        true
    ));
    assert!(generic_label_matches(
        "*",
        0,
        "Carbocycloalkenyl",
        &aromatic_order_ring,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "Carbocycloalkenyl",
        &carbon_chain_with_double,
        true
    ));

    // The source per-ring witness must occur in every member ring of the
    // fused system. One unsaturated ring is insufficient; two pass.
    assert!(!generic_label_matches(
        "*",
        0,
        "Carbocycloalkenyl",
        &fused_saturated,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "Carbocycloalkenyl",
        &fused_one_unsaturated_ring,
        true
    ));
    assert!(generic_label_matches(
        "*",
        0,
        "Carbocycloalkenyl",
        &fused_two_unsaturated_rings,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "Carbocycloalkyl",
        &fused_one_unsaturated_ring,
        true
    ));

    // The H variant keeps its pinned degree-one IsHydrogen early return.
    assert!(generic_label_matches(
        "C*",
        1,
        "CarbocycloalkenylH",
        &carbon_hydrogen,
        true
    ));
}

#[test]
fn q55_fused_ring_match_keeps_source_system_and_component_boundaries() {
    let fused_carbon_rings = topology(
        &[Element::C; 10],
        &[
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 0),
            (3, 6),
            (6, 7),
            (7, 8),
            (8, 9),
            (9, 2),
        ],
    );
    let fused_hetero_ring = topology(
        &[
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::O,
            Element::C,
            Element::C,
            Element::C,
        ],
        &[
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 0),
            (3, 6),
            (6, 7),
            (7, 8),
            (8, 9),
            (9, 2),
        ],
    );
    let carbon_ring_with_disconnected_oxygen_ring = topology(
        &[
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::O,
            Element::O,
            Element::O,
        ],
        &[
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 0),
            (6, 7),
            (7, 8),
            (8, 6),
        ],
    );
    let carbon_chain = topology(&[Element::C, Element::C, Element::C], &[(0, 1), (1, 2)]);
    let carbon_hydrogen = topology(&[Element::C, Element::H], &[(0, 1)]);

    // FusedRingMatch adds a ring only when it shares at least two atoms with
    // the accumulated ring system. Heterocyclic's witness is over the whole
    // resulting system, so a heteroatom in the second fused ring suffices.
    assert!(!generic_label_matches(
        "*",
        0,
        "Heterocyclic",
        &fused_carbon_rings,
        true
    ));
    assert!(generic_label_matches(
        "*",
        0,
        "Heterocyclic",
        &fused_hetero_ring,
        true
    ));

    // A heteroatom in another connected component is not part of a carbon
    // root's fused ring system; an acyclic root also has no ring system.
    assert!(!generic_label_matches(
        "[C]",
        0,
        "Heterocyclic",
        &carbon_ring_with_disconnected_oxygen_ring,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "Heterocyclic",
        &carbon_chain,
        true
    ));

    // HeterocyclicH retains the source degree-one H shortcut before invoking
    // FusedRingMatch; the ordinary label still rejects an acyclic H root.
    assert!(!generic_label_matches(
        "C*",
        1,
        "Heterocyclic",
        &carbon_hydrogen,
        true
    ));
    assert!(generic_label_matches(
        "C*",
        1,
        "HeterocyclicH",
        &carbon_hydrogen,
        true
    ));
}

#[test]
fn q56_carbocycloalkyl_variants_keep_source_atom_bond_and_h_rules() {
    let ring_edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)];
    let carbon_ring = [Element::C; 6];
    let cyclohexane = topology(&carbon_ring, &ring_edges);
    let cyclohexene = topology_with_bond_orders(
        &carbon_ring,
        &ring_edges,
        &[
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
        ],
    );
    let flagged_ring = |aromatic_atom: bool, aromatic_bond: bool| {
        let atoms = (0..6)
            .map(|index| {
                Atom::from_spec(
                    AtomId::new(index),
                    AtomSpec::new(Element::C).with_aromatic(aromatic_atom && index == 0),
                )
            })
            .collect();
        let bonds = ring_edges
            .iter()
            .enumerate()
            .map(|(index, &(begin, end))| {
                Bond::from_spec(
                    BondId::new(index),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single)
                        .with_aromatic(aromatic_bond && index == 0),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("fixed aromatic-flagged carbocycloalkyl ring topology is valid")
    };
    let aromatic_atom_ring = flagged_ring(true, false);
    let aromatic_bond_ring = flagged_ring(false, true);
    let carbon_oxygen_ring = topology(
        &[
            Element::C,
            Element::C,
            Element::O,
            Element::C,
            Element::C,
            Element::C,
        ],
        &ring_edges,
    );
    let carbon_chain = topology(&[Element::C, Element::C, Element::C], &[(0, 1), (1, 2)]);
    let carbon_hydrogen = topology(&[Element::C, Element::H], &[(0, 1)]);
    let isolated_hydrogen = topology(&[Element::H], &[]);

    // The source requires every fused-system atom to be nonaromatic carbon
    // and every corresponding ring bond to be nonaromatic and single.
    assert!(generic_label_matches(
        "*",
        0,
        "Carbocycloalkyl",
        &cyclohexane,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "Carbocycloalkyl",
        &cyclohexene,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "Carbocycloalkyl",
        &aromatic_atom_ring,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "Carbocycloalkyl",
        &aromatic_bond_ring,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "Carbocycloalkyl",
        &carbon_oxygen_ring,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "Carbocycloalkyl",
        &carbon_chain,
        true
    ));
    // A second mapped ring atom remains in the by-value ignore mask and
    // causes source checkAtomRing rejection; the root itself is exempt.
    assert!(!generic_label_matches(
        "C*",
        1,
        "Carbocycloalkyl",
        &cyclohexane,
        true
    ));

    // CarbocycloalkylH keeps only the source degree-one hydrogen shortcut.
    assert!(!generic_label_matches(
        "C*",
        1,
        "Carbocycloalkyl",
        &carbon_hydrogen,
        true
    ));
    assert!(generic_label_matches(
        "C*",
        1,
        "CarbocycloalkylH",
        &carbon_hydrogen,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "CarbocycloalkylH",
        &isolated_hydrogen,
        true
    ));
}

#[test]
fn q57_carbocycloalkenyl_variants_keep_source_ring_atom_and_bond_witnesses() {
    let ring_edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)];
    let carbon_ring = [Element::C; 6];
    let cyclohexane = topology(&carbon_ring, &ring_edges);
    let cyclohexene = topology_with_bond_orders(
        &carbon_ring,
        &ring_edges,
        &[
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
        ],
    );
    let aromatic_order_ring = topology_with_bond_orders(
        &carbon_ring,
        &ring_edges,
        &[
            BondOrder::Aromatic,
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
        ],
    );
    let flagged_ring = |aromatic_atom: bool, aromatic_bond: bool| {
        let atoms = (0..6)
            .map(|index| {
                Atom::from_spec(
                    AtomId::new(index),
                    AtomSpec::new(Element::C).with_aromatic(aromatic_atom && index == 0),
                )
            })
            .collect();
        let bonds = ring_edges
            .iter()
            .enumerate()
            .map(|(index, &(begin, end))| {
                Bond::from_spec(
                    BondId::new(index),
                    BondSpec::new(
                        AtomId::new(begin),
                        AtomId::new(end),
                        if index == 0 && !aromatic_bond {
                            BondOrder::Double
                        } else {
                            BondOrder::Single
                        },
                    )
                    .with_aromatic(aromatic_bond && index == 0),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("fixed carbocycloalkenyl aromatic-state ring topology is valid")
    };
    let aromatic_atom_with_double = flagged_ring(true, false);
    let aromatic_bond_flag_ring = flagged_ring(false, true);
    let carbon_oxygen_ring = topology_with_bond_orders(
        &[
            Element::C,
            Element::C,
            Element::O,
            Element::C,
            Element::C,
            Element::C,
        ],
        &ring_edges,
        &[
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
        ],
    );
    let carbon_chain_with_double = topology_with_bond_orders(
        &[Element::C, Element::C, Element::C],
        &[(0, 1), (1, 2)],
        &[BondOrder::Double, BondOrder::Single],
    );
    let carbon_hydrogen = topology(&[Element::C, Element::H], &[(0, 1)]);
    let isolated_hydrogen = topology(&[Element::H], &[]);

    // Every fused-system atom must be carbon; each ring needs a bond whose
    // source state is aromatic or whose order is double/aromatic.
    assert!(!generic_label_matches(
        "*",
        0,
        "Carbocycloalkenyl",
        &cyclohexane,
        true
    ));
    assert!(generic_label_matches(
        "*",
        0,
        "Carbocycloalkenyl",
        &cyclohexene,
        true
    ));
    assert!(generic_label_matches(
        "*",
        0,
        "Carbocycloalkenyl",
        &aromatic_order_ring,
        true
    ));
    assert!(generic_label_matches(
        "*",
        0,
        "Carbocycloalkenyl",
        &aromatic_atom_with_double,
        true
    ));
    assert!(generic_label_matches(
        "*",
        0,
        "Carbocycloalkenyl",
        &aromatic_bond_flag_ring,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "Carbocycloalkenyl",
        &carbon_oxygen_ring,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "Carbocycloalkenyl",
        &carbon_chain_with_double,
        true
    ));
    assert!(!generic_label_matches(
        "C*",
        1,
        "Carbocycloalkenyl",
        &cyclohexene,
        true
    ));

    // The H label adds only the source degree-one shortcut; isolated H still
    // falls through and fails the fused-ring matcher.
    assert!(!generic_label_matches(
        "C*",
        1,
        "CarbocycloalkenylH",
        &carbon_chain_with_double,
        true
    ));
    assert!(generic_label_matches(
        "C*",
        1,
        "CarbocycloalkenylH",
        &carbon_hydrogen,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "CarbocycloalkenylH",
        &isolated_hydrogen,
        true
    ));
}

#[test]
fn q58_carboaryl_variants_keep_source_aromatic_atom_bond_and_h_rules() {
    let ring_edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)];
    let make_ring = |hetero_first_atom: bool,
                     nonaromatic_first_atom: bool,
                     bond_order: BondOrder,
                     aromatic_bond_flag: bool| {
        let atoms = (0..6)
            .map(|index| {
                let element = if hetero_first_atom && index == 0 {
                    Element::O
                } else {
                    Element::C
                };
                Atom::from_spec(
                    AtomId::new(index),
                    AtomSpec::new(element).with_aromatic(!(nonaromatic_first_atom && index == 0)),
                )
            })
            .collect();
        let bonds = ring_edges
            .iter()
            .enumerate()
            .map(|(index, &(begin, end))| {
                Bond::from_spec(
                    BondId::new(index),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), bond_order)
                        .with_aromatic(aromatic_bond_flag),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("fixed carboaryl ring topology is valid")
    };
    let aromatic_type_ring = make_ring(false, false, BondOrder::Aromatic, false);
    let aromatic_flag_ring = make_ring(false, false, BondOrder::Single, true);
    let nonaromatic_bond_ring = make_ring(false, false, BondOrder::Single, false);
    let nonaromatic_atom_ring = make_ring(false, true, BondOrder::Single, true);
    let heteroatom_ring = make_ring(true, false, BondOrder::Single, true);
    let carbon_hydrogen = topology(&[Element::C, Element::H], &[(0, 1)]);
    let isolated_hydrogen = topology(&[Element::H], &[]);

    // Carboaryl requires aromatic carbon at every ring atom and accepts the
    // source aromatic bond flag or aromatic bond type independently.
    assert!(generic_label_matches(
        "*",
        0,
        "Carboaryl",
        &aromatic_type_ring,
        true
    ));
    assert!(generic_label_matches(
        "*",
        0,
        "Carboaryl",
        &aromatic_flag_ring,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "Carboaryl",
        &nonaromatic_bond_ring,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "Carboaryl",
        &nonaromatic_atom_ring,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "Carboaryl",
        &heteroatom_ring,
        true
    ));

    // The generic matcher retains the mapped nonroot atom in its ignore mask;
    // checkAtomRing rejects that member even when aromatic predicates pass.
    assert!(!generic_label_matches(
        "c~*",
        1,
        "Carboaryl",
        &aromatic_type_ring,
        true
    ));

    // CarboarylH accepts only the source degree-one H shortcut; isolated H
    // still reaches the ring matcher and fails.
    assert!(!generic_label_matches(
        "C*",
        1,
        "Carboaryl",
        &carbon_hydrogen,
        true
    ));
    assert!(generic_label_matches(
        "C*",
        1,
        "CarboarylH",
        &carbon_hydrogen,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "CarboarylH",
        &isolated_hydrogen,
        true
    ));
}

#[test]
fn q59_carbocyclic_variants_keep_source_carbon_ring_and_h_rules() {
    let ring_edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)];
    let carbon_ring = [Element::C; 6];
    let cyclohexane = topology(&carbon_ring, &ring_edges);
    let cyclohexene = topology_with_bond_orders(
        &carbon_ring,
        &ring_edges,
        &[
            BondOrder::Double,
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
            BondOrder::Single,
        ],
    );
    let carbon_oxygen_ring = topology(
        &[
            Element::C,
            Element::C,
            Element::O,
            Element::C,
            Element::C,
            Element::C,
        ],
        &ring_edges,
    );
    let fused_carbon_oxygen_ring = topology(
        &[
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::O,
            Element::C,
            Element::C,
            Element::C,
        ],
        &[
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 0),
            (1, 6),
            (6, 7),
            (7, 8),
            (8, 9),
            (9, 0),
        ],
    );
    let separate_heterocycle = topology(
        &[
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::O,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
        ],
        &[
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 0),
            (6, 7),
            (7, 8),
            (8, 9),
            (9, 10),
            (10, 6),
        ],
    );
    let carbon_chain = topology(&[Element::C, Element::C, Element::C], &[(0, 1), (1, 2)]);
    let carbon_hydrogen = topology(&[Element::C, Element::H], &[(0, 1)]);
    let isolated_hydrogen = topology(&[Element::H], &[]);

    // Carbocyclic tests only atomic number across the root's fused ring set;
    // source imposes no aromaticity or bond-order filter.
    assert!(generic_label_matches(
        "*",
        0,
        "Carbocyclic",
        &cyclohexane,
        true
    ));
    assert!(generic_label_matches("*", 0, "CBC", &cyclohexane, true));
    assert!(generic_label_matches(
        "*",
        0,
        "Carbocyclic",
        &cyclohexene,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "Carbocyclic",
        &carbon_oxygen_ring,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "Carbocyclic",
        &fused_carbon_oxygen_ring,
        true
    ));
    assert!(generic_label_matches(
        "*",
        0,
        "Carbocyclic",
        &separate_heterocycle,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "Carbocyclic",
        &carbon_chain,
        true
    ));

    // As in the source by-value ignore path, mapping another ring atom makes
    // checkAtomRing reject the ring while the root itself remains exempt.
    assert!(!generic_label_matches(
        "C*",
        1,
        "Carbocyclic",
        &cyclohexane,
        true
    ));

    // CarbocyclicH short-circuits only for source degree-one H; isolated H
    // still reaches CarbocyclicAtomMatcher and fails its ring-membership test.
    assert!(!generic_label_matches(
        "C*",
        1,
        "Carbocyclic",
        &carbon_hydrogen,
        true
    ));
    assert!(generic_label_matches(
        "C*",
        1,
        "CarbocyclicH",
        &carbon_hydrogen,
        true
    ));
    assert!(generic_label_matches(
        "C*",
        1,
        "CBH",
        &carbon_hydrogen,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "CarbocyclicH",
        &isolated_hydrogen,
        true
    ));
}

#[test]
fn q60_no_carbon_ring_variants_keep_source_noncarbon_and_h_rules() {
    let ring_edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)];
    let noncarbon_ring = topology(
        &[
            Element::N,
            Element::O,
            Element::N,
            Element::O,
            Element::N,
            Element::O,
        ],
        &ring_edges,
    );
    let carbon_ring = topology(&[Element::C; 6], &ring_edges);
    let mixed_ring = topology(
        &[
            Element::N,
            Element::O,
            Element::C,
            Element::O,
            Element::N,
            Element::O,
        ],
        &ring_edges,
    );
    let fused_noncarbon_carbon_ring = topology(
        &[
            Element::N,
            Element::N,
            Element::N,
            Element::N,
            Element::N,
            Element::N,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
        ],
        &[
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 0),
            (1, 6),
            (6, 7),
            (7, 8),
            (8, 9),
            (9, 0),
        ],
    );
    let separate_carbon_ring = topology(
        &[
            Element::N,
            Element::O,
            Element::N,
            Element::O,
            Element::N,
            Element::O,
            Element::C,
            Element::C,
            Element::C,
        ],
        &[
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 0),
            (6, 7),
            (7, 8),
            (8, 6),
        ],
    );
    let acyclic_oxygen = topology(&[Element::C, Element::O], &[(0, 1)]);
    let carbon_hydrogen = topology(&[Element::C, Element::H], &[(0, 1)]);
    let isolated_hydrogen = topology(&[Element::H], &[]);

    // The source requires every atom in the root's fused ring system to have
    // atomic number other than carbon; it adds no bond or aromaticity filter.
    assert!(generic_label_matches(
        "*",
        0,
        "NoCarbonRing",
        &noncarbon_ring,
        true
    ));
    assert!(generic_label_matches("*", 0, "CXX", &noncarbon_ring, true));
    assert!(!generic_label_matches(
        "*",
        0,
        "NoCarbonRing",
        &carbon_ring,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "NoCarbonRing",
        &mixed_ring,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "NoCarbonRing",
        &fused_noncarbon_carbon_ring,
        true
    ));
    assert!(generic_label_matches(
        "*",
        0,
        "NoCarbonRing",
        &separate_carbon_ring,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "NoCarbonRing",
        &acyclic_oxygen,
        true
    ));

    // A mapped nonroot member is still rejected by checkAtomRing's ignore
    // rule, while a disconnected carbon ring stays outside this ring system.
    assert!(!generic_label_matches(
        "N*",
        1,
        "NoCarbonRing",
        &noncarbon_ring,
        true
    ));

    // NoCarbonRingH uses the same degree-one H shortcut; isolated H fails.
    assert!(!generic_label_matches(
        "C*",
        1,
        "NoCarbonRing",
        &carbon_hydrogen,
        true
    ));
    assert!(generic_label_matches(
        "C*",
        1,
        "NoCarbonRingH",
        &carbon_hydrogen,
        true
    ));
    assert!(generic_label_matches(
        "C*",
        1,
        "CXH",
        &carbon_hydrogen,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "NoCarbonRingH",
        &isolated_hydrogen,
        true
    ));
}

#[test]
fn q61_heterocyclic_variants_keep_source_ring_system_witness_and_h_rules() {
    let ring_edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)];
    let carbon_ring = topology(&[Element::C; 6], &ring_edges);
    let heterocyclic_ring = topology(
        &[
            Element::C,
            Element::C,
            Element::N,
            Element::C,
            Element::C,
            Element::C,
        ],
        &ring_edges,
    );
    let hydrogen_only_noncarbon_ring = topology(
        &[
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::H,
        ],
        &ring_edges,
    );
    let fused_hetero_ring = topology(
        &[
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::N,
            Element::C,
            Element::C,
            Element::C,
        ],
        &[
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 0),
            (1, 6),
            (6, 7),
            (7, 8),
            (8, 9),
            (9, 0),
        ],
    );
    let carbon_ring_with_pendant_n = topology(
        &[
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::N,
        ],
        &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (0, 6)],
    );
    let carbon_and_disconnected_hetero_ring = topology(
        &[
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::N,
            Element::O,
            Element::N,
            Element::O,
            Element::N,
            Element::O,
        ],
        &[
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 0),
            (6, 7),
            (7, 8),
            (8, 9),
            (9, 10),
            (10, 11),
            (11, 6),
        ],
    );
    let carbon_hydrogen = topology(&[Element::C, Element::H], &[(0, 1)]);
    let isolated_hydrogen = topology(&[Element::H], &[]);

    // Source FusedRingMatch imposes no per-atom predicate here. It accepts a
    // root carbon when a non-C/non-H witness exists in the same fused system.
    assert!(!generic_label_matches(
        "C",
        0,
        "Heterocyclic",
        &carbon_ring,
        true
    ));
    assert!(generic_label_matches(
        "C",
        0,
        "Heterocyclic",
        &heterocyclic_ring,
        true
    ));
    assert!(generic_label_matches(
        "C",
        0,
        "CHC",
        &heterocyclic_ring,
        true
    ));
    assert!(generic_label_matches(
        "C",
        0,
        "Heterocyclic",
        &fused_hetero_ring,
        true
    ));
    assert!(!generic_label_matches(
        "C",
        0,
        "Heterocyclic",
        &carbon_ring_with_pendant_n,
        true
    ));
    assert!(!generic_label_matches(
        "C",
        0,
        "Heterocyclic",
        &carbon_and_disconnected_hetero_ring,
        true
    ));
    assert!(!generic_label_matches(
        "C",
        0,
        "Heterocyclic",
        &hydrogen_only_noncarbon_ring,
        true
    ));

    // Existing ring members mapped by the query remain in the source ignore
    // mask even though the hetero witness is evaluated over the full set.
    assert!(!generic_label_matches(
        "C*",
        1,
        "Heterocyclic",
        &heterocyclic_ring,
        true
    ));

    // HeterocyclicH uses IsHydrogen before ring matching: degree-one H passes;
    // isolated H does not.
    assert!(!generic_label_matches(
        "C*",
        1,
        "Heterocyclic",
        &carbon_hydrogen,
        true
    ));
    assert!(generic_label_matches(
        "C*",
        1,
        "HeterocyclicH",
        &carbon_hydrogen,
        true
    ));
    assert!(generic_label_matches(
        "C*",
        1,
        "CHH",
        &carbon_hydrogen,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "HeterocyclicH",
        &isolated_hydrogen,
        true
    ));
}

#[test]
fn q62_heteroaryl_variants_keep_source_aromatic_ring_and_witness_rules() {
    let ring_edges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)];
    let make_ring = |hetero_first_atom: bool,
                     nonaromatic_first_atom: bool,
                     bond_order: BondOrder,
                     aromatic_bond_flag: bool| {
        let atoms = (0..6)
            .map(|index| {
                let element = if hetero_first_atom && index == 0 {
                    Element::N
                } else {
                    Element::C
                };
                Atom::from_spec(
                    AtomId::new(index),
                    AtomSpec::new(element).with_aromatic(!(nonaromatic_first_atom && index == 0)),
                )
            })
            .collect();
        let bonds = ring_edges
            .iter()
            .enumerate()
            .map(|(index, &(begin, end))| {
                Bond::from_spec(
                    BondId::new(index),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), bond_order)
                        .with_aromatic(aromatic_bond_flag),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("fixed heteroaryl ring topology is valid")
    };
    let aromatic_type_hetero_ring = make_ring(true, false, BondOrder::Aromatic, false);
    let aromatic_flag_hetero_ring = make_ring(true, false, BondOrder::Single, true);
    let nonaromatic_bond_hetero_ring = make_ring(true, false, BondOrder::Single, false);
    let nonaromatic_atom_hetero_ring = make_ring(true, true, BondOrder::Aromatic, false);
    let aromatic_carbon_ring = make_ring(false, false, BondOrder::Aromatic, false);
    let fused_heteroaryl = {
        let atoms = (0..10)
            .map(|index| {
                let element = if index == 6 { Element::N } else { Element::C };
                Atom::from_spec(
                    AtomId::new(index),
                    AtomSpec::new(element).with_aromatic(true),
                )
            })
            .collect();
        let edges = [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 0),
            (1, 6),
            (6, 7),
            (7, 8),
            (8, 9),
            (9, 0),
        ];
        let bonds = edges
            .iter()
            .enumerate()
            .map(|(index, &(begin, end))| {
                Bond::from_spec(
                    BondId::new(index),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Aromatic),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("fixed fused heteroaryl topology is valid")
    };
    let fused_nonaromatic_ring = {
        let atoms = (0..10)
            .map(|index| {
                let element = if index == 6 { Element::N } else { Element::C };
                Atom::from_spec(
                    AtomId::new(index),
                    AtomSpec::new(element).with_aromatic(index != 7),
                )
            })
            .collect();
        let edges = [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 0),
            (1, 6),
            (6, 7),
            (7, 8),
            (8, 9),
            (9, 0),
        ];
        let bonds = edges
            .iter()
            .enumerate()
            .map(|(index, &(begin, end))| {
                Bond::from_spec(
                    BondId::new(index),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Aromatic),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("fixed fused nonaromatic topology is valid")
    };
    let disconnected_heteroaryl = {
        let atoms = (0..12)
            .map(|index| {
                let element = if index < 6 {
                    Element::C
                } else if index % 2 == 0 {
                    Element::N
                } else {
                    Element::O
                };
                Atom::from_spec(
                    AtomId::new(index),
                    AtomSpec::new(element).with_aromatic(true),
                )
            })
            .collect();
        let edges = [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 0),
            (6, 7),
            (7, 8),
            (8, 9),
            (9, 10),
            (10, 11),
            (11, 6),
        ];
        let bonds = edges
            .iter()
            .enumerate()
            .map(|(index, &(begin, end))| {
                Bond::from_spec(
                    BondId::new(index),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Aromatic),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("fixed disconnected heteroaryl topology is valid")
    };
    let carbon_hydrogen = topology(&[Element::C, Element::H], &[(0, 1)]);
    let isolated_hydrogen = topology(&[Element::H], &[]);

    // Source requires aromatic atoms and bonds throughout each fused ring,
    // plus a non-C/non-H atom anywhere in that root's ring system.
    assert!(generic_label_matches(
        "c",
        0,
        "Heteroaryl",
        &aromatic_type_hetero_ring,
        true
    ));
    assert!(generic_label_matches(
        "c",
        0,
        "HAR",
        &aromatic_flag_hetero_ring,
        true
    ));
    assert!(!generic_label_matches(
        "c",
        0,
        "Heteroaryl",
        &nonaromatic_bond_hetero_ring,
        true
    ));
    assert!(!generic_label_matches(
        "c",
        0,
        "Heteroaryl",
        &nonaromatic_atom_hetero_ring,
        true
    ));
    assert!(!generic_label_matches(
        "c",
        0,
        "Heteroaryl",
        &aromatic_carbon_ring,
        true
    ));
    assert!(generic_label_matches(
        "c",
        0,
        "Heteroaryl",
        &fused_heteroaryl,
        true
    ));
    assert!(!generic_label_matches(
        "c",
        0,
        "Heteroaryl",
        &fused_nonaromatic_ring,
        true
    ));
    assert!(!generic_label_matches(
        "c",
        0,
        "Heteroaryl",
        &disconnected_heteroaryl,
        true
    ));

    // Mapped ring atoms remain ignored; the H variant keeps only the source
    // degree-one shortcut, with no shortcut for isolated hydrogen.
    assert!(!generic_label_matches(
        "c*",
        1,
        "Heteroaryl",
        &aromatic_type_hetero_ring,
        true
    ));
    assert!(!generic_label_matches(
        "C*",
        1,
        "Heteroaryl",
        &carbon_hydrogen,
        true
    ));
    assert!(generic_label_matches(
        "C*",
        1,
        "HeteroarylH",
        &carbon_hydrogen,
        true
    ));
    assert!(generic_label_matches(
        "C*",
        1,
        "HAH",
        &carbon_hydrogen,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "HeteroarylH",
        &isolated_hydrogen,
        true
    ));
}

#[test]
fn q63_cyclic_variants_keep_source_ring_membership_and_h_rules() {
    let triangle_edges = [(0, 1), (1, 2), (2, 0)];
    let cyclopropane = topology(&[Element::C; 3], &triangle_edges);
    let ring_with_acyclic_oxygen = topology(
        &[Element::C, Element::C, Element::C, Element::O],
        &[(0, 1), (1, 2), (2, 0), (0, 3)],
    );
    let carbon_chain = topology(&[Element::C, Element::C, Element::C], &[(0, 1), (1, 2)]);
    let carbon_hydrogen = topology(&[Element::C, Element::H], &[(0, 1)]);
    let isolated_hydrogen = topology(&[Element::H], &[]);

    assert!(generic_label_matches("*", 0, "Cyclic", &cyclopropane, true));
    assert!(generic_label_matches("*", 0, "CYC", &cyclopropane, true));
    assert!(!generic_label_matches(
        "C",
        0,
        "Cyclic",
        &carbon_chain,
        true
    ));
    assert!(!generic_label_matches(
        "O",
        0,
        "Cyclic",
        &ring_with_acyclic_oxygen,
        true
    ));
    assert!(generic_label_matches(
        "C",
        0,
        "Cyclic",
        &ring_with_acyclic_oxygen,
        true
    ));

    // A mapped ring member remains ignored by checkAtomRing, while the root
    // itself is exempt; only the H label has the degree-one shortcut.
    assert!(!generic_label_matches(
        "C*",
        1,
        "Cyclic",
        &cyclopropane,
        true
    ));
    assert!(!generic_label_matches(
        "C*",
        1,
        "Cyclic",
        &carbon_hydrogen,
        true
    ));
    assert!(generic_label_matches(
        "C*",
        1,
        "CyclicH",
        &carbon_hydrogen,
        true
    ));
    assert!(generic_label_matches(
        "C*",
        1,
        "CYH",
        &carbon_hydrogen,
        true
    ));
    assert!(!generic_label_matches(
        "*",
        0,
        "CyclicH",
        &isolated_hydrogen,
        true
    ));

    // The source reuses an already Fast RingInfo. Exercise that detached
    // SearchTarget path as well as the absent-RingInfo path above.
    let fast_ring_info =
        cosmolkit_core::fast_find_rings(&cyclopropane).expect("fixed cycle ring finding succeeds");
    let mut query = parse_smarts("C", &SmartsParseParams::default()).expect("parse query");
    query
        .atoms_mut()
        .get_mut(0)
        .expect("cyclic query atom exists")
        .set_prop("_QueryAtomGenericLabel", "Cyclic")
        .expect("set generic group label");
    let coordinates = CoordinateBlock::default();
    let target = SearchTarget::new(
        &cyclopropane,
        &coordinates,
        &cyclopropane.stereo_groups,
        Some(&fast_ring_info),
        None,
    );
    let params = SubstructMatchParams {
        use_generic_matchers: true,
        ..SubstructMatchParams::default()
    };
    assert!(!get_substruct_matches_with_params(&target, &query, &params).is_empty());
}
