use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, BondStereo, ChiralTag, CipDescriptor, Element,
    MoleculeProperties, TopologyBlock,
};
use cosmolkit_stereo::{CipLabelOptions, CipLabelerError, assign_cip_labels};
use cosmolkit_types::BondOrder;

fn topology(atom_specs: Vec<AtomSpec>, bond_specs: Vec<BondSpec>) -> TopologyBlock {
    let atoms = atom_specs
        .into_iter()
        .enumerate()
        .map(|(index, spec)| Atom::from_spec(AtomId::new(index), spec))
        .collect();
    let bonds = bond_specs
        .into_iter()
        .enumerate()
        .map(|(index, spec)| Bond::from_spec(BondId::new(index), spec))
        .collect();
    TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap()
}

fn tetrahedral_topology(tag: ChiralTag) -> TopologyBlock {
    topology(
        vec![
            AtomSpec::new(Element::C).with_chiral_tag(tag),
            AtomSpec::new(Element::F),
            AtomSpec::new(Element::CL),
            AtomSpec::new(Element::BR),
            AtomSpec::new(Element::I),
        ],
        vec![
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            BondSpec::new(AtomId::new(0), AtomId::new(2), BondOrder::Single),
            BondSpec::new(AtomId::new(0), AtomId::new(3), BondOrder::Single),
            BondSpec::new(AtomId::new(0), AtomId::new(4), BondOrder::Single),
        ],
    )
}

fn sp2_topology(stereo: BondStereo, explicit_stereo_atoms: bool) -> TopologyBlock {
    let ranked = |element, rank: &str| {
        AtomSpec::new(element)
            .with_computed_prop("_CIPRank", rank)
            .unwrap()
    };
    let atoms = vec![
        ranked(Element::F, "10"),
        AtomSpec::new(Element::C),
        AtomSpec::new(Element::C),
        ranked(Element::CL, "20"),
        ranked(Element::BR, "30"),
        ranked(Element::I, "40"),
    ];
    let mut axis =
        BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Double).with_stereo(stereo);
    if explicit_stereo_atoms {
        axis = axis.with_stereo_atoms(AtomId::new(4), AtomId::new(5));
    }
    topology(
        atoms,
        vec![
            axis,
            BondSpec::new(AtomId::new(1), AtomId::new(0), BondOrder::Single),
            BondSpec::new(AtomId::new(1), AtomId::new(4), BondOrder::Single),
            BondSpec::new(AtomId::new(2), AtomId::new(3), BondOrder::Single),
            BondSpec::new(AtomId::new(2), AtomId::new(5), BondOrder::Single),
        ],
    )
}

fn atropisomer_topology(stereo: BondStereo, low_index_carriers_first: bool) -> TopologyBlock {
    let mut bonds =
        vec![BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Single).with_stereo(stereo)];
    let begin_carriers = if low_index_carriers_first {
        [(1, 0), (1, 4)]
    } else {
        [(1, 4), (1, 0)]
    };
    let end_carriers = if low_index_carriers_first {
        [(2, 3), (2, 5)]
    } else {
        [(2, 5), (2, 3)]
    };
    for (begin, end) in begin_carriers.into_iter().chain(end_carriers) {
        bonds.push(BondSpec::new(
            AtomId::new(begin),
            AtomId::new(end),
            BondOrder::Single,
        ));
    }
    topology(
        vec![
            AtomSpec::new(Element::F),
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::CL),
            AtomSpec::new(Element::BR),
            AtomSpec::new(Element::I),
        ],
        bonds,
    )
}

#[test]
fn cip_label_options_preserve_source_selection_and_counter_shape() {
    let default = CipLabelOptions::default();
    assert_eq!(default.atoms(), None);
    assert_eq!(default.bonds(), None);
    assert_eq!(default.max_recursive_iterations(), 0);

    let options = CipLabelOptions::default()
        .with_atoms([AtomId::new(2), AtomId::new(2)])
        .with_bonds([BondId::new(3)])
        .with_max_recursive_iterations(u32::MAX);
    assert_eq!(options.atoms(), Some(&[AtomId::new(2), AtomId::new(2)][..]));
    assert_eq!(options.bonds(), Some(&[BondId::new(3)][..]));
    assert_eq!(options.max_recursive_iterations(), u32::MAX);
}

#[test]
fn no_configuration_still_sets_computed_completion_and_preserves_user_props() {
    let properties = MoleculeProperties::default()
        .with_prop("user", "kept")
        .unwrap()
        .with_prop("_CIPComputed", "old")
        .unwrap();
    let assignment = assign_cip_labels(
        TopologyBlock::default(),
        properties,
        &CipLabelOptions::default(),
    )
    .unwrap();

    assert!(assignment.topology().atoms.is_empty());
    assert_eq!(assignment.properties().prop("user"), Some("kept"));
    assert_eq!(assignment.properties().prop("_CIPComputed"), Some("true"));
    assert!(assignment.properties().is_prop_computed("_CIPComputed"));
}

#[test]
fn full_tetrahedral_assignment_writes_exact_primary_and_neighbor_order_state() {
    for (tag, expected) in [
        (ChiralTag::TetrahedralCcw, CipDescriptor::S),
        (ChiralTag::TetrahedralCw, CipDescriptor::R),
    ] {
        let assignment = assign_cip_labels(
            tetrahedral_topology(tag),
            MoleculeProperties::default(),
            &CipLabelOptions::default(),
        )
        .unwrap();
        let center = &assignment.topology().atoms[0];
        assert_eq!(center.cip_descriptor().unwrap(), Some(expected));
        assert_eq!(center.prop("_CIPNeighborOrder"), Some("[4,3,2,1]"));
        assert!(!center.is_prop_computed("_CIPCode"));
        assert!(center.is_prop_computed("_CIPNeighborOrder"));
    }
}

#[test]
fn selection_masks_are_exact_duplicates_are_idempotent_and_repeated_calls_preserve_state() {
    let mut input = tetrahedral_topology(ChiralTag::TetrahedralCcw);
    input.atoms[0].set_prop("_CIPCode", "old").unwrap();
    input.atoms[0]
        .set_computed_prop("_CIPNeighborOrder", "[0]")
        .unwrap();

    let none_selected = assign_cip_labels(
        input,
        MoleculeProperties::default(),
        &CipLabelOptions::default().with_bonds([BondId::new(0), BondId::new(0)]),
    )
    .unwrap();
    assert_eq!(
        none_selected.topology().atoms[0].prop("_CIPCode"),
        Some("old")
    );
    assert_eq!(
        none_selected.topology().atoms[0].prop("_CIPNeighborOrder"),
        Some("[0]")
    );

    let (topology, properties) = none_selected.into_parts();
    let selected = assign_cip_labels(
        topology,
        properties,
        &CipLabelOptions::default().with_atoms([AtomId::new(0), AtomId::new(0)]),
    )
    .unwrap();
    assert_eq!(selected.topology().atoms[0].prop("_CIPCode"), Some("S"));
    assert_eq!(
        selected.topology().atoms[0].prop("_CIPNeighborOrder"),
        Some("[4,3,2,1]")
    );

    let (topology, properties) = selected.into_parts();
    let all_again = assign_cip_labels(
        topology,
        properties,
        &CipLabelOptions::default().with_atoms([]).with_bonds([]),
    )
    .unwrap();
    assert_eq!(all_again.topology().atoms[0].prop("_CIPCode"), Some("S"));
    assert!(all_again.properties().is_prop_computed("_CIPComputed"));
}

#[test]
fn invalid_typed_selection_reports_the_first_exact_category_error() {
    let error = assign_cip_labels(
        tetrahedral_topology(ChiralTag::TetrahedralCcw),
        MoleculeProperties::default(),
        &CipLabelOptions::default().with_atoms([AtomId::new(5), AtomId::new(6)]),
    )
    .unwrap_err();
    assert_eq!(
        error,
        CipLabelerError::AtomIndexOutOfRange {
            index: 5,
            atom_count: 5,
        }
    );
    assert_eq!(
        error.to_string(),
        "CIPLabeler atom index 5 is out of range for 5 atoms"
    );

    let error = assign_cip_labels(
        tetrahedral_topology(ChiralTag::TetrahedralCcw),
        MoleculeProperties::default(),
        &CipLabelOptions::default().with_bonds([BondId::new(4)]),
    )
    .unwrap_err();
    assert_eq!(
        error,
        CipLabelerError::BondIndexOutOfRange {
            index: 4,
            bond_count: 4,
        }
    );
}

#[test]
fn sp2_explicit_carriers_assign_uppercase_labels_and_normalize_source_stereo() {
    for (input_stereo, expected_stereo, expected) in [
        (BondStereo::Z, BondStereo::Cis, CipDescriptor::Z),
        (BondStereo::E, BondStereo::Trans, CipDescriptor::E),
    ] {
        let assignment = assign_cip_labels(
            sp2_topology(input_stereo, true),
            MoleculeProperties::default(),
            &CipLabelOptions::default(),
        )
        .unwrap();
        let axis = &assignment.topology().bonds[0];
        assert_eq!(axis.cip_descriptor().unwrap(), Some(expected));
        assert_eq!(axis.prop("_CIPNeighborOrder"), Some("[4,5]"));
        assert_eq!(axis.stereo(), expected_stereo);
        assert_eq!(axis.stereo_atoms(), Some([AtomId::new(4), AtomId::new(5)]));
        assert!(axis.is_prop_computed("_CIPNeighborOrder"));
        assert!(!axis.is_prop_computed("_CIPCode"));
    }
}

#[test]
fn sp2_rank_fallback_is_exact_and_missing_or_tied_rank_fails_closed() {
    let assignment = assign_cip_labels(
        sp2_topology(BondStereo::E, false),
        MoleculeProperties::default(),
        &CipLabelOptions::default(),
    )
    .unwrap();
    assert_eq!(assignment.topology().bonds[0].prop("_CIPCode"), Some("E"));
    assert_eq!(
        assignment.topology().bonds[0].stereo_atoms(),
        Some([AtomId::new(4), AtomId::new(5)])
    );

    let mut missing = sp2_topology(BondStereo::E, false);
    missing.atoms[4].clear_prop("_CIPRank");
    assert!(matches!(
        assign_cip_labels(
            missing,
            MoleculeProperties::default(),
            &CipLabelOptions::default()
        ),
        Err(CipLabelerError::IncorrectNumberOfStereoAtoms)
    ));

    let mut tied = sp2_topology(BondStereo::E, false);
    tied.atoms[0].set_prop("_CIPRank", "30").unwrap();
    assert!(matches!(
        assign_cip_labels(
            tied,
            MoleculeProperties::default(),
            &CipLabelOptions::default()
        ),
        Err(CipLabelerError::IncorrectNumberOfStereoAtoms)
    ));
}

#[test]
fn atropisomer_assignment_uses_core_carrier_order_and_preserves_axis_stereo() {
    for (stereo, expected) in [
        (BondStereo::AtropCcw, CipDescriptor::M),
        (BondStereo::AtropCw, CipDescriptor::P),
    ] {
        let assignment = assign_cip_labels(
            atropisomer_topology(stereo, false),
            MoleculeProperties::default(),
            &CipLabelOptions::default(),
        )
        .unwrap();
        let axis = &assignment.topology().bonds[0];
        assert_eq!(axis.cip_descriptor().unwrap(), Some(expected));
        assert_eq!(axis.prop("_CIPNeighborOrder"), Some("[4,5]"));
        assert_eq!(axis.stereo(), stereo);
        assert_eq!(axis.stereo_atoms(), None);
    }
}

#[test]
fn unsupported_non_tetrahedral_configuration_is_structured_and_does_not_guess() {
    let input = topology(
        vec![
            AtomSpec::new(Element::PT).with_chiral_tag(ChiralTag::SquarePlanar),
            AtomSpec::new(Element::F),
            AtomSpec::new(Element::CL),
            AtomSpec::new(Element::BR),
            AtomSpec::new(Element::I),
        ],
        vec![
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            BondSpec::new(AtomId::new(0), AtomId::new(2), BondOrder::Single),
            BondSpec::new(AtomId::new(0), AtomId::new(3), BondOrder::Single),
            BondSpec::new(AtomId::new(0), AtomId::new(4), BondOrder::Single),
        ],
    );
    let error = assign_cip_labels(
        input,
        MoleculeProperties::default(),
        &CipLabelOptions::default(),
    )
    .unwrap_err();
    assert_eq!(
        error,
        CipLabelerError::UnsupportedConfiguration {
            atom: 0,
            tag: ChiralTag::SquarePlanar,
        }
    );
}

#[test]
fn source_visible_limit_and_node_cap_diagnostics_remain_exact() {
    assert_eq!(
        CipLabelerError::MaxIterationsExceeded.to_string(),
        "Max Iterations Exceeded in CIP label calculation"
    );
    assert_eq!(
        CipLabelerError::TooManyNodes { limit: 100_000 }.to_string(),
        "Digraph generation failed: more than 100000 nodes found."
    );

    let assignment = assign_cip_labels(
        tetrahedral_topology(ChiralTag::TetrahedralCcw),
        MoleculeProperties::default(),
        &CipLabelOptions::default().with_max_recursive_iterations(1),
    )
    .unwrap();
    assert_eq!(assignment.topology().atoms[0].prop("_CIPCode"), Some("S"));
}

#[test]
fn all_twelve_persisted_descriptor_spellings_remain_typed_at_the_assignment_boundary() {
    for (spelling, expected) in [
        ("R", CipDescriptor::R),
        ("S", CipDescriptor::S),
        ("r", CipDescriptor::LowerR),
        ("s", CipDescriptor::LowerS),
        ("E", CipDescriptor::E),
        ("Z", CipDescriptor::Z),
        ("e", CipDescriptor::LowerE),
        ("z", CipDescriptor::LowerZ),
        ("M", CipDescriptor::M),
        ("P", CipDescriptor::P),
        ("m", CipDescriptor::LowerM),
        ("p", CipDescriptor::LowerP),
    ] {
        assert_eq!(spelling.parse::<CipDescriptor>().unwrap(), expected);
        assert_eq!(expected.as_str(), spelling);
    }
}
