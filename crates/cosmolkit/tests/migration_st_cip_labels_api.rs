#[path = "support/coordinate_views.rs"]
mod coordinate_views;

use std::error::Error as _;

use cosmolkit::{
    Atom, AtomId, AtomSpec, BINDING_CONTRACT, BindingExposure, BindingItem, BindingOwner,
    BindingParity, BindingSupport, BlockSet, Bond, BondId, BondOrder, BondSpec, BondStereo,
    ChiralTag, CipDescriptor, CipLabelOptions, CipLabelerError, Conformer2D, CoordinateBlock,
    Element, Molecule, MoleculeOpKind, MoleculeOpOutput, MoleculeProperties, OperationDomain,
    OperationError, ParityPolicy, StateModel, SupportStatus, TopologyBlock, TopologyEditKind,
    feature_spec, operation_invariant, operation_parity, operation_spec, support_matrix,
};

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

fn molecule(topology: TopologyBlock) -> Molecule {
    let atom_count = topology.atoms.len();
    Molecule::from_parts(
        topology,
        CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(
                17,
                (0..atom_count)
                    .map(|index| [index as f64, -(index as f64)])
                    .collect(),
            )],
            ..CoordinateBlock::default()
        },
        MoleculeProperties::default()
            .with_name("cip-label-public")
            .with_prop("source", "preserved")
            .unwrap(),
    )
    .unwrap()
}

fn tetrahedral_molecule(tag: ChiralTag) -> Molecule {
    molecule(topology(
        vec![
            AtomSpec::new(Element::C)
                .with_chiral_tag(tag)
                .with_prop("atom-user", "kept")
                .unwrap(),
            AtomSpec::new(Element::F),
            AtomSpec::new(Element::CL),
            AtomSpec::new(Element::BR),
            AtomSpec::new(Element::I),
        ],
        vec![
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single)
                .with_prop("bond-user", "kept")
                .unwrap(),
            BondSpec::new(AtomId::new(0), AtomId::new(2), BondOrder::Single),
            BondSpec::new(AtomId::new(0), AtomId::new(3), BondOrder::Single),
            BondSpec::new(AtomId::new(0), AtomId::new(4), BondOrder::Single),
        ],
    ))
}

fn sp2_molecule(stereo: BondStereo) -> Molecule {
    let ranked = |element, rank: &str| {
        AtomSpec::new(element)
            .with_computed_prop("_CIPRank", rank)
            .unwrap()
    };
    molecule(topology(
        vec![
            ranked(Element::F, "10"),
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            ranked(Element::CL, "20"),
            ranked(Element::BR, "30"),
            ranked(Element::I, "40"),
        ],
        vec![
            BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Double)
                .with_stereo(stereo)
                .with_stereo_atoms(AtomId::new(4), AtomId::new(5)),
            BondSpec::new(AtomId::new(1), AtomId::new(0), BondOrder::Single),
            BondSpec::new(AtomId::new(1), AtomId::new(4), BondOrder::Single),
            BondSpec::new(AtomId::new(2), AtomId::new(3), BondOrder::Single),
            BondSpec::new(AtomId::new(2), AtomId::new(5), BondOrder::Single),
        ],
    ))
}

fn atropisomer_molecule(stereo: BondStereo) -> Molecule {
    molecule(topology(
        vec![
            AtomSpec::new(Element::F),
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::CL),
            AtomSpec::new(Element::BR),
            AtomSpec::new(Element::I),
        ],
        vec![
            BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Single).with_stereo(stereo),
            BondSpec::new(AtomId::new(1), AtomId::new(4), BondOrder::Single),
            BondSpec::new(AtomId::new(1), AtomId::new(0), BondOrder::Single),
            BondSpec::new(AtomId::new(2), AtomId::new(5), BondOrder::Single),
            BondSpec::new(AtomId::new(2), AtomId::new(3), BondOrder::Single),
        ],
    ))
}

#[test]
fn canonical_public_signatures_and_defaults_are_exact() {
    let _: fn(&Molecule) -> Result<Molecule, OperationError> = Molecule::with_cip_labels;
    let _: for<'a, 'b> fn(&'a Molecule, &'b CipLabelOptions) -> Result<Molecule, OperationError> =
        Molecule::with_cip_labels_with_options;
    let _: fn(&mut Molecule) -> Result<(), OperationError> = Molecule::assign_cip_labels_;
    let _: for<'a, 'b> fn(&'a mut Molecule, &'b CipLabelOptions) -> Result<(), OperationError> =
        Molecule::assign_cip_labels_with_options_;

    let options = CipLabelOptions::default();
    assert_eq!(options.atoms(), None);
    assert_eq!(options.bonds(), None);
    assert_eq!(options.max_recursive_iterations(), 0);
}

#[test]
fn binding_contract_exposes_the_exact_eight_public_entries() {
    let expected = [
        "types.CipDescriptor",
        "types.CipDescriptorError",
        "types.CipLabelOptions",
        "types.CipLabelerError",
        "Molecule.with_cip_labels",
        "Molecule.with_cip_labels_with_options",
        "Molecule.assign_cip_labels_",
        "Molecule.assign_cip_labels_with_options_",
    ];
    let rows = BINDING_CONTRACT
        .iter()
        .filter(|row| expected.contains(&row.semantic_id))
        .collect::<Vec<_>>();
    assert_eq!(
        rows.iter().map(|row| row.semantic_id).collect::<Vec<_>>(),
        expected
    );
    assert!(rows.iter().all(|row| row.feature == "stereo"));
    assert!(
        rows.iter()
            .all(|row| row.exposure == BindingExposure::Public)
    );
    for row in &rows[..4] {
        assert_eq!(row.item, BindingItem::Type);
        assert_eq!(row.owner, BindingOwner::Type);
    }
    assert_eq!(rows[0].support, BindingSupport::SupportedWithRdkitParity);
    assert_eq!(rows[0].parity, BindingParity::RequiredNow);
    for row in &rows[1..4] {
        assert_eq!(row.support, BindingSupport::Supported);
        assert_eq!(row.parity, BindingParity::NotApplicable);
    }
    for row in &rows[4..] {
        assert_eq!(row.item, BindingItem::Callable);
        assert_eq!(row.owner, BindingOwner::Molecule);
        assert_eq!(row.support, BindingSupport::SupportedWithRdkitParity);
        assert_eq!(row.parity, BindingParity::RequiredNow);
    }
    assert_eq!(rows[4].javascript_name, "withCipLabels");
    assert_eq!(rows[5].javascript_name, "withCipLabelsWithOptions");
    assert_eq!(rows[6].javascript_name, "assignCipLabels");
    assert_eq!(rows[7].javascript_name, "assignCipLabelsWithOptions");
    assert_eq!(rows[4].callable.unwrap().parameters.len(), 0);
    assert_eq!(rows[5].callable.unwrap().parameters.len(), 1);
    assert_eq!(rows[6].callable.unwrap().parameters.len(), 0);
    assert_eq!(rows[7].callable.unwrap().parameters.len(), 1);
    assert_eq!(
        rows[4].callable.unwrap().state_model,
        StateModel::ValueReturning
    );
    assert_eq!(rows[6].callable.unwrap().state_model, StateModel::InPlace);
}

#[test]
fn generated_registry_and_all_four_matrices_share_one_single_output_operation() {
    let feature = feature_spec("stereo").unwrap();
    assert_eq!(feature.status, SupportStatus::SupportedWithRdkitParity);
    assert!(feature.rdkit_parity_sensitive);

    let spec = operation_spec("with_cip_labels_with_options").unwrap();
    assert_eq!(spec.output, MoleculeOpOutput::Single);
    assert_eq!(spec.result_type, "Molecule");
    assert_eq!(spec.domain, OperationDomain::Topology);
    assert_eq!(spec.kind, MoleculeOpKind::Weak);
    assert_eq!(spec.topology_edit, TopologyEditKind::Local);
    assert_eq!(spec.access.read(), BlockSet::NONE);
    assert_eq!(
        spec.access.write(),
        BlockSet::TOPOLOGY
            .union(BlockSet::PROPERTIES)
            .union(BlockSet::DERIVED_CACHE)
    );
    assert_eq!(spec.may_mutate, spec.access.write());
    assert_eq!(spec.auto_remap, BlockSet::NONE);
    assert_eq!(format!("{:?}", spec.requires_mapping), "None");
    assert_eq!(spec.derived_effects.recompute.bits(), 0);
    assert_eq!(
        spec.derived_effects.preserve.bits(),
        (1 << 0) | (1 << 1) | (1 << 2) | (1 << 3) | (1 << 5)
    );
    assert_eq!(
        spec.derived_effects.invalidate.bits(),
        (1 << 4) | (1 << 6) | (1 << 7)
    );
    assert_eq!(format!("{:?}", spec.cip_state), "Assign");
    assert_eq!(spec.support, SupportStatus::SupportedWithRdkitParity);
    assert_eq!(spec.parity, ParityPolicy::RequiredNow);
    assert!(!spec.io_roundtrip);
    assert_eq!(
        operation_invariant(spec.method).unwrap().profile,
        "weak_cip_label_assignment"
    );
    assert_eq!(
        operation_parity(spec.method).unwrap().profile,
        "assign_cip_labels_rdkit"
    );
    let support = support_matrix()
        .iter()
        .find(|row| {
            row.operation
                .is_some_and(|operation| std::ptr::eq(operation, spec))
        })
        .unwrap();
    assert!(std::ptr::eq(support.operation.unwrap(), spec));
}

#[test]
fn value_forms_assign_tetrahedral_state_preserve_props_and_keep_weak_mapping() {
    let source = tetrahedral_molecule(ChiralTag::TetrahedralCcw)
        .with_assigned_valence()
        .unwrap();
    let observer = source.clone();
    let short = source.with_cip_labels().unwrap();
    let explicit = source
        .with_cip_labels_with_options(&CipLabelOptions::default())
        .unwrap();

    assert_eq!(short, explicit);
    assert_eq!(source, observer);
    assert!(std::ptr::eq(source.topology(), observer.topology()));
    coordinate_views::assert_shared_coordinates(&source, &observer);
    assert!(std::ptr::eq(source.properties(), observer.properties()));
    assert!(!std::ptr::eq(source.topology(), short.topology()));
    assert!(!std::ptr::eq(source.properties(), short.properties()));
    coordinate_views::assert_shared_coordinates(&source, &short);

    assert_eq!(
        source.atom(AtomId::new(0)).unwrap().cip_descriptor(),
        Ok(None)
    );
    assert_eq!(
        short.atom(AtomId::new(0)).unwrap().cip_descriptor(),
        Ok(Some(CipDescriptor::S))
    );
    assert_eq!(
        short
            .atom(AtomId::new(0))
            .unwrap()
            .prop("_CIPNeighborOrder"),
        Some("[4,3,2,1]")
    );
    assert_eq!(
        short.atom(AtomId::new(0)).unwrap().prop("atom-user"),
        Some("kept")
    );
    assert_eq!(
        short.bond(BondId::new(0)).unwrap().prop("bond-user"),
        Some("kept")
    );
    assert_eq!(source.property("_CIPComputed"), None);
    assert_eq!(short.property("_CIPComputed"), Some("true"));
    assert!(short.properties().is_prop_computed("_CIPComputed"));
    assert_eq!(source.property("source"), Some("preserved"));
    assert_eq!(short.property("source"), Some("preserved"));
    assert!(format!("{source:?}").contains("derived_cache_is_empty: false"));
    assert!(format!("{short:?}").contains("derived_cache_is_empty: false"));

    assert_eq!(source.num_atoms(), short.num_atoms());
    assert_eq!(source.num_bonds(), short.num_bonds());
    assert!(
        short
            .atoms()
            .iter()
            .enumerate()
            .all(|(index, atom)| atom.id() == AtomId::new(index))
    );
    assert!(
        short
            .bonds()
            .iter()
            .enumerate()
            .all(|(index, bond)| bond.id() == BondId::new(index))
    );
}

#[test]
fn selection_and_repeated_assignment_preserve_unselected_state_exactly() {
    let mut topology = tetrahedral_molecule(ChiralTag::TetrahedralCw)
        .topology()
        .clone();
    topology.atoms[0].set_prop("_CIPCode", "old").unwrap();
    topology.atoms[0]
        .set_computed_prop("_CIPNeighborOrder", "[0]")
        .unwrap();
    let source = molecule(topology);

    let none_selected = source
        .with_cip_labels_with_options(
            &CipLabelOptions::default().with_bonds([BondId::new(0), BondId::new(0)]),
        )
        .unwrap();
    assert_eq!(
        none_selected.atom(AtomId::new(0)).unwrap().prop("_CIPCode"),
        Some("old")
    );
    assert_eq!(
        none_selected
            .atom(AtomId::new(0))
            .unwrap()
            .prop("_CIPNeighborOrder"),
        Some("[0]")
    );
    assert_eq!(none_selected.property("_CIPComputed"), Some("true"));

    let selected = none_selected
        .with_cip_labels_with_options(
            &CipLabelOptions::default().with_atoms([AtomId::new(0), AtomId::new(0)]),
        )
        .unwrap();
    assert_eq!(
        selected.atom(AtomId::new(0)).unwrap().cip_descriptor(),
        Ok(Some(CipDescriptor::R))
    );
    assert_eq!(
        selected
            .atom(AtomId::new(0))
            .unwrap()
            .prop("_CIPNeighborOrder"),
        Some("[4,3,2,1]")
    );
    assert_eq!(selected.with_cip_labels().unwrap(), selected);
}

#[test]
fn sp2_and_atropisomer_results_cross_the_public_transaction_without_reordering() {
    for (input, normalized, descriptor) in [
        (BondStereo::Z, BondStereo::Cis, CipDescriptor::Z),
        (BondStereo::E, BondStereo::Trans, CipDescriptor::E),
    ] {
        let output = sp2_molecule(input).with_cip_labels().unwrap();
        let axis = output.bond(BondId::new(0)).unwrap();
        assert_eq!(axis.cip_descriptor(), Ok(Some(descriptor)));
        assert_eq!(axis.stereo(), normalized);
        assert_eq!(axis.stereo_atoms(), Some([AtomId::new(4), AtomId::new(5)]));
        assert_eq!(axis.prop("_CIPNeighborOrder"), Some("[4,5]"));
    }

    for (stereo, descriptor) in [
        (BondStereo::AtropCcw, CipDescriptor::M),
        (BondStereo::AtropCw, CipDescriptor::P),
    ] {
        let output = atropisomer_molecule(stereo).with_cip_labels().unwrap();
        let axis = output.bond(BondId::new(0)).unwrap();
        assert_eq!(axis.cip_descriptor(), Ok(Some(descriptor)));
        assert_eq!(axis.stereo(), stereo);
        assert_eq!(axis.stereo_atoms(), None);
        assert_eq!(axis.prop("_CIPNeighborOrder"), Some("[4,5]"));
    }
}

#[test]
fn inplace_forms_match_value_results_and_commit_only_after_success() {
    let source = tetrahedral_molecule(ChiralTag::TetrahedralCcw);
    let expected = source.with_cip_labels().unwrap();

    let mut short = source.clone();
    short.assign_cip_labels_().unwrap();
    assert_eq!(short, expected);

    let mut explicit = source;
    explicit
        .assign_cip_labels_with_options_(&CipLabelOptions::default())
        .unwrap();
    assert_eq!(explicit, expected);
}

#[test]
fn typed_failures_are_atomic_for_value_and_inplace_entrypoints() {
    let source = tetrahedral_molecule(ChiralTag::TetrahedralCcw);
    let observer = source.clone();
    let invalid = CipLabelOptions::default().with_atoms([AtomId::new(5)]);
    let error = source.with_cip_labels_with_options(&invalid).unwrap_err();
    assert_eq!(
        error,
        OperationError::CipLabeler(CipLabelerError::AtomIndexOutOfRange {
            index: 5,
            atom_count: 5,
        })
    );
    assert!(error.source().is_some());
    assert_eq!(source, observer);
    assert!(std::ptr::eq(source.topology(), observer.topology()));
    coordinate_views::assert_shared_coordinates(&source, &observer);
    assert!(std::ptr::eq(source.properties(), observer.properties()));

    let mut target = source.clone();
    let target_observer = target.clone();
    assert_eq!(
        target.assign_cip_labels_with_options_(&invalid),
        Err(OperationError::CipLabeler(
            CipLabelerError::AtomIndexOutOfRange {
                index: 5,
                atom_count: 5,
            }
        ))
    );
    assert_eq!(target, source);
    assert!(std::ptr::eq(target.topology(), target_observer.topology()));
    coordinate_views::assert_shared_coordinates(&target, &target_observer);
    assert!(std::ptr::eq(
        target.properties(),
        target_observer.properties()
    ));

    let unsupported = molecule(topology(
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
    ));
    assert_eq!(
        unsupported.with_cip_labels(),
        Err(OperationError::CipLabeler(
            CipLabelerError::UnsupportedConfiguration {
                atom: 0,
                tag: ChiralTag::SquarePlanar,
            }
        ))
    );
}

#[test]
fn single_output_contract_makes_multi_candidate_semantics_not_applicable() {
    let spec = operation_spec("with_cip_labels_with_options").unwrap();
    assert_eq!(spec.output, MoleculeOpOutput::Single);
    assert_eq!(spec.result_type, "Molecule");
    assert!(
        BINDING_CONTRACT
            .iter()
            .filter(|row| row.semantic_id.starts_with("Molecule.with_cip_labels"))
            .all(|row| row.callable.unwrap().state_model == StateModel::ValueReturning)
    );
}
