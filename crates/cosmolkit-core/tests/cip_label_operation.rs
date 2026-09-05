use cosmolkit_core::{
    AtomId, AtomSpec, BlockSet, BondId, BondOrder, BondSpec, ChiralTag, CipDescriptor, CipLabelOptions,
    CipLabelerError, DerivedState, Element, MOLECULE_OPS, MappingRequirement, Molecule, MoleculeBuilder,
    MoleculeOpKind, OPERATION_INVARIANT_MATRIX, OperationDomain, OperationError, PARITY_MATRIX, ParityPolicy,
    SUPPORT_MATRIX, SupportStatus, TopologyEditKind,
};

fn add_tetrahedral_center(builder: &mut MoleculeBuilder, chiral_tag: ChiralTag) -> AtomId {
    let center = builder.add_atom(AtomSpec::new(Element::C).with_chiral_tag(chiral_tag));
    for element in [Element::F, Element::CL, Element::BR, Element::I] {
        let neighbor = builder.add_atom(AtomSpec::new(element));
        builder
            .add_bond(BondSpec::new(center, neighbor, BondOrder::Single))
            .expect("tetrahedral bond must build");
    }
    center
}

fn add_stereo_double_bond(builder: &mut MoleculeBuilder) -> BondId {
    let begin = builder.add_atom(AtomSpec::new(Element::C));
    let end = builder.add_atom(AtomSpec::new(Element::C));
    let begin_low = builder.add_atom(AtomSpec::new(Element::F));
    let begin_high = builder.add_atom(AtomSpec::new(Element::BR));
    let end_low = builder.add_atom(AtomSpec::new(Element::CL));
    let end_high = builder.add_atom(AtomSpec::new(Element::I));
    let double_bond = builder
        .add_bond(
            BondSpec::new(begin, end, BondOrder::Double)
                .with_stereo(cosmolkit_core::BondStereo::Cis)
                .with_stereo_atoms(begin_high, end_high),
        )
        .expect("stereo double bond must build");
    for (focus, neighbor) in [(begin, begin_low), (begin, begin_high), (end, end_low), (end, end_high)] {
        builder
            .add_bond(BondSpec::new(focus, neighbor, BondOrder::Single))
            .expect("double-bond substituent must build");
    }
    double_bond
}

fn mixed_stereo_molecule() -> (Molecule, AtomId, AtomId, BondId) {
    let mut builder = MoleculeBuilder::new();
    let first = add_tetrahedral_center(&mut builder, ChiralTag::TetrahedralCcw);
    let second = add_tetrahedral_center(&mut builder, ChiralTag::TetrahedralCw);
    let double_bond = add_stereo_double_bond(&mut builder);
    builder
        .add_3d_conformer(vec![[0.0, 0.0, 0.0]; 16])
        .expect("coordinates must align with all atoms");
    (
        builder.build().expect("mixed stereo molecule must build"),
        first,
        second,
        double_bond,
    )
}

#[test]
fn cip_label_registry_contract_declares_exact_capabilities_and_matrices() {
    let operation = MOLECULE_OPS
        .iter()
        .copied()
        .find(|operation| operation.method == "with_cip_labels_with_options")
        .expect("modern CIP assignment must be registered");

    assert_eq!(operation.impl_fn, "assigned_cip_labels_impl");
    assert_eq!(operation.domain, OperationDomain::Topology);
    assert_eq!(operation.kind, MoleculeOpKind::Weak);
    assert_eq!(operation.topology_edit, TopologyEditKind::Local);
    assert_eq!(operation.access.read(), BlockSet::NONE);
    assert_eq!(
        operation.access.write(),
        BlockSet::TOPOLOGY
            .union(BlockSet::PROPERTIES)
            .union(BlockSet::DERIVED_CACHE)
    );
    assert_eq!(operation.may_mutate, operation.access.write());
    assert_eq!(operation.auto_remap, BlockSet::NONE);
    assert_eq!(operation.requires_mapping, MappingRequirement::None);
    assert_eq!(operation.parity, ParityPolicy::RequiredWhenSupported);
    assert!(operation.io_roundtrip);
    assert_eq!(operation.derived_effects.recompute(), DerivedState::NONE);
    assert_eq!(operation.derived_effects.preserve(), DerivedState::NONE);
    assert_eq!(
        operation.derived_effects.invalidate(),
        DerivedState::STEREO
            .union(DerivedState::DRAWING)
            .union(DerivedState::FINGERPRINT)
    );
    assert_eq!(
        operation.support,
        SupportStatus::SupportedWithRdkitParity {
            rdkit_version: "2026.03.1",
        }
    );

    let support = SUPPORT_MATRIX
        .iter()
        .find(|entry| entry.operation == Some(operation))
        .expect("modern CIP operation must have a support entry");
    assert_eq!(support.feature.name, "stereo.cip_labeler");
    assert_eq!(
        OPERATION_INVARIANT_MATRIX
            .iter()
            .find(|entry| entry.operation == operation)
            .expect("modern CIP operation must have an invariant entry")
            .profile,
        "weak_topology_state"
    );
    assert_eq!(
        PARITY_MATRIX
            .iter()
            .find(|entry| entry.operation == operation)
            .expect("modern CIP operation must have a parity entry")
            .profile,
        "modern_ciplabeler_rdkit"
    );
}

#[test]
fn cip_label_value_and_in_place_forms_preserve_source_indices_and_coordinates() {
    let (molecule, first, second, double_bond) = mixed_stereo_molecule();
    let source = molecule.clone();
    let source_coordinates = molecule.conformers_3d().to_vec();
    let source_atom_ids = molecule.atoms().iter().map(|atom| atom.id()).collect::<Vec<_>>();
    let source_bond_ids = molecule.bonds().iter().map(|bond| bond.id()).collect::<Vec<_>>();

    let labeled = molecule.with_cip_labels().expect("full assignment must succeed");
    assert_eq!(molecule, source, "value-style assignment mutated its source");
    assert_eq!(labeled.conformers_3d(), source_coordinates);
    assert_eq!(
        labeled.atoms().iter().map(|atom| atom.id()).collect::<Vec<_>>(),
        source_atom_ids
    );
    assert_eq!(
        labeled.bonds().iter().map(|bond| bond.id()).collect::<Vec<_>>(),
        source_bond_ids
    );
    assert!(labeled.atoms()[first.index()].cip_descriptor().unwrap().is_some());
    assert!(labeled.atoms()[second.index()].cip_descriptor().unwrap().is_some());
    assert!(matches!(
        labeled.bonds()[double_bond.index()].cip_descriptor().unwrap(),
        Some(CipDescriptor::E | CipDescriptor::Z)
    ));
    assert_eq!(labeled.prop("_CIPComputed"), Some("1"));

    let shared_source = molecule.clone();
    let mut in_place = molecule;
    in_place.assign_cip_labels_().expect("in-place assignment must succeed");
    assert_eq!(shared_source, source, "COW assignment escaped to a shared clone");
    assert_eq!(in_place, labeled);
}

#[test]
fn cip_label_selection_matches_pinned_wrapper_truth_value_dispatch() {
    let (molecule, first, second, double_bond) = mixed_stereo_molecule();

    let atom_only = molecule
        .with_cip_labels_with_options(CipLabelOptions::default().with_atoms([first]))
        .expect("selected atom assignment must succeed");
    assert!(atom_only.atoms()[first.index()].cip_descriptor().unwrap().is_some());
    assert_eq!(atom_only.atoms()[second.index()].cip_descriptor().unwrap(), None);
    assert_eq!(atom_only.bonds()[double_bond.index()].cip_descriptor().unwrap(), None);

    let bond_only = molecule
        .with_cip_labels_with_options(CipLabelOptions::default().with_bonds([double_bond]))
        .expect("selected bond assignment must succeed");
    assert_eq!(bond_only.atoms()[first.index()].cip_descriptor().unwrap(), None);
    assert!(matches!(
        bond_only.bonds()[double_bond.index()].cip_descriptor().unwrap(),
        Some(CipDescriptor::E | CipDescriptor::Z)
    ));

    let simultaneous = molecule
        .with_cip_labels_with_options(
            CipLabelOptions::default()
                .with_atoms([second])
                .with_bonds([double_bond]),
        )
        .expect("simultaneous selected assignment must succeed");
    assert_eq!(simultaneous.atoms()[first.index()].cip_descriptor().unwrap(), None);
    assert!(simultaneous.atoms()[second.index()].cip_descriptor().unwrap().is_some());
    assert!(
        simultaneous.bonds()[double_bond.index()]
            .cip_descriptor()
            .unwrap()
            .is_some()
    );

    let empty = molecule
        .with_cip_labels_with_options(CipLabelOptions::default().with_atoms([]))
        .expect("two false-valued selections must dispatch to full assignment");
    assert!(
        empty
            .atoms()
            .iter()
            .filter(|atom| atom.chiral_tag() != ChiralTag::Unspecified)
            .all(|atom| atom.cip_descriptor().unwrap().is_some())
    );
    assert!(matches!(
        empty.bonds()[double_bond.index()].cip_descriptor().unwrap(),
        Some(CipDescriptor::E | CipDescriptor::Z)
    ));

    let empty_atoms_with_bond = molecule
        .with_cip_labels_with_options(CipLabelOptions::default().with_atoms([]).with_bonds([double_bond]))
        .expect("a non-empty bond selection must use the selected overload");
    assert!(
        empty_atoms_with_bond
            .atoms()
            .iter()
            .all(|atom| atom.cip_descriptor().unwrap().is_none())
    );
    assert!(
        empty_atoms_with_bond.bonds()[double_bond.index()]
            .cip_descriptor()
            .unwrap()
            .is_some()
    );
}

#[test]
fn cip_label_recursion_error_is_structured_and_in_place_storage_remains_valid() {
    let molecule = Molecule::from_smiles("C/C=C/[C@@H](/C=C/O)[C@H](C)[C@H](/C=C/C)/C=C/O")
        .expect("recursive CIP fixture must parse");
    let options = CipLabelOptions::default().with_max_recursive_iterations(1);
    let error = molecule
        .with_cip_labels_with_options(options.clone())
        .expect_err("the source recursion budget must be enforced");
    assert!(matches!(
        error,
        OperationError::CipLabeler {
            source: CipLabelerError::MaxIterationsExceeded,
            ..
        }
    ));

    let atom_count = molecule.num_atoms();
    let bond_count = molecule.num_bonds();
    let mut in_place = molecule;
    let error = in_place
        .assign_cip_labels_with_options_(options)
        .expect_err("in-place recursion limit must be enforced");
    assert!(matches!(
        error,
        OperationError::CipLabeler {
            source: CipLabelerError::MaxIterationsExceeded,
            ..
        }
    ));
    assert_eq!(in_place.num_atoms(), atom_count);
    assert_eq!(in_place.num_bonds(), bond_count);
    for bond in in_place.bonds() {
        assert!(bond.begin().index() < atom_count);
        assert!(bond.end().index() < atom_count);
    }
    in_place
        .with_assigned_valence()
        .expect("a failed in-place CIP call must leave usable complete storage");
}
