use cosmolkit_core::{
    AtomSpec, BlockSet, BondOrder, BondSpec, ChiralTag, DerivedState, Element, MOLECULE_OPS, MappingRequirement,
    Molecule, MoleculeBuilder, MoleculeOpKind, OPERATION_INVARIANT_MATRIX, OperationDomain, OperationError,
    PARITY_MATRIX, ParityPolicy, SUPPORT_MATRIX, StereoError, SupportStatus, TopologyEditKind,
};

const CENTER: usize = 0;

fn tetrahedral_star(center: AtomSpec) -> Molecule {
    let mut builder = MoleculeBuilder::new();
    let center = builder.add_atom(center);
    for element in [Element::F, Element::CL, Element::BR] {
        let neighbor = builder.add_atom(AtomSpec::new(element));
        builder
            .add_bond(BondSpec::new(center, neighbor, BondOrder::Single))
            .expect("tetrahedral star bond should be valid");
    }
    builder
        .add_3d_conformer(vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        .expect("tetrahedral star coordinates should align");
    builder.build().expect("tetrahedral star should be a valid molecule")
}

fn prepared_tetrahedral_star(center: AtomSpec) -> Molecule {
    tetrahedral_star(center)
        .with_assigned_valence()
        .expect("tetrahedral star valence should be assignable")
}

#[test]
fn registry_declares_complete_assign_chiral_tags_operation_contract() {
    let operation = MOLECULE_OPS
        .iter()
        .copied()
        .find(|operation| operation.method == "with_chiral_tags_from_structure")
        .expect("chiral-tag assignment must be registered");

    assert_eq!(operation.impl_fn, "with_chiral_tags_from_structure_impl");
    assert_eq!(operation.domain, OperationDomain::Topology);
    assert_eq!(operation.kind, MoleculeOpKind::Weak);
    assert_eq!(operation.topology_edit, TopologyEditKind::Local);
    assert_eq!(operation.access.read(), BlockSet::COORDINATES);
    assert_eq!(
        operation.access.write(),
        BlockSet::TOPOLOGY
            .union(BlockSet::PROPERTIES)
            .union(BlockSet::DERIVED_CACHE)
    );
    assert_eq!(operation.may_mutate, operation.access.write());
    assert_eq!(operation.auto_remap, BlockSet::NONE);
    assert_eq!(operation.requires_mapping, MappingRequirement::None);
    assert_eq!(operation.parity, ParityPolicy::RequiredNow);
    assert!(operation.io_roundtrip);
    assert_eq!(
        operation.derived_effects.invalidate(),
        DerivedState::STEREO
            .union(DerivedState::DRAWING)
            .union(DerivedState::FINGERPRINT)
    );
    assert_eq!(operation.derived_effects.recompute(), DerivedState::NONE);
    assert_eq!(operation.derived_effects.preserve(), DerivedState::NONE);
    assert_eq!(
        operation.support,
        SupportStatus::SupportedWithRdkitParity {
            rdkit_version: "2026.03.1",
        }
    );

    let support = SUPPORT_MATRIX
        .iter()
        .find(|entry| {
            entry
                .operation
                .is_some_and(|candidate| candidate.method == operation.method)
        })
        .expect("registered operation must have a support entry");
    assert_eq!(support.feature.name, "stereo.perception");

    let invariant = OPERATION_INVARIANT_MATRIX
        .iter()
        .find(|entry| entry.operation.method == operation.method)
        .expect("registered operation must have an invariant entry");
    assert_eq!(invariant.profile, "weak_topology_state");

    let parity = PARITY_MATRIX
        .iter()
        .find(|entry| entry.operation.method == operation.method)
        .expect("required parity operation must have a parity entry");
    assert_eq!(parity.profile, "assign_atom_chiral_tags_from_structure_rdkit");
}

#[test]
fn value_and_in_place_forms_preserve_source_coordinates_and_unrelated_properties() {
    let molecule = prepared_tetrahedral_star(AtomSpec::new(Element::C).with_prop("center_keep", "yes"))
        .with_prop("_StereochemDone", "1")
        .with_prop("molecule_keep", "yes");
    let original = molecule.clone();
    let original_coordinates = molecule.conformers_3d().to_vec();
    let original_bonds = molecule.bonds().to_vec();
    let original_neighbors = molecule.atoms()[1..].to_vec();

    let assigned = molecule
        .with_chiral_tags_from_structure(-1, true)
        .expect("value-style assignment should succeed");

    assert_eq!(molecule, original, "value-style call mutated its source");
    assert_eq!(assigned.conformers_3d(), original_coordinates);
    assert_eq!(assigned.bonds(), original_bonds);
    assert_eq!(&assigned.atoms()[1..], original_neighbors);
    assert_eq!(assigned.atoms()[CENTER].chiral_tag(), ChiralTag::TetrahedralCcw);
    assert_eq!(assigned.atoms()[CENTER].prop("_NonExplicit3DChirality"), Some("1"));
    assert_eq!(assigned.atoms()[CENTER].prop("center_keep"), Some("yes"));
    assert_eq!(assigned.prop("_StereochemDone"), None);
    assert_eq!(assigned.prop("molecule_keep"), Some("yes"));

    let shared_source = molecule.clone();
    let mut in_place = molecule.clone();
    in_place
        .assign_chiral_tags_from_structure_(-1, true)
        .expect("in-place assignment should succeed");
    assert_eq!(shared_source, original, "COW mutation escaped into a shared clone");
    assert_eq!(in_place, assigned);

    in_place
        .with_assigned_valence()
        .expect("cache-dependent valence operation should remain usable after invalidation");
}

#[test]
fn source_defined_noop_paths_and_replace_false_preserve_state() {
    let no_conformer = Molecule::from_smiles("F[C@](Cl)(Br)I").expect("tagged no-conformer molecule should parse");
    let no_conformer_result = no_conformer
        .with_chiral_tags_from_structure(-1, true)
        .expect("no conformer is a source-defined no-op");
    assert_eq!(no_conformer_result, no_conformer);

    let prepared = prepared_tetrahedral_star(AtomSpec::new(Element::C));
    let coordinates = prepared.conformers_3d()[0].coordinates().to_vec();
    let non_3d = prepared
        .with_only_3d_conformer(coordinates, false)
        .expect("non-3D conformer should be representable");
    let non_3d_result = non_3d
        .with_chiral_tags_from_structure(-1, true)
        .expect("non-3D conformer is a source-defined no-op");
    assert_eq!(non_3d_result, non_3d);

    let existing = prepared_tetrahedral_star(AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::TetrahedralCw));
    let preserve_existing = existing
        .with_chiral_tags_from_structure(-1, false)
        .expect("replace=false assignment should succeed");
    assert_eq!(preserve_existing, existing);
}

#[test]
fn structured_errors_are_transactional_for_value_and_in_place_forms() {
    let prepared = prepared_tetrahedral_star(AtomSpec::new(Element::C)).with_prop("_StereochemDone", "1");
    let original = prepared.clone();
    let error = prepared
        .with_chiral_tags_from_structure(17, true)
        .expect_err("missing conformer id must fail explicitly");
    assert!(matches!(
        error,
        OperationError::Stereo {
            operation,
            source: StereoError::ConformerNotFound { conformer_id: 17 },
        } if operation.method == "with_chiral_tags_from_structure"
    ));
    assert_eq!(prepared, original);

    let mut in_place = prepared.clone();
    let error = in_place
        .assign_chiral_tags_from_structure_(17, true)
        .expect_err("in-place missing conformer id must fail explicitly");
    assert!(matches!(
        error,
        OperationError::Stereo {
            source: StereoError::ConformerNotFound { conformer_id: 17 },
            ..
        }
    ));
    assert_eq!(in_place, original, "failed in-place call committed partial state");

    let mut missing_valence = tetrahedral_star(AtomSpec::new(Element::C)).with_prop("_StereochemDone", "1");
    let missing_valence_original = missing_valence.clone();
    let error = missing_valence
        .assign_chiral_tags_from_structure_(-1, true)
        .expect_err("missing implicit-H state must fail explicitly");
    assert!(matches!(
        error,
        OperationError::Stereo {
            source: StereoError::MissingImplicitHydrogenState,
            ..
        }
    ));
    assert_eq!(
        missing_valence, missing_valence_original,
        "kernel failure committed topology or property changes"
    );
}
