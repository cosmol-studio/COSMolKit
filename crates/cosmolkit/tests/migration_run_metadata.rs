use cosmolkit::{
    BINDING_CONTRACT, BindingItem, BindingOwner, BindingTypeRole, FeatureSpec, FunctionStatus,
    MOLECULE_OPS, OPERATION_INVARIANT_MATRIX, PARITY_MATRIX, ParityPolicy, SUPPORT_MATRIX,
    StateModel, UnsupportedFeatureError, feature_spec, feature_specs, operation_invariant,
    operation_invariant_matrix, operation_parity, operation_spec, operation_specs, parity_matrix,
    support_matrix, version,
};

fn binding_entry(semantic_id: &str) -> &'static cosmolkit::BindingContractEntry {
    BINDING_CONTRACT
        .iter()
        .find(|entry| entry.semantic_id == semantic_id)
        .unwrap_or_else(|| panic!("missing binding contract entry {semantic_id}"))
}

fn expected_feature_names() -> Vec<&'static str> {
    let mut expected = Vec::new();
    if cfg!(feature = "cap-stereoisomers") {
        expected.push("cap-stereoisomers");
    }
    if cfg!(feature = "cap-reaction") {
        expected.push("cap-reaction");
    }
    if cfg!(feature = "cap-alignment") {
        expected.push("cap-alignment");
    }
    if cfg!(feature = "cap-conformer") {
        expected.push("cap-conformer");
    }
    if cfg!(feature = "cap-tautomer") {
        expected.push("cap-tautomer");
    }
    if cfg!(feature = "cap-forcefields") {
        expected.push("cap-forcefields");
    }
    if cfg!(feature = "cap-fingerprints") {
        expected.push("cap-fingerprints");
    }
    if cfg!(feature = "cap-sanitize") {
        expected.push("cap-sanitize");
    }
    if cfg!(feature = "cap-kekulize") {
        expected.push("cap-kekulize");
    }
    if cfg!(feature = "cap-aromaticity") {
        expected.push("cap-aromaticity");
    }
    if cfg!(feature = "cap-valence") {
        expected.push("cap-valence");
    }
    if cfg!(feature = "cap-radicals") {
        expected.push("cap-radicals");
    }
    if cfg!(feature = "cap-rings") {
        expected.push("cap-rings");
    }
    if cfg!(feature = "cap-stereo") {
        expected.push("cap-stereo");
    }
    if cfg!(feature = "cap-hydrogens") {
        expected.push("cap-hydrogens");
    }
    if cfg!(feature = "cap-transforms") {
        expected.push("cap-transforms");
    }
    if cfg!(feature = "cap-depict") {
        expected.push("cap-depict");
    }
    expected
}

fn expected_operation_methods() -> Vec<&'static str> {
    let mut expected = Vec::new();
    if cfg!(feature = "cap-stereoisomers") {
        expected.extend([
            "enumerate_stereoisomers_with_options",
            "enumerate_stereoisomers_with_random_bits",
        ]);
    }
    if cfg!(feature = "cap-reaction") {
        expected.extend([
            "reaction_products_with_params",
            "run",
            "reaction_products_from_inputs",
            "apply_reaction_with_params",
        ]);
    }
    if cfg!(feature = "cap-alignment") {
        expected.extend([
            "with_alignment_to_with_params",
            "with_aligned_conformers_with_params",
        ]);
    }
    if cfg!(feature = "cap-conformer") {
        expected.extend([
            "with_3d_conformer_with_params",
            "with_3d_conformer_result_with_params",
            "with_3d_conformers_with_params",
            "with_3d_conformers_result_with_params",
        ]);
    }
    if cfg!(feature = "cap-tautomer") {
        expected.extend([
            "enumerate_tautomers_with_params",
            "canonical_tautomer_with_params",
        ]);
    }
    if cfg!(feature = "cap-forcefields") {
        expected.push("with_mmff_optimized_with_params");
        expected.push("with_mmff_optimized_confs_with_params");
    }
    if cfg!(feature = "cap-fingerprints") {
        expected.push("with_atom_pair_atom_code");
    }
    if cfg!(feature = "cap-forcefields") {
        expected.push("with_uff_optimized_with_params");
        expected.push("with_uff_optimized_confs_with_params");
    }
    if cfg!(feature = "cap-sanitize") {
        expected.push("sanitize_with_params");
    }
    if cfg!(feature = "cap-kekulize") {
        expected.push("with_kekulized_bonds_with_params");
    }
    if cfg!(feature = "cap-aromaticity") {
        expected.push("with_assigned_aromaticity_with_params");
    }
    if cfg!(feature = "cap-valence") {
        expected.push("with_assigned_valence_with_params");
    }
    if cfg!(feature = "cap-radicals") {
        expected.push("with_assigned_radicals");
    }
    if cfg!(feature = "cap-rings") {
        expected.extend([
            "with_assigned_rings",
            "with_assigned_ring_families_with_params",
        ]);
    }
    if cfg!(feature = "cap-stereo") {
        expected.extend([
            "with_chiral_tags_from_structure_with_params",
            "potential_stereo_with_params",
            "with_cip_labels_with_options",
        ]);
    }
    if cfg!(feature = "cap-hydrogens") {
        expected.extend([
            "with_hydrogens_with_params",
            "without_hydrogens_with_params",
        ]);
    }
    if cfg!(feature = "cap-transforms") {
        expected.push("with_atom_position_with_params");
        expected.extend([
            "with_2d_coordinate_block_with_params",
            "with_3d_coordinates_with_params",
            "with_added_3d_conformer_with_params",
            "with_only_3d_conformer_with_params",
            "with_cleared_3d_conformers",
        ]);
    }
    if cfg!(feature = "cap-depict") {
        expected.push("with_2d_coordinates_with_params");
    }
    expected
}

// The unchanged generator omits only explicitly NotApplicable declarations.
// Keep the full literal operation fixture and every original parity row/order.
fn expected_parity_methods() -> Vec<&'static str> {
    expected_operation_methods()
        .into_iter()
        .filter(|method| {
            !matches!(
                *method,
                "with_2d_coordinate_block_with_params"
                    | "with_3d_coordinates_with_params"
                    | "with_added_3d_conformer_with_params"
                    | "with_only_3d_conformer_with_params"
                    | "with_cleared_3d_conformers"
            )
        })
        .collect()
}

#[test]
fn status_and_parity_values_preserve_every_public_branch() {
    let statuses = [
        FunctionStatus::Parity { reference: "RDKit" },
        FunctionStatus::ParityWithDifferences {
            reference: "Gemmi",
            explanation: "Approved difference for the documented boundary",
        },
        FunctionStatus::Native,
        FunctionStatus::Experimental,
    ];
    assert_eq!(statuses.len(), 4);
    assert!(statuses.iter().enumerate().all(|(i, status)| {
        statuses
            .iter()
            .enumerate()
            .all(|(j, candidate)| (candidate == status) == (i == j))
    }));

    let policies = [
        ParityPolicy::NotApplicable,
        ParityPolicy::RequiredWhenSupported,
        ParityPolicy::RequiredNow,
    ];
    assert_eq!(policies.len(), 3);
    assert_ne!(policies[0], policies[1]);
    assert_ne!(policies[1], policies[2]);

    let unsupported = UnsupportedFeatureError {
        feature: "test-unsupported",
        reason: "not implemented",
    };
    assert_eq!(unsupported.feature, "test-unsupported");
    assert_eq!(unsupported.reason, "not implemented");
}

#[test]
fn slice_accessors_return_the_exact_generated_tables() {
    assert!(core::ptr::eq(operation_specs(), MOLECULE_OPS));
    assert!(core::ptr::eq(support_matrix(), SUPPORT_MATRIX));
    assert!(core::ptr::eq(
        operation_invariant_matrix(),
        OPERATION_INVARIANT_MATRIX
    ));
    assert!(core::ptr::eq(parity_matrix(), PARITY_MATRIX));
}

#[test]
fn metadata_binding_contract_rows_are_complete_and_read_only() {
    for (semantic_id, role) in [
        ("types.FunctionStatus", BindingTypeRole::Value),
        ("types.ParityPolicy", BindingTypeRole::Value),
        ("types.FeatureSpec", BindingTypeRole::Value),
        ("types.MoleculeOpSpec", BindingTypeRole::Result),
        ("types.SupportMatrixEntry", BindingTypeRole::Result),
        ("types.OperationInvariantEntry", BindingTypeRole::Result),
        ("types.ParityMatrixEntry", BindingTypeRole::Result),
    ] {
        let entry = binding_entry(semantic_id);
        assert_eq!(entry.item, BindingItem::Type);
        assert_eq!(entry.owner, BindingOwner::Type);
        assert_eq!(entry.feature, "metadata");
        assert_eq!(entry.status, FunctionStatus::Experimental);
        assert_eq!(entry.type_role, Some(role));
    }

    let iterator = binding_entry("types.FeatureSpecIter");
    assert_eq!(iterator.item, BindingItem::Type);
    assert_eq!(iterator.owner, BindingOwner::Type);
    assert_eq!(iterator.status, FunctionStatus::Experimental);
    assert_eq!(iterator.type_role, Some(BindingTypeRole::Result));

    for (semantic_id, python, javascript) in [
        ("module.feature_specs", "feature_specs", "featureSpecs"),
        ("module.feature_spec", "feature_spec", "featureSpec"),
        (
            "module.operation_specs",
            "operation_specs",
            "operationSpecs",
        ),
        ("module.operation_spec", "operation_spec", "operationSpec"),
        ("module.support_matrix", "support_matrix", "supportMatrix"),
        (
            "module.operation_invariant_matrix",
            "operation_invariant_matrix",
            "operationInvariantMatrix",
        ),
        ("module.parity_matrix", "parity_matrix", "parityMatrix"),
        (
            "module.operation_invariant",
            "operation_invariant",
            "operationInvariant",
        ),
        (
            "module.operation_parity",
            "operation_parity",
            "operationParity",
        ),
        ("module.version", "version", "version"),
    ] {
        let entry = binding_entry(semantic_id);
        assert_eq!(entry.item, BindingItem::Callable);
        assert_eq!(entry.owner, BindingOwner::Module);
        assert_eq!(entry.python_name, python);
        assert_eq!(entry.javascript_name, javascript);
        assert_eq!(entry.feature, "metadata");
        assert_eq!(entry.status, FunctionStatus::Experimental);
        let callable = entry.callable.expect("module callable payload");
        assert_eq!(callable.state_model, StateModel::ReadOnly);
        assert_eq!(callable.operation_semantic_id, None);
    }
}

#[test]
fn generated_tables_have_one_source_and_queries_do_not_define_parallel_rows() {
    let metadata_source = include_str!("../src/ops/metadata.rs");
    let registry_source = include_str!("../src/ops/registry.rs");
    let generator_source = include_str!("../../cosmolkit-macros/src/matrices.rs");

    for table in [
        "MOLECULE_OPS",
        "SUPPORT_MATRIX",
        "OPERATION_INVARIANT_MATRIX",
        "PARITY_MATRIX",
    ] {
        assert!(generator_source.contains(&format!("pub static {table}:")));
        assert!(!metadata_source.contains(&format!("pub const {table}:")));
        assert!(metadata_source.contains(&format!("super::runtime::registry::{table}")));
        assert!(!metadata_source.contains(&format!("super::registry::{table}")));
    }
    assert_eq!(registry_source.matches("molecule_ops!").count(), 1);
}

#[cfg(not(feature = "cap-hydrogens"))]
#[test]
fn default_configuration_exposes_empty_generated_metadata_and_exact_misses() {
    let expected_features = expected_feature_names();
    let expected_methods = expected_operation_methods();
    let features = feature_specs().collect::<Vec<_>>();
    let operations = operation_specs();
    assert_eq!(
        features
            .iter()
            .map(|feature| feature.name)
            .collect::<Vec<_>>(),
        expected_features
    );
    assert_eq!(
        operations
            .iter()
            .map(|operation| operation.method)
            .collect::<Vec<_>>(),
        expected_methods
    );
    assert_eq!(support_matrix().len(), expected_methods.len());
    assert_eq!(operation_invariant_matrix().len(), expected_methods.len());
    assert_eq!(
        parity_matrix()
            .iter()
            .map(|row| row.operation.method)
            .collect::<Vec<_>>(),
        expected_parity_methods()
    );

    if expected_methods.is_empty() {
        assert!(feature_specs().next().is_none());
        assert!(operation_specs().is_empty());
        assert!(support_matrix().is_empty());
        assert!(operation_invariant_matrix().is_empty());
        assert!(parity_matrix().is_empty());
    }
    for feature in features {
        assert!(core::ptr::eq(feature_spec(feature.name).unwrap(), feature));
    }
    for (index, operation) in operations.iter().enumerate() {
        let support = &support_matrix()[index];
        let invariant = &operation_invariant_matrix()[index];

        assert!(core::ptr::eq(
            operation_spec(operation.method).unwrap(),
            *operation
        ));
        assert!(core::ptr::eq(support.operation.unwrap(), *operation));
        assert!(core::ptr::eq(invariant.operation, *operation));

        assert!(core::ptr::eq(
            operation_invariant(operation.method).unwrap(),
            invariant
        ));
        if let Some(parity_index) = expected_parity_methods()
            .iter()
            .position(|method| *method == operation.method)
        {
            let parity = &parity_matrix()[parity_index];
            assert!(core::ptr::eq(parity.operation, *operation));
            assert!(core::ptr::eq(
                operation_parity(operation.method).unwrap(),
                parity
            ));
            assert!(core::ptr::eq(support.feature, parity.feature));
        } else {
            assert_eq!(operation_parity(operation.method), None);
            assert_eq!(operation.parity, ParityPolicy::NotApplicable);
            assert_eq!(operation.status, FunctionStatus::Native);
        }
        assert!(core::ptr::eq(
            feature_spec(support.feature.name).unwrap(),
            support.feature
        ));
    }
    #[cfg(feature = "cap-depict")]
    {
        let operation = operation_spec("with_2d_coordinates_with_params").unwrap();
        let index = expected_methods
            .iter()
            .position(|method| *method == operation.method)
            .unwrap();
        assert_eq!(operation.parity, ParityPolicy::RequiredWhenSupported);
        assert_eq!(
            operation_invariant_matrix()[index].profile,
            "coordinate_2d_layout"
        );
        let parity_index = expected_parity_methods()
            .iter()
            .position(|method| *method == operation.method)
            .unwrap();
        assert_eq!(
            parity_matrix()[parity_index].profile,
            "compute_2d_coordinates_rdkit"
        );
        assert!(core::ptr::eq(operations[index], operation));
        assert!(core::ptr::eq(
            support_matrix()[index].operation.unwrap(),
            operation
        ));
        assert!(core::ptr::eq(
            operation_invariant_matrix()[index].operation,
            operation
        ));
        assert!(core::ptr::eq(
            parity_matrix()[parity_index].operation,
            operation
        ));
    }

    for name in ["", "cap-hydrogens", "Cap-hydrogens", "unknown"] {
        assert_eq!(feature_spec(name), None);
    }
    for method in ["", "with_hydrogens", "WITH_HYDROGENS", "unknown"] {
        assert_eq!(operation_spec(method), None);
        assert_eq!(operation_invariant(method), None);
        assert_eq!(operation_parity(method), None);
    }
    assert_eq!(version(), env!("CARGO_PKG_VERSION"));
}

#[cfg(feature = "cap-hydrogens")]
#[test]
fn hydrogens_configuration_preserves_order_profiles_and_pointer_identity() {
    let features = feature_specs().collect::<Vec<_>>();
    assert_eq!(
        features
            .iter()
            .map(|feature| feature.name)
            .collect::<Vec<_>>(),
        expected_feature_names()
    );
    let hydrogens = feature_spec("cap-hydrogens").expect("hydrogens feature");
    assert_eq!(hydrogens.category, "chemistry");
    assert!(!hydrogens.docs.is_empty());

    let operations = operation_specs();
    assert_eq!(
        operations
            .iter()
            .map(|operation| operation.method)
            .collect::<Vec<_>>(),
        expected_operation_methods()
    );
    assert_eq!(support_matrix().len(), operations.len());
    assert_eq!(operation_invariant_matrix().len(), operations.len());
    assert_eq!(
        parity_matrix()
            .iter()
            .map(|row| row.operation.method)
            .collect::<Vec<_>>(),
        expected_parity_methods()
    );

    for (index, operation) in operations.iter().enumerate() {
        assert!(core::ptr::eq(
            operation_spec(operation.method).expect("operation lookup"),
            *operation
        ));
        assert!(core::ptr::eq(
            support_matrix()[index]
                .operation
                .expect("support operation"),
            *operation
        ));
        assert!(core::ptr::eq(
            operation_invariant_matrix()[index].operation,
            *operation
        ));

        assert!(core::ptr::eq(
            feature_spec(support_matrix()[index].feature.name).expect("feature lookup"),
            support_matrix()[index].feature
        ));
        let expected_parity = if matches!(
            operation.method,
            "with_2d_coordinates_with_params"
                | "with_3d_conformer_with_params"
                | "with_3d_conformer_result_with_params"
                | "with_3d_conformers_with_params"
                | "with_3d_conformers_result_with_params"
        ) {
            ParityPolicy::RequiredWhenSupported
        } else if matches!(
            operation.method,
            "with_2d_coordinate_block_with_params"
                | "with_3d_coordinates_with_params"
                | "with_added_3d_conformer_with_params"
                | "with_only_3d_conformer_with_params"
                | "with_cleared_3d_conformers"
        ) {
            ParityPolicy::NotApplicable
        } else {
            ParityPolicy::RequiredNow
        };
        assert_eq!(operation.parity, expected_parity);
        assert_eq!(
            operation_invariant(operation.method).expect("invariant lookup"),
            &operation_invariant_matrix()[index]
        );
        if let Some(parity_index) = expected_parity_methods()
            .iter()
            .position(|method| *method == operation.method)
        {
            let parity = &parity_matrix()[parity_index];
            assert!(core::ptr::eq(parity.operation, *operation));
            assert!(core::ptr::eq(
                support_matrix()[index].feature,
                parity.feature
            ));
            assert_eq!(
                operation_parity(operation.method).expect("parity lookup"),
                parity
            );
        } else {
            assert_eq!(operation_parity(operation.method), None);
            assert_eq!(operation.status, FunctionStatus::Native);
        }
    }

    let add_index = operations
        .iter()
        .position(|operation| operation.method == "with_hydrogens_with_params")
        .expect("add-hydrogen operation");
    let remove_index = operations
        .iter()
        .position(|operation| operation.method == "without_hydrogens_with_params")
        .expect("remove-hydrogen operation");
    assert_eq!(
        operation_invariant_matrix()[add_index].profile,
        "strong_topology_with_coordinates"
    );
    assert_eq!(
        operation_invariant_matrix()[remove_index].profile,
        "strong_topology_with_coordinates"
    );
    assert_eq!(parity_matrix()[add_index].profile, "add_hydrogens_rdkit");
    assert_eq!(
        parity_matrix()[remove_index].profile,
        "remove_hydrogens_rdkit"
    );
    assert_eq!(parity_matrix()[add_index].rdkit_version, None);
    assert_eq!(parity_matrix()[remove_index].rdkit_version, None);
}

#[cfg(feature = "cap-hydrogens")]
#[test]
fn hydrogens_lookups_are_exact_and_reject_unknown_or_wrong_case_names() {
    let generated_hydrogens = feature_specs()
        .find(|feature| feature.name == "cap-hydrogens")
        .expect("feature iterator must contain hydrogens");
    assert!(core::ptr::eq(
        feature_spec("cap-hydrogens").expect("feature lookup"),
        generated_hydrogens
    ));
    for name in ["", "Hydrogens", "HYDROGENS", "unknown"] {
        assert_eq!(feature_spec(name), None);
    }
    for method in ["", "With_Hydrogens", "WITH_HYDROGENS", "unknown"] {
        assert_eq!(operation_spec(method), None);
        assert_eq!(operation_invariant(method), None);
        assert_eq!(operation_parity(method), None);
    }
    assert_eq!(version(), env!("CARGO_PKG_VERSION"));
}
