use cosmolkit::{
    BINDING_CONTRACT, BindingExposure, BindingItem, BindingOwner, BindingParity, BindingSupport,
    BindingTypeRole, FeatureSpec, MOLECULE_OPS, OPERATION_INVARIANT_MATRIX, PARITY_MATRIX,
    ParityPolicy, SUPPORT_MATRIX, StateModel, SupportStatus, UnsupportedFeatureError, feature_spec,
    feature_specs, operation_invariant, operation_invariant_matrix, operation_parity,
    operation_spec, operation_specs, parity_matrix, support_matrix, version,
};

fn binding_entry(semantic_id: &str) -> &'static cosmolkit::BindingContractEntry {
    BINDING_CONTRACT
        .iter()
        .find(|entry| entry.semantic_id == semantic_id)
        .unwrap_or_else(|| panic!("missing binding contract entry {semantic_id}"))
}

fn expected_feature_names() -> Vec<&'static str> {
    let mut expected = Vec::new();
    if cfg!(feature = "sanitize") {
        expected.push("sanitize");
    }
    if cfg!(feature = "kekulize") {
        expected.push("kekulize");
    }
    if cfg!(feature = "aromaticity") {
        expected.push("aromaticity");
    }
    if cfg!(feature = "valence") {
        expected.push("valence");
    }
    if cfg!(feature = "radicals") {
        expected.push("radicals");
    }
    if cfg!(feature = "rings") {
        expected.push("rings");
    }
    if cfg!(feature = "stereo") {
        expected.push("stereo");
    }
    if cfg!(feature = "hydrogens") {
        expected.push("hydrogens");
    }
    if cfg!(feature = "transforms") {
        expected.push("transforms");
    }
    expected
}

fn expected_operation_methods() -> Vec<&'static str> {
    let mut expected = Vec::new();
    if cfg!(feature = "sanitize") {
        expected.push("sanitize_with_params");
    }
    if cfg!(feature = "kekulize") {
        expected.push("with_kekulized_bonds_with_params");
    }
    if cfg!(feature = "aromaticity") {
        expected.push("with_assigned_aromaticity_with_params");
    }
    if cfg!(feature = "valence") {
        expected.push("with_assigned_valence_with_params");
    }
    if cfg!(feature = "radicals") {
        expected.push("with_assigned_radicals");
    }
    if cfg!(feature = "rings") {
        expected.extend([
            "with_assigned_rings",
            "with_assigned_ring_families_with_params",
        ]);
    }
    if cfg!(feature = "stereo") {
        expected.extend([
            "with_chiral_tags_from_structure_with_params",
            "potential_stereo_with_params",
        ]);
    }
    if cfg!(feature = "hydrogens") {
        expected.extend([
            "with_hydrogens_with_params",
            "without_hydrogens_with_params",
        ]);
    }
    if cfg!(feature = "transforms") {
        expected.push("with_atom_position_with_params");
    }
    expected
}

#[test]
fn status_and_parity_values_preserve_every_public_branch() {
    let statuses = [
        SupportStatus::Supported,
        SupportStatus::SupportedWithRdkitParity,
        SupportStatus::PreservedOnly,
        SupportStatus::Experimental,
        SupportStatus::Unsupported { reason: "missing" },
    ];
    assert_eq!(statuses.len(), 5);
    assert!(
        statuses
            .iter()
            .all(|status| statuses.iter().any(|candidate| candidate == status))
    );

    let policies = [
        ParityPolicy::NotApplicable,
        ParityPolicy::RequiredWhenSupported,
        ParityPolicy::RequiredNow,
    ];
    assert_eq!(policies.len(), 3);
    assert_ne!(policies[0], policies[1]);
    assert_ne!(policies[1], policies[2]);

    static UNSUPPORTED: FeatureSpec = FeatureSpec {
        name: "test-unsupported",
        category: "test",
        status: SupportStatus::Unsupported {
            reason: "not implemented",
        },
        rdkit_parity_sensitive: false,
        docs: "Test-only local value.",
    };
    assert_eq!(
        UnsupportedFeatureError::from_spec(&UNSUPPORTED),
        UnsupportedFeatureError {
            feature: "test-unsupported",
            reason: "not implemented",
        }
    );
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
        ("types.SupportStatus", BindingTypeRole::Value),
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
        assert_eq!(entry.exposure, BindingExposure::Public);
        assert_eq!(entry.support, BindingSupport::Supported);
        assert_eq!(entry.parity, BindingParity::NotApplicable);
        assert_eq!(entry.type_role, Some(role));
    }

    let iterator = binding_entry("types.FeatureSpecIter");
    assert_eq!(iterator.item, BindingItem::Type);
    assert_eq!(iterator.owner, BindingOwner::Type);
    assert_eq!(iterator.exposure, BindingExposure::Registered);
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
        assert_eq!(entry.support, BindingSupport::Supported);
        assert_eq!(entry.parity, BindingParity::NotApplicable);
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
        assert!(generator_source.contains(&format!("pub const {table}:")));
        assert!(!metadata_source.contains(&format!("pub const {table}:")));
        assert!(metadata_source.contains(&format!("super::runtime::registry::{table}")));
        assert!(!metadata_source.contains(&format!("super::registry::{table}")));
    }
    assert_eq!(registry_source.matches("molecule_ops!").count(), 1);
}

#[cfg(not(feature = "hydrogens"))]
#[test]
fn default_configuration_exposes_empty_generated_metadata_and_exact_misses() {
    assert!(feature_specs().next().is_none());
    assert!(operation_specs().is_empty());
    assert!(support_matrix().is_empty());
    assert!(operation_invariant_matrix().is_empty());
    assert!(parity_matrix().is_empty());

    for name in ["", "hydrogens", "Hydrogens", "unknown"] {
        assert_eq!(feature_spec(name), None);
    }
    for method in ["", "with_hydrogens", "WITH_HYDROGENS", "unknown"] {
        assert_eq!(operation_spec(method), None);
        assert_eq!(operation_invariant(method), None);
        assert_eq!(operation_parity(method), None);
    }
    assert_eq!(version(), env!("CARGO_PKG_VERSION"));
}

#[cfg(feature = "hydrogens")]
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
    let hydrogens = feature_spec("hydrogens").expect("hydrogens feature");
    assert_eq!(hydrogens.category, "chemistry");
    assert!(hydrogens.rdkit_parity_sensitive);
    assert_eq!(hydrogens.status, SupportStatus::SupportedWithRdkitParity);

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
    assert_eq!(parity_matrix().len(), operations.len());

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
        assert!(core::ptr::eq(parity_matrix()[index].operation, *operation));
        assert!(core::ptr::eq(
            support_matrix()[index].feature,
            parity_matrix()[index].feature
        ));
        assert!(core::ptr::eq(
            feature_spec(support_matrix()[index].feature.name).expect("feature lookup"),
            support_matrix()[index].feature
        ));
        assert_eq!(operation.parity, ParityPolicy::RequiredNow);
        assert_eq!(
            operation_invariant(operation.method).expect("invariant lookup"),
            &operation_invariant_matrix()[index]
        );
        assert_eq!(
            operation_parity(operation.method).expect("parity lookup"),
            &parity_matrix()[index]
        );
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

#[cfg(feature = "hydrogens")]
#[test]
fn hydrogens_lookups_are_exact_and_reject_unknown_or_wrong_case_names() {
    let generated_hydrogens = feature_specs()
        .find(|feature| feature.name == "hydrogens")
        .expect("feature iterator must contain hydrogens");
    assert!(core::ptr::eq(
        feature_spec("hydrogens").expect("feature lookup"),
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
