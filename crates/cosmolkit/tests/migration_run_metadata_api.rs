use cosmolkit::{
    BINDING_CONTRACT, BindingDefault, BindingExposure, BindingItem, BindingKind, BindingOwner,
    BindingParity, BindingSupport, BindingTypeRole, FeatureSpec, FeatureSpecIter, MOLECULE_OPS,
    MoleculeOpSpec, OPERATION_INVARIANT_MATRIX, OperationInvariantEntry, PARITY_MATRIX,
    ParityMatrixEntry, ParityPolicy, SUPPORT_MATRIX, StateModel, SupportMatrixEntry, SupportStatus,
    feature_spec, feature_specs, operation_invariant, operation_invariant_matrix, operation_parity,
    operation_spec, operation_specs, parity_matrix, support_matrix, version,
};

fn binding_entry(semantic_id: &str) -> &'static cosmolkit::BindingContractEntry {
    BINDING_CONTRACT
        .iter()
        .find(|entry| entry.semantic_id == semantic_id)
        .unwrap_or_else(|| panic!("missing binding contract entry {semantic_id}"))
}

fn compact(text: &str) -> String {
    text.chars()
        .filter(|character| !character.is_whitespace())
        .collect()
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
    if cfg!(feature = "depict") {
        expected.push("depict");
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
            "with_cip_labels_with_options",
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
    if cfg!(feature = "depict") {
        expected.push("with_2d_coordinates_with_params");
    }
    expected
}

#[test]
fn canonical_metadata_signatures_compile_from_the_public_crate() {
    let _: fn() -> FeatureSpecIter = feature_specs;
    let _: fn(&str) -> Option<&'static FeatureSpec> = feature_spec;
    let _: fn() -> &'static [&'static MoleculeOpSpec] = operation_specs;
    let _: fn(&str) -> Option<&'static MoleculeOpSpec> = operation_spec;
    let _: fn() -> &'static [SupportMatrixEntry] = support_matrix;
    let _: fn() -> &'static [OperationInvariantEntry] = operation_invariant_matrix;
    let _: fn() -> &'static [ParityMatrixEntry] = parity_matrix;
    let _: fn(&str) -> Option<&'static OperationInvariantEntry> = operation_invariant;
    let _: fn(&str) -> Option<&'static ParityMatrixEntry> = operation_parity;
    let _: fn() -> &'static str = version;

    assert_eq!(version(), env!("CARGO_PKG_VERSION"));
    assert_ne!(SupportStatus::Supported, SupportStatus::Experimental);
    assert_ne!(
        ParityPolicy::RequiredWhenSupported,
        ParityPolicy::RequiredNow
    );
}

#[test]
fn metadata_type_rows_have_exact_public_roles_and_paths() {
    let expected = [
        (
            "types.SupportStatus",
            "crate::SupportStatus",
            BindingTypeRole::Value,
        ),
        (
            "types.ParityPolicy",
            "crate::ParityPolicy",
            BindingTypeRole::Value,
        ),
        (
            "types.FeatureSpec",
            "crate::FeatureSpec",
            BindingTypeRole::Value,
        ),
        (
            "types.MoleculeOpSpec",
            "crate::MoleculeOpSpec",
            BindingTypeRole::Result,
        ),
        (
            "types.SupportMatrixEntry",
            "crate::SupportMatrixEntry",
            BindingTypeRole::Result,
        ),
        (
            "types.OperationInvariantEntry",
            "crate::OperationInvariantEntry",
            BindingTypeRole::Result,
        ),
        (
            "types.ParityMatrixEntry",
            "crate::ParityMatrixEntry",
            BindingTypeRole::Result,
        ),
    ];

    for (semantic_id, rust_path, role) in expected {
        let entry = binding_entry(semantic_id);
        assert_eq!(entry.item, BindingItem::Type);
        assert_eq!(entry.owner, BindingOwner::Type);
        assert_eq!(compact(entry.rust_path), compact(rust_path));
        assert_eq!(entry.feature, "metadata");
        assert_eq!(entry.exposure, BindingExposure::Public);
        assert_eq!(entry.support, BindingSupport::Supported);
        assert_eq!(entry.parity, BindingParity::NotApplicable);
        assert_eq!(entry.type_role, Some(role));
        assert_eq!(entry.callable, None);
    }

    let iterator = binding_entry("types.FeatureSpecIter");
    assert_eq!(iterator.item, BindingItem::Type);
    assert_eq!(iterator.owner, BindingOwner::Type);
    assert_eq!(compact(iterator.rust_path), "crate::FeatureSpecIter");
    assert_eq!(iterator.feature, "metadata");
    assert_eq!(iterator.exposure, BindingExposure::Registered);
    assert_eq!(iterator.type_role, Some(BindingTypeRole::Result));
}

#[test]
fn metadata_callable_rows_exactly_match_names_signatures_and_defaults() {
    let expected = [
        (
            "module.feature_specs",
            "feature_specs",
            "featureSpecs",
            "crate::FeatureSpecIter",
            None,
            BindingExposure::Registered,
        ),
        (
            "module.feature_spec",
            "feature_spec",
            "featureSpec",
            "Option<&'staticcrate::FeatureSpec>",
            Some(("name", "&str")),
            BindingExposure::Registered,
        ),
        (
            "module.operation_specs",
            "operation_specs",
            "operationSpecs",
            "&'static[&'staticcrate::MoleculeOpSpec]",
            None,
            BindingExposure::Registered,
        ),
        (
            "module.operation_spec",
            "operation_spec",
            "operationSpec",
            "Option<&'staticcrate::MoleculeOpSpec>",
            Some(("method", "&str")),
            BindingExposure::Registered,
        ),
        (
            "module.support_matrix",
            "support_matrix",
            "supportMatrix",
            "&'static[crate::SupportMatrixEntry]",
            None,
            BindingExposure::Registered,
        ),
        (
            "module.operation_invariant_matrix",
            "operation_invariant_matrix",
            "operationInvariantMatrix",
            "&'static[crate::OperationInvariantEntry]",
            None,
            BindingExposure::Registered,
        ),
        (
            "module.parity_matrix",
            "parity_matrix",
            "parityMatrix",
            "&'static[crate::ParityMatrixEntry]",
            None,
            BindingExposure::Registered,
        ),
        (
            "module.operation_invariant",
            "operation_invariant",
            "operationInvariant",
            "Option<&'staticcrate::OperationInvariantEntry>",
            Some(("method", "&str")),
            BindingExposure::Registered,
        ),
        (
            "module.operation_parity",
            "operation_parity",
            "operationParity",
            "Option<&'staticcrate::ParityMatrixEntry>",
            Some(("method", "&str")),
            BindingExposure::Registered,
        ),
        (
            "module.version",
            "version",
            "version",
            "&'staticstr",
            None,
            BindingExposure::Public,
        ),
    ];

    for (semantic_id, rust_name, javascript_name, output, parameter, exposure) in expected {
        let entry = binding_entry(semantic_id);
        let callable = entry.callable.expect("metadata row must be callable");
        assert_eq!(entry.item, BindingItem::Callable);
        assert_eq!(entry.owner, BindingOwner::Module);
        assert_eq!(compact(entry.rust_path), format!("crate::{rust_name}"));
        assert_eq!(entry.python_name, rust_name);
        assert_eq!(entry.javascript_name, javascript_name);
        assert_eq!(entry.feature, "metadata");
        assert_eq!(entry.exposure, exposure);
        assert_eq!(entry.support, BindingSupport::Supported);
        assert_eq!(entry.parity, BindingParity::NotApplicable);
        assert_eq!(callable.kind, BindingKind::Module);
        assert_eq!(compact(callable.output_type), compact(output));
        assert_eq!(callable.error_type, None);
        assert_eq!(callable.state_model, StateModel::ReadOnly);
        assert_eq!(callable.operation_semantic_id, None);
        match parameter {
            Some((name, type_name)) => {
                assert_eq!(callable.parameters.len(), 1);
                assert_eq!(callable.parameters[0].name, name);
                assert_eq!(
                    compact(callable.parameters[0].type_name),
                    compact(type_name)
                );
                assert_eq!(callable.parameters[0].default, BindingDefault::Required);
            }
            None => assert!(callable.parameters.is_empty()),
        }
    }
}

#[test]
fn public_queries_preserve_generated_identity_and_fail_closed_misses() {
    assert!(core::ptr::eq(operation_specs(), MOLECULE_OPS));
    assert!(core::ptr::eq(support_matrix(), SUPPORT_MATRIX));
    assert!(core::ptr::eq(
        operation_invariant_matrix(),
        OPERATION_INVARIANT_MATRIX
    ));
    assert!(core::ptr::eq(parity_matrix(), PARITY_MATRIX));

    assert_eq!(
        feature_specs()
            .map(|feature| feature.name)
            .collect::<Vec<_>>(),
        expected_feature_names()
    );
    assert_eq!(
        operation_specs()
            .iter()
            .map(|operation| operation.method)
            .collect::<Vec<_>>(),
        expected_operation_methods()
    );
    assert_eq!(support_matrix().len(), operation_specs().len());
    assert_eq!(operation_invariant_matrix().len(), operation_specs().len());
    assert_eq!(parity_matrix().len(), operation_specs().len());

    for (index, operation) in operation_specs().iter().enumerate() {
        assert!(core::ptr::eq(
            operation_spec(operation.method).expect("registered operation"),
            *operation
        ));
        assert!(core::ptr::eq(
            support_matrix()[index]
                .operation
                .expect("support operation"),
            *operation
        ));
        assert!(core::ptr::eq(
            operation_invariant(operation.method)
                .expect("invariant row")
                .operation,
            *operation
        ));
        assert!(core::ptr::eq(
            operation_parity(operation.method)
                .expect("parity row")
                .operation,
            *operation
        ));
    }

    if cfg!(feature = "hydrogens") {
        let hydrogens = feature_spec("hydrogens").expect("hydrogen feature");
        assert!(core::ptr::eq(
            hydrogens,
            feature_specs()
                .find(|feature| feature.name == "hydrogens")
                .expect("hydrogen feature iterator row")
        ));
        for method in [
            "with_hydrogens_with_params",
            "without_hydrogens_with_params",
        ] {
            assert_eq!(operation_spec(method).unwrap().method, method);
            assert_eq!(
                operation_invariant(method).unwrap().operation.method,
                method
            );
            assert_eq!(operation_parity(method).unwrap().operation.method, method);
        }
    } else {
        assert_eq!(feature_spec("hydrogens"), None);
        assert_eq!(operation_spec("with_hydrogens_with_params"), None);
        assert_eq!(operation_invariant("with_hydrogens_with_params"), None);
        assert_eq!(operation_parity("with_hydrogens_with_params"), None);
    }
    for name in ["", "Hydrogens", "unknown"] {
        assert_eq!(feature_spec(name), None);
    }
    for method in ["", "WITH_HYDROGENS", "unknown"] {
        assert_eq!(operation_spec(method), None);
        assert_eq!(operation_invariant(method), None);
        assert_eq!(operation_parity(method), None);
    }
}

#[test]
fn metadata_projection_has_no_transaction_or_parallel_registry_escape() {
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
    for forbidden in [
        "OpParts",
        "MultiOutputOpParts",
        "MoleculeBuilder",
        "begin_topology_mut",
        "commit_topology",
        "derived_cache",
        "Arc<",
    ] {
        assert!(
            !metadata_source.contains(forbidden),
            "metadata must not expose runtime authority: {forbidden}"
        );
    }
}
