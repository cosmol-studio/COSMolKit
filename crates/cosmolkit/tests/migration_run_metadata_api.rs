use cosmolkit::{
    BINDING_CONTRACT, BindingDefault, BindingItem, BindingKind, BindingOwner, BindingTypeRole,
    FeatureSpec, FeatureSpecIter, FunctionStatus, MOLECULE_OPS, MoleculeOpSpec,
    OPERATION_INVARIANT_MATRIX, OperationInvariantEntry, PARITY_MATRIX, ParityMatrixEntry,
    ParityPolicy, SUPPORT_MATRIX, StateModel, SupportMatrixEntry, feature_spec, feature_specs,
    operation_invariant, operation_invariant_matrix, operation_parity, operation_spec,
    operation_specs, parity_matrix, support_matrix, version,
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
    if cfg!(feature = "cap-stereoisomers") {
        expected.push("cap-stereoisomers");
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
    assert_ne!(FunctionStatus::Native, FunctionStatus::Experimental);
    assert_ne!(
        ParityPolicy::RequiredWhenSupported,
        ParityPolicy::RequiredNow
    );
}

#[test]
fn metadata_type_rows_have_exact_public_roles_and_paths() {
    let expected = [
        (
            "types.FunctionStatus",
            "crate::FunctionStatus",
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
        assert_eq!(entry.status, FunctionStatus::Experimental);
        assert_eq!(entry.type_role, Some(role));
        assert_eq!(entry.callable, None);
    }

    let iterator = binding_entry("types.FeatureSpecIter");
    assert_eq!(iterator.item, BindingItem::Type);
    assert_eq!(iterator.owner, BindingOwner::Type);
    assert_eq!(compact(iterator.rust_path), "crate::FeatureSpecIter");
    assert_eq!(iterator.feature, "metadata");
    assert_eq!(iterator.status, FunctionStatus::Experimental);
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
            FunctionStatus::Experimental,
        ),
        (
            "module.feature_spec",
            "feature_spec",
            "featureSpec",
            "Option<&'staticcrate::FeatureSpec>",
            Some(("name", "&str")),
            FunctionStatus::Experimental,
        ),
        (
            "module.operation_specs",
            "operation_specs",
            "operationSpecs",
            "&'static[&'staticcrate::MoleculeOpSpec]",
            None,
            FunctionStatus::Experimental,
        ),
        (
            "module.operation_spec",
            "operation_spec",
            "operationSpec",
            "Option<&'staticcrate::MoleculeOpSpec>",
            Some(("method", "&str")),
            FunctionStatus::Experimental,
        ),
        (
            "module.support_matrix",
            "support_matrix",
            "supportMatrix",
            "&'static[crate::SupportMatrixEntry]",
            None,
            FunctionStatus::Experimental,
        ),
        (
            "module.operation_invariant_matrix",
            "operation_invariant_matrix",
            "operationInvariantMatrix",
            "&'static[crate::OperationInvariantEntry]",
            None,
            FunctionStatus::Experimental,
        ),
        (
            "module.parity_matrix",
            "parity_matrix",
            "parityMatrix",
            "&'static[crate::ParityMatrixEntry]",
            None,
            FunctionStatus::Experimental,
        ),
        (
            "module.operation_invariant",
            "operation_invariant",
            "operationInvariant",
            "Option<&'staticcrate::OperationInvariantEntry>",
            Some(("method", "&str")),
            FunctionStatus::Experimental,
        ),
        (
            "module.operation_parity",
            "operation_parity",
            "operationParity",
            "Option<&'staticcrate::ParityMatrixEntry>",
            Some(("method", "&str")),
            FunctionStatus::Experimental,
        ),
        (
            "module.version",
            "version",
            "version",
            "&'staticstr",
            None,
            FunctionStatus::Experimental,
        ),
    ];

    for (semantic_id, rust_name, javascript_name, output, parameter, status) in expected {
        let entry = binding_entry(semantic_id);
        let callable = entry.callable.expect("metadata row must be callable");
        assert_eq!(entry.item, BindingItem::Callable);
        assert_eq!(entry.owner, BindingOwner::Module);
        assert_eq!(compact(entry.rust_path), format!("crate::{rust_name}"));
        assert_eq!(entry.python_name, rust_name);
        assert_eq!(entry.javascript_name, javascript_name);
        assert_eq!(entry.feature, "metadata");
        assert_eq!(entry.status, status);
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
    assert_eq!(
        parity_matrix()
            .iter()
            .map(|row| row.operation.method)
            .collect::<Vec<_>>(),
        expected_parity_methods()
    );

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
        if let Some(parity_index) = expected_parity_methods()
            .iter()
            .position(|method| *method == operation.method)
        {
            let parity = operation_parity(operation.method).expect("original parity row");
            assert!(core::ptr::eq(parity, &parity_matrix()[parity_index]));
            assert!(core::ptr::eq(parity.operation, *operation));
        } else {
            assert_eq!(operation_parity(operation.method), None);
            assert_eq!(operation.parity, cosmolkit::ParityPolicy::NotApplicable);
            assert_eq!(operation.status, FunctionStatus::Native);
        }
    }

    if cfg!(feature = "cap-hydrogens") {
        let hydrogens = feature_spec("cap-hydrogens").expect("hydrogen feature");
        assert!(core::ptr::eq(
            hydrogens,
            feature_specs()
                .find(|feature| feature.name == "cap-hydrogens")
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
        assert_eq!(feature_spec("cap-hydrogens"), None);
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
