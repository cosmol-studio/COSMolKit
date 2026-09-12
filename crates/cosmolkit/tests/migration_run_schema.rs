use std::collections::HashSet;

use cosmolkit::{
    BINDING_CONTRACT, BindingCallableContract, BindingDefault, BindingExposure, BindingItem,
    BindingKind, BindingOwner, BindingParameterContract, BindingParity, BindingSupport,
    BindingTypeRole, StateModel,
};

fn entry(semantic_id: &str) -> &'static cosmolkit::BindingContractEntry {
    BINDING_CONTRACT
        .iter()
        .find(|entry| entry.semantic_id == semantic_id)
        .unwrap_or_else(|| panic!("missing binding contract entry {semantic_id}"))
}

#[test]
fn canonical_registry_preserves_order_and_feature_local_subsets() {
    let mut expected = vec![
        "types.Molecule",
        "types.MoleculeBuilder",
        "Molecule.new",
        "Molecule.from_parts",
        "Molecule.to_builder",
        "Molecule.num_atoms",
        "Molecule.num_bonds",
        "Molecule.atoms",
        "Molecule.bonds",
        "Molecule.atom",
        "Molecule.bond",
        "Molecule.topology",
        "Molecule.coordinates",
        "Molecule.conformers",
        "Molecule.properties",
        "Molecule.property",
        "MoleculeBuilder.new",
        "MoleculeBuilder.from_parts",
        "MoleculeBuilder.build",
        "MoleculeBuilder.add_atom",
        "MoleculeBuilder.add_bond",
        "MoleculeBuilder.set_atom_formal_charge",
        "MoleculeBuilder.set_bond_order",
        "MoleculeBuilder.remove_bond_between_atoms",
        "MoleculeBuilder.degree",
        "MoleculeBuilder.neighbor_bonds",
        "MoleculeBuilder.bond_between_atoms",
        "MoleculeBuilder.atoms",
        "MoleculeBuilder.bonds",
        "MoleculeBuilder.substance_groups",
        "MoleculeBuilder.stereo_groups",
        "MoleculeBuilder.coordinates",
        "MoleculeBuilder.properties",
        "MoleculeBuilder.set_2d_coordinates",
        "MoleculeBuilder.add_2d_conformer",
        "MoleculeBuilder.add_3d_conformer",
        "MoleculeBuilder.add_substance_group",
        "MoleculeBuilder.add_stereo_group",
        "MoleculeBuilder.with_name",
        "MoleculeBuilder.with_property",
        "MoleculeBuilder.with_sdf_data_field",
        "MoleculeBuilder.with_properties",
        "types.OperationError",
        "types.SupportStatus",
        "types.ParityPolicy",
        "types.FeatureSpec",
        "types.MoleculeOpSpec",
        "types.SupportMatrixEntry",
        "types.OperationInvariantEntry",
        "types.ParityMatrixEntry",
        "types.FeatureSpecIter",
        "module.feature_specs",
        "module.feature_spec",
        "module.operation_specs",
        "module.operation_spec",
        "module.support_matrix",
        "module.operation_invariant_matrix",
        "module.parity_matrix",
        "module.operation_invariant",
        "module.operation_parity",
        "module.version",
    ];
    if cfg!(feature = "descriptors") {
        expected.extend([
            "Molecule.molecular_weight",
            "Molecule.exact_molecular_weight",
            "Molecule.molecular_formula",
        ]);
    }
    if cfg!(feature = "hydrogens") {
        expected.extend([
            "Molecule.with_hydrogens",
            "Molecule.add_hydrogens_",
            "Molecule.without_hydrogens",
            "Molecule.remove_hydrogens_",
        ]);
    }

    assert_eq!(
        BINDING_CONTRACT
            .iter()
            .map(|entry| entry.semantic_id)
            .collect::<Vec<_>>(),
        expected
    );
    assert_eq!(
        BINDING_CONTRACT
            .iter()
            .map(|entry| entry.semantic_id)
            .collect::<HashSet<_>>()
            .len(),
        BINDING_CONTRACT.len(),
        "semantic identities must be unique"
    );
}

#[test]
fn always_present_entries_have_exact_type_and_module_payloads() {
    let molecule = entry("types.Molecule");
    assert_eq!(molecule.item, BindingItem::Type);
    assert_eq!(molecule.owner, BindingOwner::Type);
    assert!(molecule.rust_path.ends_with("Molecule"));
    assert_eq!(molecule.python_name, "Molecule");
    assert_eq!(molecule.javascript_name, "Molecule");
    assert_eq!(molecule.feature, "runtime");
    assert_eq!(molecule.exposure, BindingExposure::Public);
    assert_eq!(molecule.support, BindingSupport::Supported);
    assert_eq!(molecule.parity, BindingParity::NotApplicable);
    assert_eq!(molecule.callable, None);
    assert_eq!(molecule.type_role, Some(BindingTypeRole::Value));

    let error = entry("types.OperationError");
    assert_eq!(error.item, BindingItem::Type);
    assert_eq!(error.owner, BindingOwner::Type);
    assert!(error.rust_path.ends_with("OperationError"));
    assert_eq!(error.type_role, Some(BindingTypeRole::Error));
    assert_eq!(error.callable, None);

    let version = entry("module.version");
    assert_eq!(version.item, BindingItem::Callable);
    assert_eq!(version.owner, BindingOwner::Module);
    assert_eq!(version.rust_path.replace(' ', ""), "crate::version");
    assert_eq!(version.python_name, "version");
    assert_eq!(version.javascript_name, "version");
    assert_eq!(version.feature, "metadata");
    assert_eq!(version.exposure, BindingExposure::Public);
    assert_eq!(version.support, BindingSupport::Supported);
    assert_eq!(version.parity, BindingParity::NotApplicable);
    assert_eq!(version.type_role, None);
    let callable = version.callable.expect("version callable metadata");
    assert_eq!(callable.kind, BindingKind::Module);
    assert!(callable.parameters.is_empty());
    assert_eq!(callable.output_type.replace(' ', ""), "&'staticstr");
    assert_eq!(callable.error_type, None);
    assert_eq!(callable.state_model, StateModel::ReadOnly);
    assert_eq!(callable.operation_semantic_id, None);
}

#[cfg(feature = "descriptors")]
#[test]
fn descriptor_entries_are_exact_canonical_read_only_methods() {
    for (semantic_id, rust_name, javascript, output) in [
        (
            "Molecule.molecular_weight",
            "molecular_weight",
            "molecularWeight",
            "f64",
        ),
        (
            "Molecule.exact_molecular_weight",
            "exact_molecular_weight",
            "exactMolecularWeight",
            "f64",
        ),
        (
            "Molecule.molecular_formula",
            "molecular_formula",
            "molecularFormula",
            "String",
        ),
    ] {
        let entry = entry(semantic_id);
        assert_eq!(entry.item, BindingItem::Callable);
        assert_eq!(entry.owner, BindingOwner::Molecule);
        assert!(entry.rust_path.replace(' ', "").ends_with(rust_name));
        assert_eq!(entry.python_name, rust_name);
        assert_eq!(entry.javascript_name, javascript);
        assert_eq!(entry.feature, "descriptors");
        assert_eq!(entry.exposure, BindingExposure::Public);
        assert_eq!(entry.support, BindingSupport::Experimental);
        assert_eq!(entry.parity, BindingParity::RequiredWhenSupported);
        assert_eq!(entry.type_role, None);
        let callable = entry.callable.expect("descriptor callable metadata");
        assert_eq!(callable.kind, BindingKind::Instance);
        assert!(callable.parameters.is_empty());
        assert_eq!(callable.output_type, output);
        assert_eq!(
            callable.error_type.map(|name| name.replace(' ', "")),
            Some("crate::OperationError".to_owned())
        );
        assert_eq!(callable.state_model, StateModel::ReadOnly);
        assert_eq!(callable.operation_semantic_id, None);
    }
}

#[cfg(feature = "hydrogens")]
#[test]
fn hydrogen_entries_bind_value_and_in_place_names_without_claiming_support() {
    for (semantic_id, rust_name, javascript, output, state) in [
        (
            "Molecule.with_hydrogens",
            "with_hydrogens",
            "withHydrogens",
            "crate::Molecule",
            StateModel::ValueReturning,
        ),
        (
            "Molecule.add_hydrogens_",
            "add_hydrogens_",
            "addHydrogens",
            "()",
            StateModel::InPlace,
        ),
        (
            "Molecule.without_hydrogens",
            "without_hydrogens",
            "withoutHydrogens",
            "crate::Molecule",
            StateModel::ValueReturning,
        ),
        (
            "Molecule.remove_hydrogens_",
            "remove_hydrogens_",
            "removeHydrogens",
            "()",
            StateModel::InPlace,
        ),
    ] {
        let entry = entry(semantic_id);
        assert_eq!(entry.owner, BindingOwner::Molecule);
        assert!(entry.rust_path.replace(' ', "").ends_with(rust_name));
        assert_eq!(entry.python_name, rust_name);
        assert_eq!(entry.javascript_name, javascript);
        assert_eq!(entry.feature, "hydrogens");
        assert_eq!(entry.exposure, BindingExposure::Public);
        assert_eq!(entry.support, BindingSupport::Unsupported);
        assert_eq!(entry.parity, BindingParity::RequiredWhenSupported);
        let callable = entry.callable.expect("hydrogen callable metadata");
        assert_eq!(callable.kind, BindingKind::Instance);
        assert!(callable.parameters.is_empty());
        assert_eq!(callable.output_type.replace(' ', ""), output);
        assert_eq!(
            callable.error_type.map(|name| name.replace(' ', "")),
            Some("crate::OperationError".to_owned())
        );
        assert_eq!(callable.state_model, state);
        assert_eq!(callable.operation_semantic_id, Some(rust_name));
    }
}

#[test]
fn schema_represents_every_closed_enum_and_parameter_payload_branch() {
    static PARAMETERS: [BindingParameterContract; 2] = [
        BindingParameterContract {
            name: "required_input",
            type_name: "&str",
            default: BindingDefault::Required,
        },
        BindingParameterContract {
            name: "configured_input",
            type_name: "bool",
            default: BindingDefault::Value("false"),
        },
    ];
    let callable = BindingCallableContract {
        kind: BindingKind::Static,
        parameters: &PARAMETERS,
        output_type: "u32",
        error_type: Some("OperationError"),
        state_model: StateModel::ValueReturning,
        operation_semantic_id: Some("from_example"),
    };
    assert_eq!(callable.kind, BindingKind::Static);
    assert_eq!(callable.parameters[0].default, BindingDefault::Required);
    assert_eq!(
        callable.parameters[1].default,
        BindingDefault::Value("false")
    );
    assert_eq!(callable.error_type, Some("OperationError"));
    assert_eq!(callable.state_model, StateModel::ValueReturning);

    assert_eq!(
        [
            BindingOwner::Molecule,
            BindingOwner::Module,
            BindingOwner::Type
        ]
        .len(),
        3
    );
    assert_eq!(
        [
            BindingKind::Instance,
            BindingKind::Static,
            BindingKind::Module
        ]
        .len(),
        3
    );
    assert_eq!(
        [BindingExposure::Registered, BindingExposure::Public].len(),
        2
    );
    assert_eq!(
        [
            BindingTypeRole::Value,
            BindingTypeRole::Parameter,
            BindingTypeRole::Result,
            BindingTypeRole::Error
        ]
        .len(),
        4
    );
    assert_eq!(
        [
            BindingSupport::Unsupported,
            BindingSupport::PreservedOnly,
            BindingSupport::Experimental,
            BindingSupport::Supported,
            BindingSupport::SupportedWithRdkitParity
        ]
        .len(),
        5
    );
    assert_eq!(
        [
            BindingParity::NotApplicable,
            BindingParity::RequiredWhenSupported,
            BindingParity::RequiredNow
        ]
        .len(),
        3
    );
    assert_eq!(
        [
            StateModel::ValueReturning,
            StateModel::InPlace,
            StateModel::ReadOnly
        ]
        .len(),
        3
    );
}

#[test]
fn registry_source_has_one_declaration_and_no_legacy_schema_or_cfg_products() {
    let registry = include_str!("../src/binding_contract/registry.rs");
    let types = include_str!("../src/binding_contract/types.rs");

    assert_eq!(registry.matches("binding_contract!").count(), 1);
    assert_eq!(registry.matches("pub static BINDING_CONTRACT").count(), 1);
    assert!(!registry.contains("cfg(all"));
    assert!(!registry.contains("cfg(not"));
    assert!(!registry.contains("ReturnKind"));
    assert!(!types.contains("ReturnKind"));
    assert!(!types.contains("type Binding"));
}
