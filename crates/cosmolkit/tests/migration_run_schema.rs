use cosmolkit::{
    BINDING_CONTRACT, BindingCallableContract, BindingDefault, BindingItem, BindingKind,
    BindingOwner, BindingParameterContract, BindingTypeRole, FunctionStatus, StateModel,
};

// Fixed public inventory of the owned force-field API, independent of generated metadata.
const PERSISTENT_FORCEFIELD_IDS: &[&str] = &[
    "types.ForceFieldDistanceConstraintParams",
    "ForceFieldDistanceConstraintParams.new",
    "MolecularForceField.add_distance_constraint_",
    "types.ForceFieldAngleConstraintParams",
    "ForceFieldAngleConstraintParams.new",
    "MolecularForceField.add_angle_constraint_",
    "types.ForceFieldTorsionConstraintParams",
    "ForceFieldTorsionConstraintParams.new",
    "MolecularForceField.add_torsion_constraint_",
    "types.ForceFieldPositionConstraintParams",
    "ForceFieldPositionConstraintParams.new",
    "MolecularForceField.add_position_constraint_",
    "types.MolecularForceFieldErrorKind",
    "types.MmffForceFieldParams",
    "MmffForceFieldParams.new",
    "types.UffForceFieldParams",
    "UffForceFieldParams.new",
    "types.ForceFieldMinimizeParams",
    "ForceFieldMinimizeParams.new",
    "types.ForceFieldEnergyGradient",
    "types.ForceFieldMinimizeOutcome",
    "types.ForceFieldError",
    "types.MmffForceFieldError",
    "types.UffForceFieldError",
    "types.MolecularForceField",
    "Molecule.mmff_force_field",
    "Molecule.mmff_force_field_with_params",
    "Molecule.uff_force_field",
    "Molecule.uff_force_field_with_params",
    "MolecularForceField.position",
    "MolecularForceField.set_position_",
    "MolecularForceField.positions",
    "MolecularForceField.set_positions_",
    "MolecularForceField.fixed_atoms",
    "MolecularForceField.set_fixed_atoms_",
    "MolecularForceField.energy",
    "MolecularForceField.gradient",
    "MolecularForceField.gradient_unconstrained",
    "MolecularForceField.energy_gradient",
    "MolecularForceField.minimize_",
    "MolecularForceField.minimize_with_params_",
];

fn entry(semantic_id: &str) -> &'static cosmolkit::BindingContractEntry {
    BINDING_CONTRACT
        .iter()
        .find(|entry| entry.semantic_id == semantic_id)
        .unwrap_or_else(|| panic!("missing binding contract entry {semantic_id}"))
}

fn canonical_svg_identity_status() -> FunctionStatus {
    FunctionStatus::ParityWithDifferences {
        reference: "RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8 MolDraw2DSVG",
        explanation: "ROOT-SVG-CANONICAL-METADATA-20261005: public SVG declares ck=https://kit.cosmol.org/ instead of the pinned source renderer identity; all other drawing bytes retain their source comparison.",
    }
}

#[cfg(feature = "cap-io")]
#[test]
fn sdf_reader_contracts_distinguish_concrete_and_preserving_results() {
    for (id, output) in [
        ("Molecule.from_sdf", "crate::Molecule"),
        ("Molecule.from_sdf_with_params", "crate::Molecule"),
        ("SdfRecord.from_sdf", "crate::SdfRecord"),
        ("SdfRecord.from_sdf_with_params", "crate::SdfRecord"),
    ] {
        let contract = entry(id);
        assert_eq!(contract.feature, "cap-io");
        assert_eq!(contract.status, FunctionStatus::Experimental);
        let callable = contract.callable.unwrap();
        assert_eq!(callable.output_type.replace(' ', ""), output);
        assert_eq!(
            callable.error_type.unwrap().replace(' ', ""),
            "crate::SdfError"
        );
        assert_eq!(callable.state_model, StateModel::ValueReturning);
        assert_eq!(callable.kind, BindingKind::Static);
    }
    for (id, role) in [
        ("types.SdfRecord", BindingTypeRole::Result),
        ("types.SdfGraph", BindingTypeRole::Result),
        (
            "types.SdfCoordinateMode",
            BindingTypeRole::ParameterSelector,
        ),
        ("types.SdfReadParams", BindingTypeRole::Parameter),
        ("types.SdfError", BindingTypeRole::Error),
    ] {
        let contract = entry(id);
        assert_eq!(contract.type_role, Some(role));
        assert_eq!(contract.feature, "cap-io");
        assert_eq!(contract.status, FunctionStatus::Experimental);
    }
    for (id, output, error) in [
        ("SdfRecord.graph", "&crate::SdfGraph", None),
        (
            "SdfRecord.molecule",
            "&crate::Molecule",
            Some("crate::SdfError"),
        ),
        (
            "SdfRecord.query_graph",
            "&crate::QueryGraph",
            Some("crate::SdfError"),
        ),
        (
            "SdfRecord.data_fields",
            "&[(crate::PropertyText,crate::PropertyText)]",
            None,
        ),
        ("SdfRecord.properties", "&crate::MoleculeProperties", None),
        (
            "SdfRecord.substance_groups",
            "&[crate::SubstanceGroup]",
            None,
        ),
        (
            "SdfRecord.source_coordinate_dim",
            "Option<crate::CoordinateDimension>",
            None,
        ),
    ] {
        let contract = entry(id);
        assert_eq!(contract.owner, BindingOwner::Type);
        assert_eq!(contract.status, FunctionStatus::Experimental);
        let callable = contract.callable.unwrap();
        assert_eq!(callable.kind, BindingKind::Instance);
        assert_eq!(callable.state_model, StateModel::ReadOnly);
        assert_eq!(callable.operation_semantic_id, None);
        assert_eq!(callable.output_type.replace(' ', ""), output);
        assert_eq!(
            callable.error_type.map(|name| name.replace(' ', "")),
            error.map(str::to_owned)
        );
    }
    for id in [
        "types.MolBlockReadParams",
        "types.MolBlockError",
        "Molecule.from_molblock",
        "Molecule.from_molblock_with_params",
    ] {
        assert!(
            BINDING_CONTRACT.iter().all(|entry| entry.semantic_id != id),
            "unimplemented interface must not appear in the public registry: {id}",
        );
    }
}

#[cfg(feature = "cap-forcefields")]
#[test]
fn uff_public_registry_has_exact_value_result_contract() {
    for (id, role) in [
        ("types.UffOptimizationParams", BindingTypeRole::Parameter),
        ("types.UffOptimizationResult", BindingTypeRole::Result),
        ("types.UffOptimizationError", BindingTypeRole::Error),
        ("types.UffOptimizationErrorKind", BindingTypeRole::Value),
    ] {
        let row = entry(id);
        assert_eq!(row.status, FunctionStatus::Experimental);
        assert_eq!(row.feature, "cap-forcefields");
        assert_eq!(row.type_role, Some(role));
    }
    for (method, javascript, parameter_count) in [
        ("with_uff_optimized", "withUffOptimized", 0),
        (
            "with_uff_optimized_with_params",
            "withUffOptimizedWithParams",
            1,
        ),
    ] {
        let row = entry(&format!("Molecule.{method}"));
        assert_eq!(row.status, FunctionStatus::Experimental);
        assert_eq!(row.python_name, method);
        assert_eq!(row.javascript_name, javascript);
        let callable = row.callable.unwrap();
        assert_eq!(callable.kind, BindingKind::Instance);
        assert_eq!(callable.state_model, StateModel::ValueReturning);
        assert_eq!(callable.parameters.len(), parameter_count);
        assert_eq!(
            callable.output_type.replace(' ', ""),
            "crate::UffOptimizationResult"
        );
        assert_eq!(
            callable.error_type.map(|name| name.replace(' ', "")),
            Some("crate::OperationError".into())
        );
        assert_eq!(callable.operation_semantic_id, Some(method));
    }
    let _: fn(
        &cosmolkit::Molecule,
    ) -> Result<cosmolkit::UffOptimizationResult, cosmolkit::OperationError> =
        cosmolkit::Molecule::with_uff_optimized;
    let _: fn(
        &cosmolkit::Molecule,
        &cosmolkit::UffOptimizationParams,
    ) -> Result<cosmolkit::UffOptimizationResult, cosmolkit::OperationError> =
        cosmolkit::Molecule::with_uff_optimized_with_params;
    let _: fn(&cosmolkit::UffOptimizationError) -> cosmolkit::UffOptimizationErrorKind =
        cosmolkit::UffOptimizationError::kind;
    let params = cosmolkit::UffOptimizationParams::default();
    assert_eq!(params.max_iterations, 1000);
    assert_eq!(params.vdw_threshold.to_bits(), 10.0_f64.to_bits());
    assert!(params.ignore_interfragment_interactions);
    assert_eq!(params.conformer_id, None);
}

#[cfg(feature = "cap-forcefields")]
#[test]
fn uff_param_p10_registry_has_exact_read_only_query_contract() {
    let forcefields_ids = BINDING_CONTRACT
        .iter()
        .filter(|row| row.feature == "cap-forcefields")
        .map(|row| row.semantic_id)
        .collect::<Vec<_>>();
    assert_eq!(
        forcefields_ids,
        vec![
            "types.UffEvaluationParams",
            "types.UffEnergyGradient",
            "Molecule.uff_energy_gradient",
            "Molecule.uff_energy_gradient_with_params",
            "UffEnergyGradient.energy",
            "UffEnergyGradient.gradient",
            "UffOptimizationResult.molecule",
            "UffOptimizationResult.status_code",
            "UffOptimizationResult.needs_more",
            "UffOptimizationResult.energy",
            "UffConformerResult.conformer_id",
            "UffConformerResult.status_code",
            "UffConformerResult.needs_more",
            "UffConformerResult.energy",
            "UffConformerOptimizationResult.molecule",
            "UffConformerOptimizationResult.conformer_results",
            "MmffAtomProperties.atom_type",
            "MmffAtomProperties.formal_charge",
            "MmffAtomProperties.partial_charge",
            "types.MmffEvaluationParams",
            "types.MmffEnergyGradient",
            "Molecule.mmff_energy_gradient",
            "Molecule.mmff_energy_gradient_with_params",
            "MmffEnergyGradient.energy",
            "MmffEnergyGradient.gradient",
            "MmffOptimizeMoleculeResult.molecule",
            "MmffOptimizeMoleculeResult.needs_more",
            "MmffOptimizeMoleculeResult.status_code",
            "MmffOptimizeMoleculeConfsResult.molecule",
            "MmffOptimizeMoleculeConfsResult.conformer_results",
            "MmffOptimizeMoleculeConfResult.needs_more",
            "MmffOptimizeMoleculeConfResult.status_code",
            "MmffOptimizeMoleculeConfResult.energy",
            "MmffProperties.is_valid",
            "MmffProperties.variant",
            "MmffProperties.atoms",
            "MmffProperties.atom_type",
            "MmffProperties.formal_charge",
            "MmffProperties.partial_charge",
            "types.MmffOptimizationParams",
            "types.MmffConformerOptimizationParams",
            "types.MmffOptimizeMoleculeResult",
            "types.MmffOptimizeMoleculeConfResult",
            "types.MmffOptimizeMoleculeConfsResult",
            "types.MmffOptimizationError",
            "Molecule.with_mmff_optimized",
            "Molecule.with_mmff_optimized_with_params",
            "Molecule.with_mmff_optimized_conformers",
            "Molecule.with_mmff_optimized_conformers_with_params",
            "types.UffParameterQueryError",
            "types.UffParameterError",
            "types.UffParameterErrorKind",
            "UffParameterError.kind",
            "types.MmffProperties",
            "types.MmffPropertiesParams",
            "types.MmffAtomProperties",
            "types.MmffVariant",
            "types.MmffMolPropertiesError",
            "Molecule.mmff_has_all_molecule_params",
            "Molecule.mmff_properties",
            "Molecule.mmff_properties_with_params",
            "Molecule.uff_has_all_molecule_params",
            "types.UffOptimizationParams",
            "types.UffOptimizationResult",
            "types.UffOptimizationError",
            "types.UffOptimizationErrorKind",
            "types.UffConformerOptimizationParams",
            "types.UffConformerOptimizationResult",
            "types.UffConformerResult",
            "UffOptimizationError.kind",
            "Molecule.with_uff_optimized",
            "Molecule.with_uff_optimized_with_params",
            "Molecule.with_uff_optimized_conformers",
            "Molecule.with_uff_optimized_conformers_with_params",
        ]
        .into_iter()
        .chain(PERSISTENT_FORCEFIELD_IDS.iter().copied())
        .collect::<Vec<_>>(),
        "complete forcefields exposes the exact UFF and MMFF identities"
    );

    for (semantic_id, rust_name, python_name, javascript_name, role) in [
        (
            "types.UffParameterQueryError",
            "UffParameterQueryError",
            "UffParameterQueryError",
            "UffParameterQueryError",
            BindingTypeRole::Error,
        ),
        (
            "types.UffParameterError",
            "UffParameterError",
            "UffParameterError",
            "UffParameterError",
            BindingTypeRole::Error,
        ),
        (
            "types.UffParameterErrorKind",
            "UffParameterErrorKind",
            "UffParameterErrorKind",
            "UffParameterErrorKind",
            BindingTypeRole::Value,
        ),
    ] {
        let row = entry(semantic_id);
        assert_eq!(row.item, BindingItem::Type);
        assert_eq!(row.owner, BindingOwner::Type);
        assert!(row.rust_path.ends_with(rust_name));
        assert_eq!(row.python_name, python_name);
        assert_eq!(row.javascript_name, javascript_name);
        assert_eq!(row.feature, "cap-forcefields");
        assert_eq!(row.status, FunctionStatus::Experimental);
        assert_eq!(row.callable, None);
        assert_eq!(row.type_role, Some(role));
    }

    let kind = entry("UffParameterError.kind");
    assert_eq!(kind.item, BindingItem::Callable);
    assert_eq!(kind.owner, BindingOwner::Type);
    assert_eq!(
        kind.rust_path.replace(' ', ""),
        "crate::UffParameterError::kind"
    );
    assert_eq!(kind.python_name, "kind");
    assert_eq!(kind.javascript_name, "kind");
    assert_eq!(kind.feature, "cap-forcefields");
    assert_eq!(kind.status, FunctionStatus::Experimental);
    assert_eq!(kind.type_role, None);
    let kind_callable = kind.callable.expect("kind accessor callable metadata");
    assert_eq!(kind_callable.kind, BindingKind::Instance);
    assert!(kind_callable.parameters.is_empty());
    assert_eq!(
        kind_callable.output_type.replace(' ', ""),
        "crate::UffParameterErrorKind"
    );
    assert_eq!(kind_callable.error_type, None);
    assert_eq!(kind_callable.state_model, StateModel::ReadOnly);
    assert_eq!(kind_callable.operation_semantic_id, None);
    let _kind_signature: fn(&cosmolkit::UffParameterError) -> cosmolkit::UffParameterErrorKind =
        cosmolkit::UffParameterError::kind;
    let kind_name = |value| match value {
        cosmolkit::UffParameterErrorKind::Preparation => "Preparation",
        cosmolkit::UffParameterErrorKind::ParameterTable => "ParameterTable",
        cosmolkit::UffParameterErrorKind::Typing => "Typing",
    };
    assert_eq!(
        [
            kind_name(cosmolkit::UffParameterErrorKind::Preparation),
            kind_name(cosmolkit::UffParameterErrorKind::ParameterTable),
            kind_name(cosmolkit::UffParameterErrorKind::Typing),
        ],
        ["Preparation", "ParameterTable", "Typing"]
    );

    let query = entry("Molecule.uff_has_all_molecule_params");
    assert_eq!(query.item, BindingItem::Callable);
    assert_eq!(query.owner, BindingOwner::Molecule);
    assert_eq!(
        query.rust_path.replace(' ', ""),
        "crate::Molecule::uff_has_all_molecule_params"
    );
    assert_eq!(query.python_name, "uff_has_all_molecule_params");
    assert_eq!(query.javascript_name, "uffHasAllMoleculeParams");
    assert_eq!(query.feature, "cap-forcefields");
    assert_eq!(query.status, FunctionStatus::Experimental);
    assert_eq!(query.type_role, None);
    let query_callable = query.callable.expect("query callable metadata");
    assert_eq!(query_callable.kind, BindingKind::Instance);
    assert!(query_callable.parameters.is_empty());
    assert_eq!(query_callable.output_type, "bool");
    assert_eq!(
        query_callable.error_type.map(|name| name.replace(' ', "")),
        Some("crate::UffParameterQueryError".to_owned())
    );
    assert_eq!(query_callable.state_model, StateModel::ReadOnly);
    assert_eq!(query_callable.operation_semantic_id, None);
    let _query_signature: fn(
        &cosmolkit::Molecule,
    ) -> Result<bool, cosmolkit::UffParameterQueryError> =
        cosmolkit::Molecule::uff_has_all_molecule_params;
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
    assert_eq!(molecule.status, FunctionStatus::Experimental);
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
    assert_eq!(version.status, FunctionStatus::Experimental);
    assert_eq!(version.type_role, None);
    let callable = version.callable.expect("version callable metadata");
    assert_eq!(callable.kind, BindingKind::Module);
    assert!(callable.parameters.is_empty());
    assert_eq!(callable.output_type.replace(' ', ""), "&'staticstr");
    assert_eq!(callable.error_type, None);
    assert_eq!(callable.state_model, StateModel::ReadOnly);
    assert_eq!(callable.operation_semantic_id, None);
}

#[cfg(feature = "cap-depict")]
#[test]
fn drawing_entries_have_exact_experimental_shared_query_contracts() {
    let error = entry("types.DrawingError");
    assert_eq!(error.item, BindingItem::Type);
    assert_eq!(error.owner, BindingOwner::Type);
    assert_eq!(error.rust_path.replace(' ', ""), "crate::DrawingError");
    assert_eq!(error.python_name, "DrawingError");
    assert_eq!(error.javascript_name, "DrawingError");
    assert_eq!(error.feature, "cap-depict");
    assert_eq!(error.status, FunctionStatus::Experimental);
    assert_eq!(error.type_role, Some(BindingTypeRole::Error));
    assert_eq!(error.callable, None);
    for (name, javascript, output) in [
        ("to_svg", "toSvg", "String"),
        ("to_png", "toPng", "Vec<u8>"),
    ] {
        let row = entry(&format!("Molecule.{name}"));
        assert_eq!(row.item, BindingItem::Callable);
        assert_eq!(row.owner, BindingOwner::Molecule);
        assert_eq!(
            row.rust_path.replace(' ', ""),
            format!("crate::Molecule::{name}")
        );
        assert_eq!(row.python_name, name);
        assert_eq!(row.javascript_name, javascript);
        assert_eq!(row.feature, "cap-depict");
        assert_eq!(
            row.status,
            if name == "to_svg" {
                canonical_svg_identity_status()
            } else {
                FunctionStatus::Experimental
            }
        );
        let callable = row.callable.unwrap();
        assert_eq!(callable.kind, BindingKind::Instance);
        assert_eq!(callable.receiver, Some(cosmolkit::BindingReceiver::Shared));
        assert_eq!(callable.state_model, StateModel::ReadOnly);
        assert_eq!(callable.operation_semantic_id, None);
        assert_eq!(callable.output_type.replace(' ', ""), output);
        assert_eq!(
            callable.error_type.unwrap().replace(' ', ""),
            "crate::DrawingError"
        );
        assert_eq!(callable.parameters.len(), 2);
        for (parameter, expected_name) in callable.parameters.iter().zip(["width", "height"]) {
            assert_eq!(parameter.name, expected_name);
            assert_eq!(parameter.type_name, "u32");
            assert_eq!(parameter.default, BindingDefault::Required);
        }
    }
}

#[cfg(feature = "cap-descriptors")]
#[test]
fn descriptor_entries_are_exact_canonical_read_only_methods() {
    for (semantic_id, rust_name, javascript, output, error) in [
        (
            "Molecule.molecular_weight",
            "molecular_weight",
            "molecularWeight",
            "f64",
            "crate::DescriptorReadError",
        ),
        (
            "Molecule.exact_molecular_weight",
            "exact_molecular_weight",
            "exactMolecularWeight",
            "f64",
            "crate::DescriptorReadError",
        ),
        (
            "Molecule.molecular_formula",
            "molecular_formula",
            "molecularFormula",
            "String",
            "crate::DescriptorReadError",
        ),
        (
            "Molecule.num_heavy_atoms",
            "num_heavy_atoms",
            "numHeavyAtoms",
            "u32",
            "crate::DescriptorReadError",
        ),
        (
            "Molecule.total_atom_count",
            "total_atom_count",
            "totalAtomCount",
            "u32",
            "crate::DescriptorReadError",
        ),
        (
            "Molecule.lipinski_hba",
            "lipinski_hba",
            "lipinskiHba",
            "u32",
            "crate::DescriptorReadError",
        ),
        (
            "Molecule.lipinski_hbd",
            "lipinski_hbd",
            "lipinskiHbd",
            "u32",
            "crate::DescriptorReadError",
        ),
        (
            "Molecule.fraction_csp3",
            "fraction_csp3",
            "fractionCsp3",
            "f64",
            "crate::DescriptorReadError",
        ),
    ] {
        let entry = entry(semantic_id);
        assert_eq!(entry.item, BindingItem::Callable);
        assert_eq!(entry.owner, BindingOwner::Molecule);
        assert!(entry.rust_path.replace(' ', "").ends_with(rust_name));
        assert_eq!(entry.python_name, rust_name);
        assert_eq!(entry.javascript_name, javascript);
        assert_eq!(entry.feature, "cap-descriptors");
        assert_eq!(entry.status, FunctionStatus::Experimental);
        assert_eq!(entry.type_role, None);
        let callable = entry.callable.expect("descriptor callable metadata");
        assert_eq!(callable.kind, BindingKind::Instance);
        assert!(callable.parameters.is_empty());
        assert_eq!(callable.output_type, output);
        assert_eq!(
            callable.error_type.map(|name| name.replace(' ', "")),
            Some(error.replace(' ', ""))
        );
        assert_eq!(callable.state_model, StateModel::ReadOnly);
        assert_eq!(callable.operation_semantic_id, None);
    }
}

#[cfg(feature = "cap-descriptors")]
#[test]
fn descriptor_error_type_entries_carry_error_roles_and_exact_names() {
    for (semantic_id, rust_type, projection) in [
        (
            "types.DescriptorError",
            "DescriptorError",
            "DescriptorError",
        ),
        (
            "types.DescriptorReadError",
            "DescriptorReadError",
            "DescriptorReadError",
        ),
    ] {
        let entry = entry(semantic_id);
        assert_eq!(entry.item, BindingItem::Type);
        assert_eq!(entry.owner, BindingOwner::Type);
        assert!(entry.rust_path.replace(' ', "").ends_with(rust_type));
        assert_eq!(entry.python_name, projection);
        assert_eq!(entry.javascript_name, projection);
        assert_eq!(entry.feature, "cap-descriptors");
        assert_eq!(entry.status, FunctionStatus::Experimental);
        assert_eq!(entry.type_role, Some(BindingTypeRole::Error));
        assert!(entry.callable.is_none());
    }
}

#[cfg(feature = "cap-hydrogens")]
#[test]
fn hydrogen_entries_bind_value_and_in_place_names_with_parity_support() {
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
            "addHydrogens_",
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
            "removeHydrogens_",
            "()",
            StateModel::InPlace,
        ),
    ] {
        let entry = entry(semantic_id);
        assert_eq!(entry.owner, BindingOwner::Molecule);
        assert!(entry.rust_path.replace(' ', "").ends_with(rust_name));
        assert_eq!(entry.python_name, rust_name);
        assert_eq!(entry.javascript_name, javascript);
        assert_eq!(entry.feature, "cap-hydrogens");
        assert_eq!(entry.status, FunctionStatus::Experimental);
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
        receiver: None,
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
            FunctionStatus::Parity { reference: "RDKit" },
            FunctionStatus::ParityWithDifferences {
                reference: "Gemmi",
                explanation: "Approved difference for the documented boundary"
            },
            FunctionStatus::Native,
            FunctionStatus::Experimental
        ]
        .len(),
        4
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
    assert_eq!(registry.matches("pub static BINDING_CONTRACT =").count(), 1);
    assert!(!registry.contains("BINDING_CONTRACT_PYTHON_ALIASES"));
    assert!(registry.contains("convenience: bio_from_mmcif"));
    // The approved wasm32 platform exclusion is the only compound gate:
    // native archive APIs must not be emitted even when `full` is selected.
    let native_archive =
        "#[cfg(all(feature = \"cap-serialization\", not(target_arch = \"wasm32\")))]";
    // Conveniences inherit the gate from these three native declarations.
    assert_eq!(registry.matches(native_archive).count(), 3);
    let remaining = registry.replace(native_archive, "");
    assert!(!remaining.contains("cfg(all"));
    assert!(!remaining.contains("cfg(not"));
    assert!(!registry.contains("ReturnKind"));
    assert!(!types.contains("ReturnKind"));
    assert!(!types.contains("type Binding"));
    for old_field in ["exposure:", "support:", "parity:"] {
        assert!(!registry.contains(old_field));
    }
}
