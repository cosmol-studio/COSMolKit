use cosmolkit::{
    BINDING_CONTRACT, BindingDefault, BindingItem, BindingKind, BindingOwner, BindingTypeRole,
    FunctionStatus, StateModel,
};

// Fixed public inventory of the owned force-field API, independent of generated metadata.
const PERSISTENT_FORCEFIELD_IDS: &[&str] = &[
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

fn compact(text: &str) -> String {
    text.chars()
        .filter(|character| !character.is_whitespace())
        .collect()
}

#[cfg(feature = "cap-forcefields")]
#[test]
fn uff_conformer_binding_schema_has_the_exact_order_and_compiled_shapes() {
    let forcefields_ids = BINDING_CONTRACT
        .iter()
        .filter(|row| row.feature == "cap-forcefields")
        .map(|row| row.semantic_id)
        .collect::<Vec<_>>();
    assert_eq!(
        forcefields_ids,
        [
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
        .collect::<Vec<_>>()
    );

    for (semantic_id, rust_name, role) in [
        (
            "types.UffConformerOptimizationParams",
            "UffConformerOptimizationParams",
            BindingTypeRole::Parameter,
        ),
        (
            "types.UffConformerOptimizationResult",
            "UffConformerOptimizationResult",
            BindingTypeRole::Result,
        ),
        (
            "types.UffConformerResult",
            "UffConformerResult",
            BindingTypeRole::Value,
        ),
    ] {
        let row = entry(semantic_id);
        assert_eq!(row.item, BindingItem::Type);
        assert_eq!(row.owner, BindingOwner::Type);
        assert_eq!(compact(row.rust_path), format!("crate::{rust_name}"));
        assert_eq!(row.python_name, rust_name);
        assert_eq!(row.javascript_name, rust_name);
        assert_eq!(row.feature, "cap-forcefields");
        assert_eq!(row.status, FunctionStatus::Experimental);
        assert_eq!(row.callable, None);
        assert_eq!(row.type_role, Some(role));
    }

    for (semantic_id, rust_name, python_name, javascript_name, parameter_count) in [
        (
            "Molecule.with_uff_optimized_conformers",
            "with_uff_optimized_conformers",
            "with_uff_optimized_conformers",
            "withUffOptimizedConformers",
            0,
        ),
        (
            "Molecule.with_uff_optimized_conformers_with_params",
            "with_uff_optimized_conformers_with_params",
            "with_uff_optimized_conformers_with_params",
            "withUffOptimizedConformersWithParams",
            1,
        ),
    ] {
        let row = entry(semantic_id);
        assert_eq!(row.item, BindingItem::Callable);
        assert_eq!(row.owner, BindingOwner::Molecule);
        assert_eq!(
            compact(row.rust_path),
            format!("crate::Molecule::{rust_name}")
        );
        assert_eq!(row.python_name, python_name);
        assert_eq!(row.javascript_name, javascript_name);
        assert_eq!(row.feature, "cap-forcefields");
        assert_eq!(row.status, FunctionStatus::Experimental);
        assert_eq!(row.type_role, None);
        let callable = row.callable.expect("operation callable metadata");
        assert_eq!(callable.kind, BindingKind::Instance);
        assert_eq!(callable.state_model, StateModel::ValueReturning);
        assert_eq!(callable.parameters.len(), parameter_count);
        assert_eq!(
            compact(callable.output_type),
            "crate::UffConformerOptimizationResult"
        );
        assert_eq!(
            callable.error_type.map(|name| compact(name)),
            Some("crate::OperationError".to_owned())
        );
        assert_eq!(callable.operation_semantic_id, Some(rust_name));
        if parameter_count == 1 {
            let parameter = callable.parameters[0];
            assert_eq!(parameter.name, "params");
            assert_eq!(
                compact(parameter.type_name),
                "&crate::UffConformerOptimizationParams"
            );
            assert_eq!(parameter.default, BindingDefault::Required);
        }
    }

    let _: fn(
        &cosmolkit::Molecule,
    ) -> Result<cosmolkit::UffConformerOptimizationResult, cosmolkit::OperationError> =
        cosmolkit::Molecule::with_uff_optimized_conformers;
    let _: fn(
        &cosmolkit::Molecule,
        &cosmolkit::UffConformerOptimizationParams,
    ) -> Result<cosmolkit::UffConformerOptimizationResult, cosmolkit::OperationError> =
        cosmolkit::Molecule::with_uff_optimized_conformers_with_params;

    let kind_name = |kind| match kind {
        cosmolkit::UffOptimizationErrorKind::MissingConformer { .. } => "MissingConformer",
        cosmolkit::UffOptimizationErrorKind::Rings => "Rings",
        cosmolkit::UffOptimizationErrorKind::Optimization => "Optimization",
        cosmolkit::UffOptimizationErrorKind::ConformerOptimization => "ConformerOptimization",
        cosmolkit::UffOptimizationErrorKind::Evaluation => "Evaluation",
    };
    assert_eq!(
        kind_name(cosmolkit::UffOptimizationErrorKind::ConformerOptimization),
        "ConformerOptimization"
    );

    assert_eq!(
        kind_name(cosmolkit::UffOptimizationErrorKind::Evaluation),
        "Evaluation"
    );

    let params = cosmolkit::UffConformerOptimizationParams::default();
    assert_eq!(params.max_iterations, 1000);
    assert_eq!(params.vdw_threshold.to_bits(), 10.0_f64.to_bits());
    assert!(params.ignore_interfragment_interactions);
}
