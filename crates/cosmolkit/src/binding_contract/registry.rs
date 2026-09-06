//! Canonical cross-language naming registry.

#[cfg(any(feature = "hydrogens", feature = "descriptors"))]
use cosmolkit_macros::binding_contract;

#[cfg(all(feature = "hydrogens", feature = "descriptors"))]
binding_contract! {
    pub static BINDING_CONTRACT = [
        {
            semantic_id: "Molecule.molecular_weight", owner: molecule, kind: instance,
            rust: molecular_weight, python: "molecular_weight", javascript: "molecularWeight",
            feature: "descriptors", returns: result_f64, state: read_only, signature: read_only_f64,
        },
        {
            semantic_id: "Molecule.exact_molecular_weight", owner: molecule, kind: instance,
            rust: exact_molecular_weight, python: "exact_molecular_weight", javascript: "exactMolecularWeight",
            feature: "descriptors", returns: result_f64, state: read_only, signature: read_only_f64,
        },
        {
            semantic_id: "Molecule.molecular_formula", owner: molecule, kind: instance,
            rust: molecular_formula, python: "molecular_formula", javascript: "molecularFormula",
            feature: "descriptors", returns: result_string, state: read_only, signature: read_only_string,
        },
        {
            semantic_id: "Molecule.without_hydrogens", owner: molecule, kind: instance,
            rust: without_hydrogens, python: "without_hydrogens", javascript: "withoutHydrogens",
            feature: "hydrogens", returns: result_molecule, state: value_returning, signature: value_molecule,
        },
        {
            semantic_id: "Molecule.remove_hydrogens_", owner: molecule, kind: instance,
            rust: remove_hydrogens_, python: "remove_hydrogens_", javascript: "removeHydrogens",
            feature: "hydrogens", returns: result_unit, state: in_place, signature: in_place_unit,
        },
    ];
}

#[cfg(all(feature = "descriptors", not(feature = "hydrogens")))]
binding_contract! {
    pub static BINDING_CONTRACT = [
        {
            semantic_id: "Molecule.molecular_weight", owner: molecule, kind: instance,
            rust: molecular_weight, python: "molecular_weight", javascript: "molecularWeight",
            feature: "descriptors", returns: result_f64, state: read_only, signature: read_only_f64,
        },
        {
            semantic_id: "Molecule.exact_molecular_weight", owner: molecule, kind: instance,
            rust: exact_molecular_weight, python: "exact_molecular_weight", javascript: "exactMolecularWeight",
            feature: "descriptors", returns: result_f64, state: read_only, signature: read_only_f64,
        },
        {
            semantic_id: "Molecule.molecular_formula", owner: molecule, kind: instance,
            rust: molecular_formula, python: "molecular_formula", javascript: "molecularFormula",
            feature: "descriptors", returns: result_string, state: read_only, signature: read_only_string,
        },
    ];
}

#[cfg(all(feature = "hydrogens", not(feature = "descriptors")))]
binding_contract! {
    pub static BINDING_CONTRACT = [
        {
            semantic_id: "Molecule.without_hydrogens", owner: molecule, kind: instance,
            rust: without_hydrogens, python: "without_hydrogens", javascript: "withoutHydrogens",
            feature: "hydrogens", returns: result_molecule, state: value_returning, signature: value_molecule,
        },
        {
            semantic_id: "Molecule.remove_hydrogens_", owner: molecule, kind: instance,
            rust: remove_hydrogens_, python: "remove_hydrogens_", javascript: "removeHydrogens",
            feature: "hydrogens", returns: result_unit, state: in_place, signature: in_place_unit,
        },
    ];
}

#[cfg(not(any(feature = "hydrogens", feature = "descriptors")))]
pub static BINDING_CONTRACT: &[super::BindingContractEntry] = &[];

#[cfg(all(test, feature = "descriptors"))]
mod descriptor_tests {
    use super::BINDING_CONTRACT;
    use crate::{BindingKind, BindingOwner, ReturnKind, StateModel};

    #[test]
    fn descriptor_names_are_canonical_molecule_methods() {
        let expected = [
            (
                "Molecule.molecular_weight",
                "molecular_weight",
                "molecularWeight",
            ),
            (
                "Molecule.exact_molecular_weight",
                "exact_molecular_weight",
                "exactMolecularWeight",
            ),
            (
                "Molecule.molecular_formula",
                "molecular_formula",
                "molecularFormula",
            ),
        ];
        for (semantic_id, rust_name, javascript_name) in expected {
            let entry = BINDING_CONTRACT
                .iter()
                .find(|entry| entry.semantic_id == semantic_id)
                .expect("descriptor registry entry");
            assert_eq!(entry.owner, BindingOwner::Molecule);
            assert_eq!(entry.kind, BindingKind::Instance);
            assert_eq!(entry.rust_name, rust_name);
            assert_eq!(entry.javascript_name, javascript_name);
            assert_eq!(
                entry.return_kind,
                if semantic_id.ends_with("formula") {
                    ReturnKind::ResultString
                } else {
                    ReturnKind::ResultF64
                }
            );
            assert_eq!(entry.state_model, StateModel::ReadOnly);
        }
    }
}

#[cfg(all(test, feature = "hydrogens"))]
mod hydrogen_tests {
    use super::BINDING_CONTRACT;
    use crate::{Molecule, OperationError};

    #[test]
    fn hydrogen_entries_remain_bound_to_the_runtime_surface() {
        assert!(
            BINDING_CONTRACT
                .iter()
                .any(|entry| entry.semantic_id == "Molecule.without_hydrogens")
        );
        let error = Molecule::new()
            .without_hydrogens()
            .expect_err("hydrogen migration skeleton remains unsupported");
        assert_eq!(
            error,
            OperationError::Unsupported {
                operation: "remove_hydrogens"
            }
        );
    }
}
