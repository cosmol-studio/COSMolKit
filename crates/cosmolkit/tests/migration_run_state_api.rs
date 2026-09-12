use std::fmt::Debug;

use cosmolkit::{
    BINDING_CONTRACT, BindingExposure, BindingItem, BindingOwner, BindingParity, BindingSupport,
    BindingTypeRole, MOLECULE_OPS, Molecule,
};

#[cfg(not(feature = "op-contracts-strict"))]
compile_error!("RUN-state public validation requires op-contracts-strict");

fn molecule_entry() -> &'static cosmolkit::BindingContractEntry {
    let entries = BINDING_CONTRACT
        .iter()
        .filter(|entry| entry.semantic_id == "types.Molecule")
        .collect::<Vec<_>>();
    assert_eq!(entries.len(), 1, "Molecule must have one canonical entry");
    entries[0]
}

fn assert_public_value_traits<T: Clone + Debug + Default + PartialEq>() {}

#[test]
fn canonical_molecule_type_is_the_real_public_entry() {
    assert_public_value_traits::<Molecule>();

    let empty = Molecule::default();
    let clone = empty.clone();
    assert_eq!(clone, empty);
    assert!(format!("{empty:?}").starts_with("Molecule"));
}

#[test]
fn binding_row_exactly_describes_the_run_state_type_surface() {
    let entry = molecule_entry();
    assert_eq!(entry.item, BindingItem::Type);
    assert_eq!(entry.owner, BindingOwner::Type);
    assert_eq!(entry.rust_path.replace(' ', ""), "crate::Molecule");
    assert_eq!(entry.python_name, "Molecule");
    assert_eq!(entry.javascript_name, "Molecule");
    assert_eq!(entry.feature, "runtime");
    assert_eq!(entry.exposure, BindingExposure::Public);
    assert_eq!(entry.support, BindingSupport::Supported);
    assert_eq!(entry.parity, BindingParity::NotApplicable);
    assert_eq!(entry.type_role, Some(BindingTypeRole::Value));
}

#[test]
fn lifecycle_type_has_no_callable_defaults_or_operation_contract() {
    let entry = molecule_entry();
    assert_eq!(entry.callable, None);
    assert!(
        MOLECULE_OPS
            .iter()
            .all(|operation| operation.method != entry.semantic_id),
        "a value type must not be projected as a molecule operation"
    );

    // RUN-state is only the live value container. Transformation mapping,
    // cache/stereo/property effects, failure atomicity, and multiple-output
    // validation belong to registered RUN operation units, so inventing an
    // empty operation contract here would create a second source of truth.
    assert!(entry.type_role.is_some());
}

#[test]
fn public_projection_keeps_state_cache_and_install_authority_private() {
    let root = include_str!("../src/lib.rs");
    let state = include_str!("../src/molecule.rs");

    assert_eq!(root.matches("pub use molecule::Molecule;").count(), 1);
    assert_eq!(state.matches("pub struct Molecule {").count(), 1);
    assert!(state.contains("struct MoleculeState {"));
    assert!(state.contains("state: Arc<MoleculeState>"));
    assert!(state.contains("pub(crate) struct DerivedCacheBlock {"));
    assert!(!state.contains("pub struct MoleculeState"));
    assert!(!state.contains("pub state:"));
    assert!(!state.contains("pub derived_cache:"));
    assert!(!state.contains("pub fn topology_mut"));
    assert!(!state.contains("pub fn coordinates_mut"));
    assert!(!state.contains("pub fn properties_mut"));
}
