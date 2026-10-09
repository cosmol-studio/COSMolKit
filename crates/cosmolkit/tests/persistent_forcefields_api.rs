//! Public persistent force-field API regressions: defaults, ownership and error contracts.
#![cfg(all(
    feature = "cap-forcefields",
    feature = "cap-smiles",
    feature = "cap-valence"
))]
use cosmolkit::{
    AtomId, Conformer3D, CoordinateBlock, ForceFieldMinimizeParams, MmffForceFieldParams,
    MolecularForceFieldErrorKind, Molecule, UffForceFieldParams,
};
use std::error::Error;

fn fixture() -> Molecule {
    let template = Molecule::from_smiles("CC.CC").unwrap();
    Molecule::from_parts(
        template.topology().clone(),
        CoordinateBlock {
            conformers_3d: vec![
                Conformer3D::new(
                    7,
                    vec![[0., 0., 0.], [1.9, 0.2, 0.], [5., 1., 0.], [6.7, 1.1, 0.3]],
                    true,
                ),
                Conformer3D::new(
                    42,
                    vec![[0., 0., 0.], [2.1, 0.3, 0.], [5., 2., 0.], [7., 2.1, 0.4]],
                    true,
                ),
            ],
            ..Default::default()
        },
        template.properties().clone(),
    )
    .unwrap()
    .with_assigned_valence()
    .unwrap()
}
#[test]
fn immutable_defaults_and_registered_properties() {
    let mmff = MmffForceFieldParams::default();
    assert_eq!(mmff.conformer_id(), None);
    assert_eq!(mmff.mmff_variant(), "MMFF94");
    assert_eq!(mmff.non_bonded_threshold(), 100.);
    assert!(mmff.ignore_interfragment_interactions());
    let uff = UffForceFieldParams::default();
    assert_eq!(uff.conformer_id(), None);
    assert_eq!(uff.vdw_threshold(), 10.);
    assert!(uff.ignore_interfragment_interactions());
    let minimize = ForceFieldMinimizeParams::default();
    assert_eq!(minimize.max_iterations(), 200);
    assert_eq!(minimize.force_tolerance(), 1e-4);
    assert_eq!(minimize.energy_tolerance(), 1e-6);
    for (name, count) in [
        ("MmffForceFieldParams", 4),
        ("UffForceFieldParams", 3),
        ("ForceFieldMinimizeParams", 3),
        ("ForceFieldEnergyGradient", 2),
        ("ForceFieldMinimizeOutcome", 3),
    ] {
        assert_eq!(
            cosmolkit::BINDING_CONTRACT_PROPERTIES
                .iter()
                .filter(|p| p.type_semantic_id == format!("types.{name}"))
                .count(),
            count
        );
    }
    assert_eq!(
        cosmolkit::BINDING_CONTRACT_KEYWORDS
            .iter()
            .filter(|p| matches!(
                p.semantic_id,
                "Molecule.mmff_force_field"
                    | "Molecule.uff_force_field"
                    | "MolecularForceField.minimize_"
            ))
            .count(),
        3
    );
}
#[test]
fn four_factories_select_existing_owned_coordinates() {
    let mol = fixture();
    let original = mol.conformers_3d().to_vec();
    let clone = mol.clone();
    let mut a = mol.mmff_force_field().unwrap();
    let b = mol
        .mmff_force_field_with_params(&MmffForceFieldParams::default())
        .unwrap();
    let mut c = mol.uff_force_field().unwrap();
    let d = mol
        .uff_force_field_with_params(&UffForceFieldParams::new(Some(42), 10., true))
        .unwrap();
    assert_eq!(a.positions(), original[0].coordinates());
    assert_eq!(a.energy().unwrap(), b.energy().unwrap());
    assert_eq!(d.positions(), original[1].coordinates());
    a.set_position_(AtomId::new(0), [0.1, 0., 0.]).unwrap();
    c.minimize_().unwrap();
    assert_eq!(mol.conformers_3d(), original);
    assert_eq!(clone.conformers_3d(), original);
    drop(mol);
    drop(clone);
    assert!(a.energy().unwrap().is_finite());
    assert!(c.energy().unwrap().is_finite());
}
#[test]
fn snapshots_cache_and_joint_evaluation() {
    let mol = fixture();
    let mut handle = mol.uff_force_field().unwrap();
    let mut snapshot = handle.positions();
    snapshot[1] = [2.3, 0.2, 0.];
    assert_ne!(handle.position(AtomId::new(1)).unwrap(), snapshot[1]);
    let before = handle.energy().unwrap();
    handle.set_positions_(&snapshot).unwrap();
    assert_ne!(handle.energy().unwrap(), before);
    let result = handle.energy_gradient().unwrap();
    assert_eq!(result.energy(), handle.energy().unwrap());
    assert_eq!(result.gradient(), handle.gradient().unwrap());
    handle
        .set_position_(AtomId::new(1), [2.4, 0.3, 0.])
        .unwrap();
    assert_eq!(snapshot[1], [2.3, 0.2, 0.]);
    assert_ne!(result.energy(), handle.energy().unwrap());
}
#[test]
fn unconstrained_gradient_is_registered_and_preserves_fixed_state() {
    let mol = fixture();
    for mut handle in [
        mol.mmff_force_field().unwrap(),
        mol.uff_force_field().unwrap(),
    ] {
        let expected = handle.gradient().unwrap();
        handle.set_fixed_atoms_(&[AtomId::new(0)]).unwrap();
        let before = handle.positions();
        let mut snapshot = handle.gradient_unconstrained().unwrap();
        assert_eq!(
            snapshot
                .iter()
                .flatten()
                .map(|x| x.to_bits())
                .collect::<Vec<_>>(),
            expected
                .iter()
                .flatten()
                .map(|x| x.to_bits())
                .collect::<Vec<_>>()
        );
        assert_ne!(snapshot[0], [0.; 3]);
        assert_eq!(handle.gradient().unwrap()[0], [0.; 3]);
        snapshot[0] = [123.; 3];
        assert_ne!(handle.gradient_unconstrained().unwrap()[0], snapshot[0]);
        assert_eq!(handle.fixed_atoms(), [AtomId::new(0)]);
        assert_eq!(handle.positions(), before);
    }
}
#[test]
fn invalid_updates_are_atomic_and_errors_retain_context() {
    let mut handle = fixture().mmff_force_field().unwrap();
    handle.set_fixed_atoms_(&[AtomId::new(1)]).unwrap();
    let original = handle.positions();
    let energy = handle.energy().unwrap();
    let e = handle.set_position_(AtomId::new(9), [0.; 3]).unwrap_err();
    assert_eq!(e.kind(), MolecularForceFieldErrorKind::InvalidAtomIndex);
    assert_eq!(e.atom_index(), Some(9));
    assert_eq!(e.expected(), Some(4));
    let e = handle
        .set_position_(AtomId::new(0), [0., f64::NAN, 0.])
        .unwrap_err();
    assert_eq!(e.component(), Some(1));
    let e = handle.set_positions_(&original[..3]).unwrap_err();
    assert_eq!(e.kind(), MolecularForceFieldErrorKind::CoordinateCount);
    assert_eq!(e.actual(), Some(3));
    assert_eq!(e.expected(), Some(4));
    let mut invalid = original.clone();
    invalid[3][2] = f64::INFINITY;
    assert!(handle.set_positions_(&invalid).is_err());
    assert!(
        handle
            .set_fixed_atoms_(&[AtomId::new(0), AtomId::new(99)])
            .is_err()
    );
    assert_eq!(handle.positions(), original);
    assert_eq!(handle.fixed_atoms(), [AtomId::new(1)]);
    assert_eq!(handle.energy().unwrap(), energy);
    let e = handle
        .minimize_with_params_(&ForceFieldMinimizeParams::new(10, f64::NAN, 1e-6))
        .unwrap_err();
    assert_eq!(e.kind(), MolecularForceFieldErrorKind::InvalidTolerance);
    assert_eq!(e.component(), Some(0));
    assert!(e.source().is_some());
    assert_eq!(handle.positions(), original);
}
#[test]
fn fixed_set_replaces_and_setters_update_anchors() {
    let mut handle = fixture().uff_force_field().unwrap();
    handle
        .set_fixed_atoms_(&[AtomId::new(1), AtomId::new(0), AtomId::new(1)])
        .unwrap();
    assert_eq!(handle.fixed_atoms(), [AtomId::new(0), AtomId::new(1)]);
    handle.set_fixed_atoms_(&[AtomId::new(0)]).unwrap();
    assert_eq!(handle.fixed_atoms(), [AtomId::new(0)]);
    handle
        .set_position_(AtomId::new(0), [0.2, 0.4, 0.6])
        .unwrap();
    assert_eq!(handle.gradient().unwrap()[0], [0.; 3]);
    let outcome = handle
        .minimize_with_params_(&ForceFieldMinimizeParams::new(5, 1e-4, 1e-6))
        .unwrap();
    assert!(outcome.iterations() <= 5);
    assert_eq!(outcome.energy(), handle.energy().unwrap());
    assert_eq!(handle.position(AtomId::new(0)).unwrap(), [0.2, 0.4, 0.6]);
    handle.set_fixed_atoms_(&[]).unwrap();
    assert!(handle.fixed_atoms().is_empty());
}
#[test]
fn minimization_returns_actual_count_and_restarts_history() {
    let mol = fixture();
    let mut first = mol.mmff_force_field().unwrap();
    let original = first.positions();
    let zero = first
        .minimize_with_params_(&ForceFieldMinimizeParams::new(0, 1e-4, 1e-6))
        .unwrap();
    assert_eq!(zero.iterations(), 0);
    assert!(!zero.converged());
    assert_eq!(first.positions(), original);
    let one = first
        .minimize_with_params_(&ForceFieldMinimizeParams::new(1, 1e-4, 1e-6))
        .unwrap();
    assert_eq!(one.iterations(), 1);
    let mut fresh = mol.mmff_force_field().unwrap();
    fresh.set_positions_(&first.positions()).unwrap();
    let options = ForceFieldMinimizeParams::new(2, 1e-4, 1e-6);
    let a = first.minimize_with_params_(&options).unwrap();
    let b = fresh.minimize_with_params_(&options).unwrap();
    assert_eq!(a, b);
    assert_eq!(first.positions(), fresh.positions());
}
#[test]
fn missing_conformer_and_invalid_parameterization_fail_without_embedding() {
    let mol = Molecule::from_smiles("CC").unwrap();
    assert!(mol.conformers_3d().is_empty());
    assert_eq!(
        mol.uff_force_field().unwrap_err().kind(),
        MolecularForceFieldErrorKind::MissingConformer
    );
    assert_eq!(
        mol.mmff_force_field().unwrap_err().kind(),
        MolecularForceFieldErrorKind::MissingConformer
    );
    assert!(mol.conformers_3d().is_empty());
    let mol = fixture();
    let uff = mol
        .uff_force_field_with_params(&UffForceFieldParams::new(Some(100), 10., true))
        .unwrap_err();
    assert_eq!(uff.requested(), Some(100));
    let mmff = mol
        .mmff_force_field_with_params(&MmffForceFieldParams::new(
            Some(100),
            "MMFF94".into(),
            100.,
            true,
        ))
        .unwrap_err();
    assert_eq!(mmff.kind(), MolecularForceFieldErrorKind::MissingConformer);
    assert_eq!(mmff.requested(), Some(100));
    // Pinned AtomTyper.cpp constructor: d_mmffs(mmffVariant == "MMFF94s").
    let fallback = mol
        .mmff_force_field_with_params(&MmffForceFieldParams::new(
            None,
            "unknown".into(),
            100.,
            true,
        ))
        .unwrap();
    assert_eq!(
        fallback.energy().unwrap(),
        mol.mmff_force_field().unwrap().energy().unwrap()
    );
    let helium = Molecule::from_smiles("[He]").unwrap();
    let helium = Molecule::from_parts(
        helium.topology().clone(),
        CoordinateBlock {
            conformers_3d: vec![Conformer3D::new(0, vec![[0.; 3]], true)],
            ..Default::default()
        },
        helium.properties().clone(),
    )
    .unwrap();
    let mmff = helium.mmff_force_field().unwrap_err();
    assert_eq!(
        mmff.kind(),
        MolecularForceFieldErrorKind::InvalidParameterization
    );
    assert!(mmff.source().is_some());
}
