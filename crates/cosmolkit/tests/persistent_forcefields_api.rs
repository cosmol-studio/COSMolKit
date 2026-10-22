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

#[test]
fn persistent_constraints_are_owned_registered_and_failure_preserving() {
    use cosmolkit::{
        ForceFieldAngleConstraintParams, ForceFieldDistanceConstraintParams,
        ForceFieldPositionConstraintParams, ForceFieldTorsionConstraintParams,
    };
    let mol = fixture();
    let source_rows = mol.conformers_3d()[0].coordinates().to_vec();
    let mut field = mol.uff_force_field().unwrap();
    let other = mol.uff_force_field().unwrap();
    let initial = other.energy().unwrap();
    field
        .add_distance_constraint_(
            AtomId::new(0),
            AtomId::new(1),
            &ForceFieldDistanceConstraintParams::new(false, 1., 1., 100.),
        )
        .unwrap();
    field
        .add_angle_constraint_(
            AtomId::new(0),
            AtomId::new(1),
            AtomId::new(2),
            &ForceFieldAngleConstraintParams::new(true, -1., 1., 10.),
        )
        .unwrap();
    field
        .add_torsion_constraint_(
            AtomId::new(0),
            AtomId::new(1),
            AtomId::new(2),
            AtomId::new(3),
            &ForceFieldTorsionConstraintParams::new(true, -1., 1., 10.),
        )
        .unwrap();
    field
        .add_position_constraint_(
            AtomId::new(0),
            &ForceFieldPositionConstraintParams::new(0.2, 100.),
        )
        .unwrap();
    let before = field.energy().unwrap();
    assert!(before > initial);
    let error = field
        .add_distance_constraint_(
            AtomId::new(0),
            AtomId::new(99),
            &ForceFieldDistanceConstraintParams::new(false, 1., 2., 10.),
        )
        .unwrap_err();
    assert_eq!(
        (
            error.kind(),
            error.atom_index(),
            error.component(),
            error.expected()
        ),
        (
            MolecularForceFieldErrorKind::InvalidAtomIndex,
            Some(99),
            Some(1),
            Some(4)
        )
    );
    assert!(error.source().is_some());
    assert_eq!(field.energy().unwrap().to_bits(), before.to_bits());
    let error = field
        .add_position_constraint_(
            AtomId::new(usize::MAX),
            &ForceFieldPositionConstraintParams::new(0., 1.),
        )
        .unwrap_err();
    assert_eq!(error.atom_index(), Some(usize::MAX));
    assert_eq!(field.energy().unwrap().to_bits(), before.to_bits());
    // Source validates the first representable out-of-range ID before a later
    // native usize ID that cannot fit source u32. Torsion checks bounds first.
    let invalid = AtomId::new(mol.num_atoms());
    let huge = AtomId::new(usize::MAX);
    let error = field
        .add_distance_constraint_(
            invalid,
            huge,
            &ForceFieldDistanceConstraintParams::new(false, 5., 1., 1.),
        )
        .unwrap_err();
    assert_eq!((error.atom_index(), error.component()), (Some(4), Some(0)));
    assert!(error.source().is_some());
    let error = field
        .add_angle_constraint_(
            invalid,
            huge,
            AtomId::new(0),
            &ForceFieldAngleConstraintParams::new(false, 190., 1., 1.),
        )
        .unwrap_err();
    assert_eq!((error.atom_index(), error.component()), (Some(4), Some(0)));
    let error = field
        .add_torsion_constraint_(
            invalid,
            huge,
            AtomId::new(0),
            AtomId::new(1),
            &ForceFieldTorsionConstraintParams::new(false, 5., 1., 1.),
        )
        .unwrap_err();
    assert_eq!(
        error.kind(),
        MolecularForceFieldErrorKind::InvalidConstraintBounds
    );
    let error = field
        .add_torsion_constraint_(
            invalid,
            huge,
            AtomId::new(0),
            AtomId::new(1),
            &ForceFieldTorsionConstraintParams::new(false, 1., 5., 1.),
        )
        .unwrap_err();
    assert_eq!((error.atom_index(), error.component()), (Some(4), Some(0)));
    assert_eq!(field.energy().unwrap().to_bits(), before.to_bits());
    field.set_position_(AtomId::new(0), [1., 0., 0.]).unwrap();
    assert_eq!(mol.conformers_3d()[0].coordinates(), source_rows);
    assert_eq!(other.energy().unwrap().to_bits(), initial.to_bits());
}

#[test]
#[cfg(feature = "cap-hydrogens")]
fn uff_optimizer_returns_the_energy_of_its_own_force_field() {
    use cosmolkit::UffOptimizationParams;
    // UFF.h::UFFOptimizeMolecule / FFConvenience.h::OptimizeMolecule returns
    // the minimizing field's calcEnergy(), not a rebuilt field's energy.
    let mut mol = Molecule::from_smiles("CP(Br)(c1ccccc1)(c1ccccc1)c1ccccc1")
        .unwrap()
        .with_hydrogens()
        .unwrap();
    let mut state = 17_u64;
    let positions = (0..mol.num_atoms())
        .map(|_| {
            std::array::from_fn(|_| {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                (((state >> 32) % 257) as f64 - 128.) / 16.
            })
        })
        .collect();
    mol = mol.with_added_3d_conformer(positions).unwrap();
    for expected in [0x42e342814525ed28_u64, 0x42b30f2267ad645e] {
        let result = mol
            .with_uff_optimized_with_params(&UffOptimizationParams {
                max_iterations: 2,
                vdw_threshold: 100.,
                ..Default::default()
            })
            .unwrap();
        assert_eq!(result.energy().unwrap().to_bits(), expected);
        mol = result.molecule;
    }
}

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
#[cfg(feature = "cap-hydrogens")]
fn mmff_observation_does_not_precondition_the_next_optimizer() {
    use cosmolkit::{MmffOptimizationParams, MmffPropertiesParams};
    // AtomTyper.cpp::MMFFMolProperties prepares aromaticity on its input.
    // Queries/factories must keep that preparation detached: on this fixed
    // macrocycle, repeated preparation is observably non-idempotent upstream.
    let mol = Molecule::from_smiles(concat!(
        "CC1=C(/C=C2\\C(C)=C3/C(C)=C(/C=C4\\C(C)=C5/C(CCC(=O)O)=",
        "C(C(C)=C5N4)C=C6\\C(CCC(=O)O)=C(C)C(=C6N2)C=C1)N3)C=C(\\C=C)C"
    ))
    .unwrap()
    .with_hydrogens()
    .unwrap();
    let mut state = 17_u64;
    let coordinates = (0..mol.num_atoms())
        .map(|_| {
            std::array::from_fn(|_| {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                (((state >> 32) % 257) as f64 - 128.) / 16.
            })
        })
        .collect();
    let mol = mol.with_added_3d_conformer(coordinates).unwrap();
    for variant in ["MMFF94", "MMFF94s"] {
        let mut queried = mol.clone();
        let mut direct = mol.clone();
        for _ in 0..2 {
            let before = queried.clone();
            queried
                .mmff_properties_with_params(&MmffPropertiesParams {
                    mmff_variant: cosmolkit::MmffVariant::from_rdkit_constructor_arg(variant),
                })
                .unwrap();
            queried
                .mmff_force_field_with_params(&MmffForceFieldParams::new(
                    None,
                    cosmolkit::MmffVariant::from_rdkit_constructor_arg(variant),
                    100.,
                    true,
                ))
                .unwrap()
                .energy()
                .unwrap();
            assert_eq!(queried, before, "{variant}: observing must not write back");
            let params = MmffOptimizationParams {
                mmff_variant: cosmolkit::MmffVariant::from_rdkit_constructor_arg(variant),
                max_iterations: 2,
                ..Default::default()
            };
            queried = queried
                .with_mmff_optimized_with_params(&params)
                .unwrap()
                .molecule;
            direct = direct
                .with_mmff_optimized_with_params(&params)
                .unwrap()
                .molecule;
            assert_eq!(queried.topology(), direct.topology());
            assert_eq!(queried.properties(), direct.properties());
            let bits = |m: &Molecule| {
                m.conformers_3d()[0]
                    .coordinates()
                    .iter()
                    .flatten()
                    .map(|v| v.to_bits())
                    .collect::<Vec<_>>()
            };
            assert_eq!(bits(&queried), bits(&direct));
        }
    }
}
#[test]
fn immutable_defaults_and_registered_properties() {
    let mmff = MmffForceFieldParams::default();
    assert_eq!(mmff.conformer_id(), None);
    assert_eq!(mmff.mmff_variant(), cosmolkit::MmffVariant::Mmff94);
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
            cosmolkit::MmffVariant::Mmff94,
            100.,
            true,
        ))
        .unwrap_err();
    assert_eq!(mmff.kind(), MolecularForceFieldErrorKind::MissingConformer);
    assert_eq!(mmff.requested(), Some(100));
    // The public configuration represents the default variant explicitly.
    // RDKit's arbitrary-string fallback is confined to the source constructor.
    let fallback = mol
        .mmff_force_field_with_params(&MmffForceFieldParams::new(
            None,
            cosmolkit::MmffVariant::Mmff94,
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
