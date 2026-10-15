// Run: cargo run -p cosmolkit --example forcefield_optimization
// Requires the conformer capability bundle (included in default full).
use cosmolkit::{
    AtomId, ForceFieldAngleConstraintParams, ForceFieldDistanceConstraintParams,
    ForceFieldPositionConstraintParams, ForceFieldTorsionConstraintParams, MmffForceFieldParams,
    MmffOptimizationParams, Molecule, UffOptimizationParams,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let molecule = Molecule::from_smiles("CCO")?.with_hydrogens()?.sanitize()?;

    let molecule = molecule.with_added_3d_conformer(vec![
        vec![0.000, 0.000, 0.000],
        vec![1.540, 0.000, 0.000],
        vec![2.100, 1.200, 0.000],
        vec![-0.600, 0.900, 0.000],
        vec![-0.600, -0.900, 0.000],
        vec![0.000, 0.000, 1.000],
        vec![1.900, -0.900, 0.000],
        vec![1.700, 0.000, 1.000],
        vec![2.900, 1.200, 0.000],
    ])?;

    if molecule.uff_has_all_molecule_params()? {
        let result = molecule.with_uff_optimized_with_params(&UffOptimizationParams {
            max_iterations: 200,
            vdw_threshold: 10.0,
            conformer_id: None,
            ignore_interfragment_interactions: true,
        })?;
        println!(
            "UFF needs_more={} energy={:.6}",
            result.status > 0,
            result.energy
        );
        println!(
            "optimized first atom: {:?}",
            result.molecule.conformers_3d()[0].coordinates()[0]
        );
    }

    if molecule.mmff_has_all_molecule_params()? {
        let result = molecule.with_mmff_optimized_with_params(&MmffOptimizationParams {
            mmff_variant: "MMFF94".into(),
            max_iterations: 200,
            non_bonded_threshold: 100.0,
            conformer_id: None,
            ignore_interfragment_interactions: true,
        })?;
        println!("MMFF94 needs_more={}", result.needs_more());
        println!(
            "optimized first atom: {:?}",
            result.molecule.conformers_3d()[0].coordinates()[0]
        );
    }

    // Persistent evaluators own their coordinates; they never write back to molecule.
    // Every constraint value below is an explicit example choice, not an API default.
    let original = molecule.conformers_3d()[0].coordinates().to_vec();
    for model in ["UFF", "MMFF94", "MMFF94s"] {
        let mut field = if model == "UFF" {
            molecule.uff_force_field()?
        } else {
            molecule.mmff_force_field_with_params(&MmffForceFieldParams::new(
                None,
                model.into(),
                100.0,
                true,
            ))?
        };
        // Absolute C-C distance limits in Angstroms.
        field.add_distance_constraint_(
            AtomId::new(0),
            AtomId::new(1),
            &ForceFieldDistanceConstraintParams::new(false, 1.4, 1.6, 100.0),
        )?;
        // Relative angular limits in degrees capture the geometry now, not at minimize().
        field.add_angle_constraint_(
            AtomId::new(0),
            AtomId::new(1),
            AtomId::new(2),
            &ForceFieldAngleConstraintParams::new(true, -5.0, 5.0, 1.0),
        )?;
        field.add_torsion_constraint_(
            AtomId::new(3),
            AtomId::new(0),
            AtomId::new(1),
            AtomId::new(2),
            &ForceFieldTorsionConstraintParams::new(true, -10.0, 10.0, 1.0),
        )?;
        // A soft position restraint anchors O here, with a 0.25 Angstrom allowance.
        field.add_position_constraint_(
            AtomId::new(2),
            &ForceFieldPositionConstraintParams::new(0.25, 10.0),
        )?;
        // A fixed atom is different: its coordinates cannot move during minimization.
        field.set_fixed_atoms_(&[AtomId::new(0)])?;
        let pinned = field.position(AtomId::new(0))?;
        let oxygen = field.position(AtomId::new(2))?;
        field.set_position_(
            AtomId::new(2),
            [oxygen[0] + 0.3, oxygen[1] + 0.1, oxygen[2] + 0.2],
        )?;
        let initial = field.energy_gradient()?;
        assert_eq!(initial.gradient().len(), molecule.num_atoms());
        let result = field.minimize_()?;
        assert!(result.energy().is_finite());
        assert!(result.energy() <= initial.energy());
        assert_eq!(field.position(AtomId::new(0))?, pinned);
        assert_eq!(
            molecule.conformers_3d()[0].coordinates(),
            original.as_slice()
        );

        // Explicit value-style coordinate installation creates a new molecule.
        let positions = field.positions();
        let relaxed =
            molecule.with_3d_coordinates(positions.iter().map(|row| row.to_vec()).collect())?;
        assert_eq!(
            relaxed.conformers_3d()[0].coordinates(),
            positions.as_slice()
        );
        assert_eq!(
            molecule.conformers_3d()[0].coordinates(),
            original.as_slice()
        );
        println!(
            "{model} restrained: energy={:.6} -> {:.6}, converged={}, iterations={}",
            initial.energy(),
            result.energy(),
            result.converged(),
            result.iterations()
        );
    }
    Ok(())
}
