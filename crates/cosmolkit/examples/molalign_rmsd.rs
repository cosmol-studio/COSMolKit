use cosmolkit::{AlignmentAtomMap, AlignmentParameters, Molecule};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let reference_coordinates = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0]];
    let probe_coordinates = reference_coordinates
        .iter()
        .map(|point| [point[0] + 3.0, point[1] - 2.0, point[2] + 1.0])
        .collect();

    let reference = Molecule::from_smiles("CCC")?.with_only_3d_conformer(reference_coordinates, true)?;
    let probe = Molecule::from_smiles("CCC")?.with_only_3d_conformer(probe_coordinates, true)?;
    let params = AlignmentParameters {
        atom_map: Some(
            (0..3)
                .map(|index| AlignmentAtomMap {
                    probe_atom: index,
                    reference_atom: index,
                })
                .collect(),
        ),
        ..AlignmentParameters::default()
    };

    let measurement = probe.alignment_transform_to(&reference, &params)?;
    let (aligned, applied) = probe.with_alignment_to(&reference, &params)?;

    println!("read-only RMSD: {}", measurement.rmsd);
    println!("applied RMSD: {}", applied.rmsd);
    println!(
        "source first coordinate: {:?}",
        probe.conformers_3d()[0].coordinates()[0]
    );
    println!(
        "aligned first coordinate: {:?}",
        aligned.conformers_3d()[0].coordinates()[0]
    );
    Ok(())
}
