use cosmolkit_experiment_model::{Conformer3D, TopologyBlock};

pub fn generate_conformer(topology: &TopologyBlock) -> Result<Conformer3D, String> {
    let coordinates = (0..topology.atom_count())
        .map(|index| [index as f64, 0.0, 0.0])
        .collect();
    Ok(Conformer3D::new(coordinates))
}
