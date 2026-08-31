use cosmolkit_experiment_model::{CoordinateBlock, MoleculeProperties, TopologyBlock};

pub fn enumerate_stereoisomers(
    topology: &TopologyBlock,
    coordinates: &CoordinateBlock,
    properties: &MoleculeProperties,
) -> Result<Vec<(TopologyBlock, CoordinateBlock, MoleculeProperties)>, String> {
    Ok(vec![(
        topology.clone(),
        coordinates.clone(),
        properties.clone(),
    )])
}
