use cosmolkit_experiment_model::{
    Atom, BondOrder, Conformer3D, CoordinateBlock, Element, MoleculeProperties, TopologyBlock,
};

pub fn add_hydrogens(
    mut topology: TopologyBlock,
    coordinates: CoordinateBlock,
    properties: MoleculeProperties,
) -> Result<(TopologyBlock, CoordinateBlock, MoleculeProperties), String> {
    // Architecture-only behavior: the production source port will live here.
    let existing = topology.atom_count();
    let hydrogen = topology.add_atom(Atom::new(Element::H));
    if existing > 0 {
        topology
            .add_bond(
                cosmolkit_experiment_model::AtomId::from_index(existing - 1),
                hydrogen,
                BondOrder::Single,
            )
            .map_err(|error| error.to_string())?;
    }
    Ok((topology, coordinates, properties))
}

pub fn remove_hydrogens(
    topology: TopologyBlock,
    coordinates: CoordinateBlock,
    properties: MoleculeProperties,
) -> Result<(TopologyBlock, CoordinateBlock, MoleculeProperties), String> {
    Ok((topology, coordinates, properties))
}

pub fn enumerate_tautomers(
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

pub fn enumerate_stereoisomers(
    topology: &TopologyBlock,
    coordinates: &CoordinateBlock,
    properties: &MoleculeProperties,
) -> Result<Vec<(TopologyBlock, CoordinateBlock, MoleculeProperties)>, String> {
    enumerate_tautomers(topology, coordinates, properties)
}

pub fn generate_conformer(topology: &TopologyBlock) -> Result<Conformer3D, String> {
    let coordinates = (0..topology.atom_count())
        .map(|index| [index as f64, 0.0, 0.0])
        .collect();
    Ok(Conformer3D::new(coordinates))
}
