//! Shared source-backed chemistry implementations for the crate-split experiment.
//!
//! This crate operates on detached model values. It deliberately does not
//! define `Molecule`, operation contracts, cache authority, or public facade
//! methods. Closely coupled foundational chemistry such as valence handling,
//! sanitization, and hydrogen transforms belongs here; independently useful
//! families may live in sibling domain crates.

use cosmolkit_experiment_model::{
    Atom, AtomId, BondOrder, CoordinateBlock, Element, MoleculeProperties, TopologyBlock,
};

/// Architecture-only hydrogen transform boundary.
pub fn add_hydrogens(
    mut topology: TopologyBlock,
    coordinates: CoordinateBlock,
    properties: MoleculeProperties,
) -> Result<(TopologyBlock, CoordinateBlock, MoleculeProperties), String> {
    let existing = topology.atom_count();
    let hydrogen = topology.add_atom(Atom::new(Element::H));
    if existing > 0 {
        topology
            .add_bond(
                AtomId::from_index(existing - 1),
                hydrogen,
                BondOrder::Single,
            )
            .map_err(|error| error.to_string())?;
    }
    Ok((topology, coordinates, properties))
}

/// Architecture-only removal boundary; the source-backed implementation will
/// share this core crate with sanitization and valence handling.
pub fn remove_hydrogens(
    topology: TopologyBlock,
    coordinates: CoordinateBlock,
    properties: MoleculeProperties,
) -> Result<(TopologyBlock, CoordinateBlock, MoleculeProperties), String> {
    Ok((topology, coordinates, properties))
}
