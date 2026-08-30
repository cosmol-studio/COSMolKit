use cosmolkit_experiment_model::{
    Atom, BondOrder, CoordinateBlock, Element, MoleculeProperties, QueryAst, TopologyBlock,
};

fn parts_from_atom_count(count: usize) -> (TopologyBlock, CoordinateBlock, MoleculeProperties) {
    let mut topology = TopologyBlock::new();
    let mut previous = None;
    for _ in 0..count {
        let current = topology.add_atom(Atom::new(Element::C));
        if let Some(previous) = previous {
            let _ = topology.add_bond(previous, current, BondOrder::Single);
        }
        previous = Some(current);
    }
    (
        topology,
        CoordinateBlock::default(),
        MoleculeProperties::new(),
    )
}

pub fn parse_smiles(
    source: &str,
) -> Result<(TopologyBlock, CoordinateBlock, MoleculeProperties), String> {
    if source.is_empty() {
        return Err("empty SMILES".to_owned());
    }
    Ok(parts_from_atom_count(
        source
            .chars()
            .filter(|character| *character == 'C')
            .count()
            .max(1),
    ))
}

pub fn parse_smarts(source: &str) -> Result<QueryAst, String> {
    if source.is_empty() {
        return Err("empty SMARTS".to_owned());
    }
    Ok(QueryAst::new(source))
}

pub fn read_sdf(
    source: &str,
) -> Result<Vec<(TopologyBlock, CoordinateBlock, MoleculeProperties)>, String> {
    if source.is_empty() {
        return Err("empty SDF".to_owned());
    }
    Ok(vec![parts_from_atom_count(1)])
}

pub fn read_pdb(
    source: &str,
) -> Result<(TopologyBlock, CoordinateBlock, MoleculeProperties), String> {
    if source.is_empty() {
        return Err("empty PDB".to_owned());
    }
    Ok(parts_from_atom_count(1))
}

pub fn read_xyz(
    source: &str,
) -> Result<(TopologyBlock, CoordinateBlock, MoleculeProperties), String> {
    if source.is_empty() {
        return Err("empty XYZ".to_owned());
    }
    Ok(parts_from_atom_count(1))
}
