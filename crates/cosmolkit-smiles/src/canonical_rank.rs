// RDKit marker convention defined in dev/source_reproduction_protocol.md.

use cosmolkit_core::{CanonicalRankError, rank_fragment_atoms};
use cosmolkit_model::{AdjacencyList, AtomId, BondId, TopologyBlock};

pub(crate) fn rank_component_atoms(
    topology: &TopologyBlock,
    component_atoms: &[usize],
) -> Result<Vec<usize>, CanonicalRankError> {
    // BEGIN RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles fragment ranking boundary
    // RDKit✔️✔️:   auto mols =
    // RDKit✔️✔️:       MolOps::getMolFrags(mol, false, nullptr, &fragsMolAtomMapping, false);
    // RDKit✔️✔️:   for (unsigned fragIdx = 0; fragIdx < mols.size(); fragIdx++) {
    // RDKit✔️✔️:     ROMol *tmol = mols[fragIdx].get();
    // RDKit✔️✔️:     if (params.canonical) {
    // RDKit✔️✔️:       Canon::rankMolAtoms(*tmol, ranks, breakTies, includeChirality,
    // RDKit✔️✔️:                           includeIsotopes, includeAtomMaps,
    // RDKit✔️✔️:                           includeChiralPresence, includeStereoGroups,
    // RDKit✔️✔️:                           useNonStereoRanks);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION SmilesWrite::detail::MolToSmiles fragment ranking boundary
    let mut old_to_new = vec![None; topology.atoms.len()];
    let mut atoms = Vec::with_capacity(component_atoms.len());
    for (new_index, &old_index) in component_atoms.iter().enumerate() {
        old_to_new[old_index] = Some(AtomId::new(new_index));
        atoms.push(
            topology.atoms[old_index]
                .clone()
                .with_id(AtomId::new(new_index)),
        );
    }
    let mut bonds = Vec::new();
    for bond in &topology.bonds {
        let (Some(begin), Some(end)) = (
            old_to_new[bond.begin().index()],
            old_to_new[bond.end().index()],
        ) else {
            continue;
        };
        let stereo_atoms = bond.stereo_atoms().and_then(|stereo_atoms| {
            Some([
                old_to_new[stereo_atoms[0].index()]?,
                old_to_new[stereo_atoms[1].index()]?,
            ])
        });
        bonds.push(
            bond.clone()
                .remapped(BondId::new(bonds.len()), begin, end, stereo_atoms),
        );
    }
    let fragment = TopologyBlock {
        adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
        atoms,
        bonds,
        ..TopologyBlock::default()
    };
    let atom_mask = vec![true; fragment.atoms.len()];
    let bond_mask = vec![true; fragment.bonds.len()];
    let local_ranks = rank_fragment_atoms(&fragment, &atom_mask, &bond_mask)?;
    let mut ranks = vec![usize::MAX; topology.atoms.len()];
    for (new_index, &old_index) in component_atoms.iter().enumerate() {
        ranks[old_index] = local_ranks[new_index];
    }
    Ok(ranks)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parse_smiles;

    #[test]
    fn detached_rank_mol_port_matches_rdkit_canonical_rank_atoms() {
        for (smiles, expected) in [
            ("OCC", vec![1, 2, 0]),
            ("OC(C)C", vec![2, 3, 0, 1]),
            ("c1ccncc1", vec![0, 1, 3, 5, 4, 2]),
            ("C12(CCCCC1)CCCCC2", vec![10, 6, 2, 0, 3, 7, 8, 4, 1, 5, 9]),
            ("C1C2C3C1C2C3", vec![0, 2, 4, 3, 5, 1]),
            ("[13CH3]CO", vec![0, 2, 1]),
            ("[CH3:7]CO", vec![2, 1, 0]),
        ] {
            let record = parse_smiles(smiles, &Default::default()).unwrap();
            let atoms = (0..record.topology.atoms.len()).collect::<Vec<_>>();
            assert_eq!(
                rank_component_atoms(&record.topology, &atoms).unwrap(),
                expected,
                "{smiles}"
            );
        }
    }
}
