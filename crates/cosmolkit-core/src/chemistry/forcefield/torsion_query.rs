//! Shared source-backed torsion-bond query selection for RDKit force fields.

use crate::search::smarts_parse::{SmartsParseParams, mol_from_smarts};
use crate::{Molecule, SubstructMatchError, SubstructMatchParams};

pub(crate) const DEFAULT_TORSION_BOND_SMARTS: &str = "[!$(*#*)&!D1]~[!$(*#*)&!D1]";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct TorsionBondMatch {
    pub(crate) begin_atom_index: usize,
    pub(crate) end_atom_index: usize,
    pub(crate) bond_index: usize,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub(crate) enum TorsionBondQueryError {
    #[error("failed to build torsion bond SMARTS query {smarts:?}: {detail}")]
    QueryBuild { smarts: String, detail: String },
    #[error(
        "torsion bond SMARTS must describe exactly two atoms joined by one bond, got {atoms} atoms and {bonds} bonds: {smarts:?}"
    )]
    QueryShape { smarts: String, atoms: usize, bonds: usize },
    #[error(transparent)]
    SubstructMatch(#[from] SubstructMatchError),
    #[error("torsion bond SMARTS match did not map both query atoms: {smarts:?}")]
    IncompleteAtomMapping { smarts: String },
    #[error("torsion bond SMARTS matched two atoms without a molecular bond: {smarts:?}")]
    MissingMatchedBond { smarts: String },
}

pub(crate) fn match_torsion_bonds(
    mol: &Molecule,
    torsion_bond_smarts: &str,
) -> Result<Vec<TorsionBondMatch>, TorsionBondQueryError> {
    // BEGIN RDKIT CPP HELPER UFF/MMFF::Tools::addTorsions query selection
    // RDKit✔️✔️:   std::vector<MatchVectType> matchVect;
    // RDKit✔️✔️:   const ROMol *defaultQuery = DefaultTorsionBondSmarts::query();
    // RDKit✔️✔️:   const ROMol *query = (torsionBondSmarts == DefaultTorsionBondSmarts::string())
    // RDKit✔️✔️:                            ? defaultQuery
    // RDKit✔️✔️:                            : SmartsToMol(torsionBondSmarts);
    // RDKit✔️✔️:   TEST_ASSERT(query);
    // RDKit✔️✔️:   unsigned int nHits = SubstructMatch(mol, *query, matchVect);
    // RDKit✔️✔️:   if (query != defaultQuery) {
    // RDKit✔️✔️:     delete query;
    // RDKit✔️✔️:   }
    // END RDKIT CPP HELPER UFF/MMFF::Tools::addTorsions query selection
    if torsion_bond_smarts == DEFAULT_TORSION_BOND_SMARTS {
        return Ok(match_default_torsion_bonds(mol));
    }
    match_parsed_torsion_bonds(mol, torsion_bond_smarts)
}

fn match_default_torsion_bonds(mol: &Molecule) -> Vec<TorsionBondMatch> {
    // The fixed two-atom query is symmetric. This direct traversal preserves
    // RDKit's first accepted orientation, atom-set uniquification, and default
    // maxMatches=1000 without constructing a general query graph.
    const MAX_MATCHES: usize = 1000;
    let mut matches = Vec::new();
    for atom_index in 0..mol.num_atoms() {
        if !matches_default_torsion_atom(mol, atom_index) {
            continue;
        }
        for neighbor in mol.topology_block().adjacency.neighbors_of(atom_index) {
            if neighbor.atom_index <= atom_index || !matches_default_torsion_atom(mol, neighbor.atom_index) {
                continue;
            }
            matches.push(TorsionBondMatch {
                begin_atom_index: atom_index,
                end_atom_index: neighbor.atom_index,
                bond_index: neighbor.bond.index(),
            });
            if matches.len() == MAX_MATCHES {
                return matches;
            }
        }
    }
    matches
}

fn matches_default_torsion_atom(mol: &Molecule, atom_index: usize) -> bool {
    let neighbors = mol.topology_block().adjacency.neighbors_of(atom_index);
    neighbors.len() != 1
        && neighbors
            .iter()
            .all(|neighbor| mol.bonds()[neighbor.bond.index()].order() != crate::BondOrder::Triple)
}

fn match_parsed_torsion_bonds(
    mol: &Molecule,
    torsion_bond_smarts: &str,
) -> Result<Vec<TorsionBondMatch>, TorsionBondQueryError> {
    // RDKit✔️✔️:   const ROMol *query = (torsionBondSmarts == DefaultTorsionBondSmarts::string())
    // RDKit✔️✔️:                            ? defaultQuery
    // RDKit✔️✔️:                            : SmartsToMol(torsionBondSmarts);
    // RDKit✔️✔️:   TEST_ASSERT(query);
    // RDKit✔️✔️:   unsigned int nHits = SubstructMatch(mol, *query, matchVect);
    // Local complexity review: canonical compilation is linear in the SMARTS
    // input and matching retains the existing VF2 graph traversal. This path
    // no longer performs a compatibility conversion or alternate parse.
    let query = mol_from_smarts(torsion_bond_smarts, &SmartsParseParams::default()).map_err(|error| {
        TorsionBondQueryError::QueryBuild {
            smarts: torsion_bond_smarts.to_owned(),
            detail: error.to_string(),
        }
    })?;
    if query.num_atoms() != 2 || query.num_bonds() != 1 {
        return Err(TorsionBondQueryError::QueryShape {
            smarts: torsion_bond_smarts.to_owned(),
            atoms: query.num_atoms(),
            bonds: query.num_bonds(),
        });
    }

    crate::try_get_substruct_matches_with_params(mol, &query, &SubstructMatchParams::default())?
        .into_iter()
        .map(|matched| match_to_torsion_bond(mol, torsion_bond_smarts, &matched.atom_mapping))
        .collect()
}

fn match_to_torsion_bond(
    mol: &Molecule,
    torsion_bond_smarts: &str,
    atom_mapping: &[usize],
) -> Result<TorsionBondMatch, TorsionBondQueryError> {
    let [begin_atom_index, end_atom_index] = atom_mapping else {
        return Err(TorsionBondQueryError::IncompleteAtomMapping {
            smarts: torsion_bond_smarts.to_owned(),
        });
    };
    let bond_index = mol
        .topology_block()
        .adjacency
        .neighbors_of(*begin_atom_index)
        .iter()
        .find(|neighbor| neighbor.atom_index == *end_atom_index)
        .map(|neighbor| neighbor.bond.index())
        .ok_or_else(|| TorsionBondQueryError::MissingMatchedBond {
            smarts: torsion_bond_smarts.to_owned(),
        })?;
    Ok(TorsionBondMatch {
        begin_atom_index: *begin_atom_index,
        end_atom_index: *end_atom_index,
        bond_index,
    })
}

#[cfg(test)]
mod tests {
    use super::{DEFAULT_TORSION_BOND_SMARTS, TorsionBondMatch, TorsionBondQueryError, match_torsion_bonds};
    use crate::{AtomSpec, BondOrder, BondSpec, Element, Molecule, MoleculeBuilder};

    fn chain(bond_orders: &[BondOrder]) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..=bond_orders.len())
            .map(|_| builder.add_atom(AtomSpec::new(Element::C)))
            .collect();
        for (atom_index, &bond_order) in bond_orders.iter().enumerate() {
            builder
                .add_bond(BondSpec::new(atoms[atom_index], atoms[atom_index + 1], bond_order))
                .expect("chain bond should build");
        }
        builder.build().expect("chain should build")
    }

    #[test]
    fn default_literal_matches_rdkit_2026_03_1() {
        assert_eq!(DEFAULT_TORSION_BOND_SMARTS, "[!$(*#*)&!D1]~[!$(*#*)&!D1]");
    }

    #[test]
    fn default_query_matches_only_the_internal_chain_bond() {
        let matches = match_torsion_bonds(&chain(&[BondOrder::Single; 3]), DEFAULT_TORSION_BOND_SMARTS)
            .expect("default query should match");

        assert_eq!(
            matches,
            vec![TorsionBondMatch {
                begin_atom_index: 1,
                end_atom_index: 2,
                bond_index: 1,
            }]
        );
    }

    #[test]
    fn default_query_preserves_rdkit_ring_match_order() {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..4).map(|_| builder.add_atom(AtomSpec::new(Element::C))).collect();
        for (begin, end) in [(0, 1), (1, 2), (2, 3), (3, 0)] {
            builder
                .add_bond(BondSpec::new(atoms[begin], atoms[end], BondOrder::Single))
                .expect("ring bond should build");
        }
        let mol = builder.build().expect("ring should build");

        let matches = match_torsion_bonds(&mol, DEFAULT_TORSION_BOND_SMARTS).expect("default query should match");
        let atom_pairs: Vec<_> = matches
            .iter()
            .map(|matched| (matched.begin_atom_index, matched.end_atom_index))
            .collect();

        assert_eq!(atom_pairs, vec![(0, 1), (0, 3), (1, 2), (2, 3)]);
    }

    #[test]
    fn default_query_excludes_atoms_adjacent_to_triple_bonds() {
        let matches = match_torsion_bonds(
            &chain(&[BondOrder::Triple, BondOrder::Single, BondOrder::Single]),
            DEFAULT_TORSION_BOND_SMARTS,
        )
        .expect("default query should evaluate");

        assert!(matches.is_empty());
    }

    #[test]
    fn smarts_consumer_torsion_query() {
        let matches = match_torsion_bonds(
            &chain(&[BondOrder::Single, BondOrder::Double, BondOrder::Single]),
            "[*:1]-[*:2]",
        )
        .expect("custom query should match only single bonds");
        let atom_pairs: Vec<_> = matches
            .iter()
            .map(|matched| (matched.begin_atom_index, matched.end_atom_index))
            .collect();

        assert_eq!(atom_pairs, vec![(0, 1), (2, 3)]);
    }

    #[test]
    fn default_query_stops_at_rdkit_default_match_limit() {
        let mol = chain(&vec![BondOrder::Single; 1004]);

        let matches =
            match_torsion_bonds(&mol, DEFAULT_TORSION_BOND_SMARTS).expect("default query should match the long chain");

        assert_eq!(matches.len(), 1000);
        assert_eq!(matches.first().map(|matched| matched.bond_index), Some(1));
        assert_eq!(matches.last().map(|matched| matched.bond_index), Some(1000));
    }

    #[test]
    fn custom_query_rejects_non_bond_shapes_explicitly() {
        let error = match_torsion_bonds(&chain(&[BondOrder::Single; 3]), "[*:1]-[*:2]-[*:3]")
            .expect_err("three-atom query is not a torsion-bond selector");

        assert!(matches!(
            error,
            TorsionBondQueryError::QueryShape { atoms: 3, bonds: 2, .. }
        ));
    }

    #[test]
    fn invalid_custom_query_preserves_parse_failure() {
        let error = match_torsion_bonds(&chain(&[BondOrder::Single]), "[")
            .expect_err("invalid SMARTS must not be treated as no matches");

        assert!(matches!(error, TorsionBondQueryError::QueryBuild { .. }));
    }
}
