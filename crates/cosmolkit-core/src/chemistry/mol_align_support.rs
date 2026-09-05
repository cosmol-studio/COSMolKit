//! Shared clone-only query adaptation used by MolAlign and distgeom pruning.

use std::{collections::BTreeMap, sync::OnceLock};

use crate::{
    BondOrder, BondQueryPredicate, DerivedState, Molecule, QueryGraph, QueryNode, SmartsParseParams,
    SubstructMatchParams, mol_from_smarts, try_get_substruct_matches_with_params,
};

fn terminal_atom_query() -> Result<&'static QueryGraph, &'static str> {
    static QUERY: OnceLock<Result<QueryGraph, &'static str>> = OnceLock::new();
    QUERY
        .get_or_init(|| {
            let params = SmartsParseParams {
                replacements: BTreeMap::from([("{atomPattern}".to_string(), "O,N;D1".to_string())]),
                ..SmartsParseParams::default()
            };
            mol_from_smarts(
                "[{atomPattern};$([{atomPattern}]-[*]=[{atomPattern}]),$([{atomPattern}]=[*]-[{atomPattern}])]~[*]",
                &params,
            )
            .map_err(|_| "bad terminal-group query pattern")
        })
        .as_ref()
        .map_err(|message| *message)
}

pub(crate) fn symmetrize_terminal_atoms(molecule: &Molecule) -> Result<Molecule, &'static str> {
    // BEGIN RDKIT CPP FUNCTION MolAlign::details::symmetrizeTerminalAtoms
    // RDKit✔️✔️: void symmetrizeTerminalAtoms(RWMol &mol) {
    // RDKit✔️✔️:   static const std::string qsmarts =
    // RDKit✔️✔️:       "[{atomPattern};$([{atomPattern}]-[*]=[{atomPattern}]),$([{atomPattern}]=[*]-[{atomPattern}])]~[*]";
    // RDKit✔️✔️:   static std::map<std::string, std::string> replacements = {
    // RDKit✔️✔️:       {"{atomPattern}", "O,N;D1"}};
    // RDKit✔️✔️:   auto matches = SubstructMatch(mol, *qry);
    // RDKit✔️✔️:   QueryBond qb;
    // RDKit✔️✔️:   qb.setQuery(makeSingleOrDoubleBondQuery());
    // RDKit✔️✔️:   for (const auto &match : matches) {
    // RDKit✔️✔️:     mol.getAtomWithIdx(match[0].second)->setFormalCharge(0);
    // RDKit✔️✔️:     auto obond = mol.getBondBetweenAtoms(match[0].second, match[1].second);
    // RDKit✔️✔️:     mol.replaceBond(obond->getIdx(), &qb);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION MolAlign::details::symmetrizeTerminalAtoms
    let query = terminal_atom_query()?;
    let matches = try_get_substruct_matches_with_params(molecule, query, &SubstructMatchParams::default())
        .map_err(|_| "terminal-group query matching is unsupported")?;
    let matched_atoms_and_bonds: Vec<_> = matches
        .into_iter()
        .map(|matched| {
            let terminal = matched.atom_mapping[0];
            let neighbor = matched.atom_mapping[1];
            let bond = molecule
                .topology_block()
                .adjacency
                .neighbors_of(terminal)
                .iter()
                .find(|entry| entry.atom_index == neighbor)
                .map(|entry| entry.bond.index())
                .ok_or("could not find expected terminal bond")?;
            Ok((terminal, bond))
        })
        .collect::<Result<_, &'static str>>()?;
    let mut symmetrized = molecule.clone();
    if matched_atoms_and_bonds.is_empty() {
        return Ok(symmetrized);
    }
    {
        let topology = symmetrized.topology_block_mut();
        for (terminal, bond) in matched_atoms_and_bonds {
            topology.atoms[terminal].set_formal_charge(0);
            topology.bonds[bond].set_query(Some(QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
                BondOrder::Single,
                BondOrder::Double,
            ]))));
        }
    }
    symmetrized
        .derived_cache_mut()
        .invalidate(DerivedState::VALENCE | DerivedState::AROMATICITY | DerivedState::STEREO);
    Ok(symmetrized)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn terminal_query_is_compiled_once_and_shared() {
        let first = terminal_atom_query().expect("terminal query");
        let second = terminal_atom_query().expect("cached terminal query");
        assert!(std::ptr::eq(first, second));
    }

    #[test]
    fn terminal_symmetrization_is_clone_only_and_invalidates_clone_caches() {
        let molecule = Molecule::from_smiles("[O-]C=O").expect("formate");
        assert!(molecule.derived_cache().valence.is_some());
        let original_atoms = molecule.atoms().to_vec();
        let original_bonds = molecule.bonds().to_vec();

        let symmetrized = symmetrize_terminal_atoms(&molecule).expect("symmetrized query clone");

        assert_eq!(molecule.atoms(), original_atoms);
        assert_eq!(molecule.bonds(), original_bonds);
        assert!(molecule.derived_cache().valence.is_some());
        assert!(symmetrized.derived_cache().valence.is_none());
        assert!(symmetrized.bonds().iter().all(|bond| bond.query().is_some()));
    }
}
