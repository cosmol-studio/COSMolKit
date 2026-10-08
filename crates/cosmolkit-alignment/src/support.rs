//! Private source-defined temporary query adaptation. No live state is changed.
use crate::{AlignmentError, AlignmentInput};
use cosmolkit_model::{
    AtomQueryPredicate, Bond, BondQueryPredicate, BondSpec, QueryAtom, QueryBond, QueryGraph,
    QueryNode,
};
use cosmolkit_search::{
    QueryMatchContext, SearchTarget, SearchTargetAccess, SmartsParseError, SmartsParseParams,
    SubstructMatchParams,
};
use cosmolkit_types::BondOrder;
use std::{collections::BTreeMap, sync::OnceLock};

pub(super) fn query_for(input: &AlignmentInput<'_>) -> Result<QueryGraph, AlignmentError> {
    // Ordinary source ROMol query rows retain their non-query origin. SEARCH
    // dispatches its existing Atom::Match/Bond::Match carrier comparisons.
    // No perception, SMARTS serialization or chemistry reconstruction occurs.
    // Performance (behavior reproduced, extra allocation): the source ordinary
    // matcher borrows ROMol; the detached QueryGraph
    // below clones O(V+E) atoms/bonds and stereo groups, an extra allocation cost.
    Ok(QueryGraph::from_parts(
        input
            .topology
            .atoms
            .iter()
            .cloned()
            .map(|atom| {
                QueryAtom::from_carrier_parts(atom, QueryNode::predicate(AtomQueryPredicate::Any))
            })
            .collect(),
        input
            .topology
            .bonds
            .iter()
            .cloned()
            .map(|bond| {
                QueryBond::from_carrier_parts(bond, QueryNode::predicate(BondQueryPredicate::Any))
            })
            .collect(),
        BTreeMap::new(),
        Vec::new(),
        Vec::new(),
        input.topology.stereo_groups.clone(),
    )?)
}

fn terminal_atom_query() -> Result<&'static QueryGraph, AlignmentError> {
    // Pinned RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8: Code/GraphMol/MolAlign/AlignMolecules.cpp
    // BEGIN VERBATIM CPP symmetrizeTerminalAtoms
    // RDKit✔️✔️: void symmetrizeTerminalAtoms(RWMol &mol) {
    // RDKit✔️✔️:   // clang-format off
    // RDKit✔️✔️:   static const std::string qsmarts =
    // RDKit✔️✔️:       "[{atomPattern};$([{atomPattern}]-[*]=[{atomPattern}]),$([{atomPattern}]=[*]-[{atomPattern}])]~[*]";
    // RDKit✔️✔️:   static std::map<std::string, std::string> replacements = {
    // RDKit✔️✔️:       {"{atomPattern}", "O,N;D1"}};
    // RDKit✔️✔️:   // clang-format on
    // RDKit✔️✔️:   static SmartsParserParams ps;
    // RDKit✔️✔️:   ps.replacements = &replacements;
    // RDKit✔️✔️:   static const std::unique_ptr<RWMol> qry{SmartsToMol(qsmarts, ps)};
    // RDKit✔️✔️:   CHECK_INVARIANT(qry, "bad query pattern");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto matches = SubstructMatch(mol, *qry);
    // RDKit✔️✔️:   if (matches.empty()) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   QueryBond qb;
    // RDKit✔️✔️:   qb.setQuery(makeSingleOrDoubleBondQuery());
    // RDKit✔️✔️:   for (const auto &match : matches) {
    // RDKit✔️✔️:     mol.getAtomWithIdx(match[0].second)->setFormalCharge(0);
    // RDKit✔️✔️:     auto obond = mol.getBondBetweenAtoms(match[0].second, match[1].second);
    // RDKit✔️✔️:     CHECK_INVARIANT(obond, "could not find expected bond");
    // RDKit✔️✔️:     mol.replaceBond(obond->getIdx(), &qb);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END VERBATIM CPP symmetrizeTerminalAtoms

    static QUERY: OnceLock<Result<QueryGraph, SmartsParseError>> = OnceLock::new();
    QUERY.get_or_init(|| {
        let params = SmartsParseParams {
            replacements: BTreeMap::from([("{atomPattern}".into(), "O,N;D1".into())]),
            ..Default::default()
        };
        cosmolkit_search::parse_smarts(
            "[{atomPattern};$([{atomPattern}]-[*]=[{atomPattern}]),$([{atomPattern}]=[*]-[{atomPattern}])]~[*]", &params,
        )
    }).as_ref().map_err(|cause| AlignmentError::QueryParse(cause.clone()))
}

pub(super) fn symmetrize_terminal_atoms(
    input: &AlignmentInput<'_>,
) -> Result<QueryGraph, AlignmentError> {
    symmetrize_terminal_query_impl(query_for(input)?, &input.search_target(), None)
}

/// Source terminal rewrite for detached conformer pruning with explicit match context.
pub fn symmetrize_terminal_query_with_context(
    query: QueryGraph,
    target: &SearchTarget<'_>,
    context: &QueryMatchContext,
) -> Result<QueryGraph, AlignmentError> {
    symmetrize_terminal_query_impl(query, target, Some(context))
}

fn symmetrize_terminal_query_impl(
    mut symmetrized: QueryGraph,
    target: &SearchTarget<'_>,
    context: Option<&QueryMatchContext>,
) -> Result<QueryGraph, AlignmentError> {
    // Pinned RDKit 351f8f378f8ad6bbd517980c38896e66bf907af8: Code/GraphMol/MolAlign/AlignMolecules.cpp
    // BEGIN VERBATIM CPP symmetrizeTerminalAtoms
    // RDKit✔️✔️: void symmetrizeTerminalAtoms(RWMol &mol) {
    // RDKit✔️✔️:   // clang-format off
    // RDKit✔️✔️:   static const std::string qsmarts =
    // RDKit✔️✔️:       "[{atomPattern};$([{atomPattern}]-[*]=[{atomPattern}]),$([{atomPattern}]=[*]-[{atomPattern}])]~[*]";
    // RDKit✔️✔️:   static std::map<std::string, std::string> replacements = {
    // RDKit✔️✔️:       {"{atomPattern}", "O,N;D1"}};
    // RDKit✔️✔️:   // clang-format on
    // RDKit✔️✔️:   static SmartsParserParams ps;
    // RDKit✔️✔️:   ps.replacements = &replacements;
    // RDKit✔️✔️:   static const std::unique_ptr<RWMol> qry{SmartsToMol(qsmarts, ps)};
    // RDKit✔️✔️:   CHECK_INVARIANT(qry, "bad query pattern");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto matches = SubstructMatch(mol, *qry);
    // RDKit✔️✔️:   if (matches.empty()) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   QueryBond qb;
    // RDKit✔️✔️:   qb.setQuery(makeSingleOrDoubleBondQuery());
    // RDKit✔️✔️:   for (const auto &match : matches) {
    // RDKit✔️✔️:     mol.getAtomWithIdx(match[0].second)->setFormalCharge(0);
    // RDKit✔️✔️:     auto obond = mol.getBondBetweenAtoms(match[0].second, match[1].second);
    // RDKit✔️✔️:     CHECK_INVARIANT(obond, "could not find expected bond");
    // RDKit✔️✔️:     mol.replaceBond(obond->getIdx(), &qb);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END VERBATIM CPP symmetrizeTerminalAtoms

    let query = terminal_atom_query()?;
    let matches = if let Some(context) = context {
        cosmolkit_search::try_get_substruct_matches_with_params_and_context(
            target,
            query,
            &SubstructMatchParams::default(),
            context,
        )?
    } else {
        cosmolkit_search::try_get_substruct_matches_with_params(
            target,
            query,
            &SubstructMatchParams::default(),
        )?
    };
    for matched in matches {
        let terminal = matched.atom_mapping[0];
        let neighbor = matched.atom_mapping[1];
        let bond = target
            .topology_block()
            .adjacency
            .neighbors_of(terminal)
            .iter()
            .find(|entry| entry.atom_index == neighbor)
            .map(|entry| entry.bond.index())
            .ok_or(AlignmentError::TerminalGroupSymmetrization {
                message: "could not find expected terminal bond",
            })?;
        symmetrized.atoms_mut()[terminal].set_formal_charge(0);
        let original = &target.topology_block().bonds[bond];
        // RWMol::replaceBond supplies index and endpoints to the default QueryBond;
        // the remaining carrier state is default, not the previous bond's state.
        let carrier = Bond::from_spec(
            original.id(),
            BondSpec::new(original.begin(), original.end(), BondOrder::Unspecified),
        );
        symmetrized.bonds_mut()[bond] = QueryBond::from_parts(
            carrier,
            QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
                BondOrder::Single,
                BondOrder::Double,
            ])),
        );
    }
    Ok(symmetrized)
}

#[cfg(test)]
mod tests {
    use super::*;
    use cosmolkit_model::{Atom, AtomSpec, TopologyBlock};
    use cosmolkit_model::{AtomId, BondId};
    use cosmolkit_types::Element;
    #[test]
    fn terminal_query_is_compiled_once_and_shared() {
        let first = terminal_atom_query().expect("terminal query");
        let second = terminal_atom_query().expect("cached terminal query");
        assert!(std::ptr::eq(first, second));
    }
    #[test]
    fn terminal_symmetrization_is_clone_only_and_invalidates_clone_caches() {
        // The domain's detached query has no derived-cache storage. It cannot
        // retain a stale cloned valence/aromaticity/stereo cache; the borrowed
        // source topology and its externally owned derived state are untouched.
        let mut atoms = vec![
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::O)),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(2), AtomSpec::new(Element::O)),
        ];
        atoms[0].set_formal_charge(-1);
        let bonds = vec![
            Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            ),
            Bond::from_spec(
                BondId::new(1),
                BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Double),
            ),
        ];
        let topology = TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap();
        let original = topology.clone();
        let coordinates = Default::default();
        let input = AlignmentInput {
            topology: &topology,
            coordinates: &coordinates,
            rings: None,
            valence: None,
        };
        let symmetrized = symmetrize_terminal_atoms(&input).expect("symmetrized query clone");
        assert_eq!(topology, original);
        assert_eq!(symmetrized.atoms()[0].formal_charge(), 0);
        assert_eq!(symmetrized.atoms()[2].formal_charge(), 0);
        assert!(
            symmetrized
                .bonds()
                .iter()
                .all(|b| !b.predicate_is_carrier_derived())
        );
        assert!(
            symmetrized
                .bonds()
                .iter()
                .all(|b| b.bond().order() == BondOrder::Unspecified)
        );
    }
}
