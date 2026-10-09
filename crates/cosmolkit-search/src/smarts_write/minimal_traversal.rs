//! Existing standalone SEARCH traversal, without the optional SMILES owner.
//! Token serialization remains shared with the full source writer.

use super::*;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn standalone_public_writer_keeps_whole_recursive_and_fragment_routes() {
        // RDKit MolToSmarts(MolFromSmarts(...)) and MolFragmentToSmarts(...,
        // [0, 1]); fixed source query trees, not inferred chemical equivalence.
        for (input, expected) in [
            ("[#6]-[#8]", "[#6]-[#8]"),
            ("[#6;$([#6]-[#8])]", "[#6&$([#6]-[#8])]"),
        ] {
            let graph = crate::parse_smarts(input, &Default::default()).unwrap();
            assert_eq!(
                query_graph_to_smarts(&graph, &Default::default()).unwrap(),
                PropertyText::from(expected)
            );
        }
        let graph = crate::parse_smarts("[#6]-[#8]-[#7]", &Default::default()).unwrap();
        assert_eq!(
            query_graph_fragment_to_smarts(
                &graph,
                &Default::default(),
                &[AtomId::new(0), AtomId::new(1)],
                None
            )
            .unwrap(),
            PropertyText::from("[#6]-[#8]")
        );
    }
}

pub(super) fn write(
    query: &QueryGraph,
    params: &SmartsWriteParams,
    atom_selection: Option<&[AtomId]>,
    bond_selection: Option<&[BondId]>,
    include_cx: bool,
) -> Result<SmartsWriteResult, SmartsWriteError> {
    #[cfg(not(feature = "smiles-integration"))]
    assert!(
        !include_cx,
        "CX entrypoints require the SMILES integration feature"
    );
    // RDKit✔️✔️:   PRECONDITION(!atomsToUse.empty(), "no atoms provided");
    // RDKit✔️✔️:   PRECONDITION(!bondsToUse || !bondsToUse->empty(), "no bonds provided");
    if atom_selection.is_some_and(<[AtomId]>::is_empty) {
        return Err(SmartsWriteError::EmptyAtomSelection);
    }
    if bond_selection.is_some_and(<[BondId]>::is_empty) {
        return Err(SmartsWriteError::EmptyBondSelection);
    }
    // RDKit✔️✔️:   SmilesWriteParams ps(params);
    // RDKit✔️✔️:   ps.rootedAtAtom = -1;
    // RDKit✔️✔️:   return molToSmarts(mol, ps, std::move(colors), atomsInPlay,
    // RDKit✔️✔️:                      bondsInPlay.get());
    let mut effective_params = *params;
    if atom_selection.is_some() {
        effective_params.rooted_at_atom = None;
    }
    // RDKit✔️✔️:   SmilesWriteParams ps(params);
    // RDKit✔️✔️:   ps.includeDativeBonds = false;
    // RDKit✔️✔️:   auto res = MolToSmarts(mol, ps);
    // The fragment wrapper calls MolFragmentToSmarts with the caller's
    // parameters unchanged, so this source override applies only to the
    // whole-graph CXSMARTS path.
    if include_cx && atom_selection.is_none() {
        effective_params.include_dative_bonds = false;
    }
    let atoms = atom_selection.map_or_else(
        || (0..query.num_atoms()).map(AtomId::new).collect::<Vec<_>>(),
        <[AtomId]>::to_vec,
    );
    if atoms.is_empty() || query.num_atoms() == 0 {
        return Ok(SmartsWriteResult::default());
    }
    for atom in &atoms {
        if atom.index() >= query.num_atoms() {
            return Err(SmartsWriteError::FragmentAtomOutOfRange { atom: atom.index() });
        }
    }
    if let Some(root) = effective_params.rooted_at_atom {
        if root >= query.num_atoms() || !atoms.iter().any(|atom| atom.index() == root) {
            return Err(SmartsWriteError::RootedAtomOutOfRange { atom: root });
        }
    }
    let selected = atoms
        .iter()
        .map(|atom| atom.index())
        .collect::<BTreeSet<_>>();
    let allowed_bonds = bond_selection.map_or_else(
        || {
            query
                .bonds()
                .iter()
                .filter(|bond| {
                    selected.contains(&bond.begin().index())
                        && selected.contains(&bond.end().index())
                })
                .map(|bond| bond.id().index())
                .collect::<BTreeSet<_>>()
        },
        |bonds| {
            bonds
                .iter()
                .map(|bond| bond.index())
                .collect::<BTreeSet<_>>()
        },
    );
    for bond in &allowed_bonds {
        if *bond >= query.num_bonds() {
            return Err(SmartsWriteError::FragmentBondOutOfRange { bond: *bond });
        }
    }
    validate_writer_carrier_valence_lists(query)?;
    let mut visited = vec![false; query.num_atoms()];
    let mut seen_bonds = BTreeSet::new();
    let mut tree_children = vec![Vec::<(BondId, AtomId)>::new(); query.num_atoms()];
    let mut ring_edges = Vec::<(BondId, AtomId, AtomId, usize)>::new();
    let mut next_ring = 1usize;
    let mut starts = atoms;
    // RDKit✔️✔️: if (params.rootedAtAtom > -1 &&
    // RDKit✔️✔️:     colors[params.rootedAtAtom] == Canon::WHITE_NODE) {
    // RDKit✔️✔️:   nextAtomIdx = params.rootedAtAtom;
    // RDKit✔️✔️: } else {
    // RDKit✔️✔️:   // Try to find a non-chiral atom we have not processed yet.
    // RDKit✔️✔️:   // If we can't find non-chiral atom, use the chiral atom with
    // RDKit✔️✔️:   // the lowest rank (we are guaranteed to find an unprocessed atom).
    // RDKit✔️✔️:   unsigned nextRank = nAtoms + 1;
    // RDKit✔️✔️:   for (auto atom : mol.atoms()) {
    // RDKit✔️✔️:     if (colors[atom->getIdx()] == Canon::WHITE_NODE) {
    // RDKit✔️✔️:       if (atom->getChiralTag() != Atom::CHI_TETRAHEDRAL_CCW &&
    // RDKit✔️✔️:           atom->getChiralTag() != Atom::CHI_TETRAHEDRAL_CW) {
    // RDKit✔️✔️:         nextAtomIdx = atom->getIdx();
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (ranks[atom->getIdx()] < nextRank) {
    // RDKit✔️✔️:         nextRank = ranks[atom->getIdx()];
    // RDKit✔️✔️:         nextAtomIdx = atom->getIdx();
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // Complexity review: the detached ordering key performs O(V log V) work
    // once, comparable to the source's repeated O(V) scans across components.
    // The DFS still visits every selected atom and allowed bond once.
    starts.sort_by_key(|atom| {
        let is_root = effective_params.rooted_at_atom == Some(atom.index());
        let is_tetrahedral = query.atom(atom.index()).is_some_and(|query_atom| {
            matches!(
                query_atom.chiral_tag(),
                ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
            )
        });
        (
            usize::from(!is_root),
            usize::from(!is_root && is_tetrahedral),
            atom.index(),
        )
    });
    for start in &starts {
        if visited[start.index()] {
            continue;
        }
        classify_query_graph(
            query,
            *start,
            None,
            &selected,
            &allowed_bonds,
            &mut visited,
            &mut seen_bonds,
            &mut tree_children,
            &mut ring_edges,
            &mut next_ring,
        )?;
    }
    renumber_ring_edges(&mut ring_edges, &tree_children, &starts, query.num_atoms());

    let mut result = SmartsWriteResult::default();
    visited.fill(false);
    let mut component_count = 0usize;
    for start in &starts {
        if visited[start.index()] {
            continue;
        }
        if component_count > 0 {
            result.smarts.push_byte(b'.');
        }
        component_count += 1;
        emit_query_graph(
            query,
            *start,
            &tree_children,
            &ring_edges,
            &mut visited,
            &effective_params,
            &mut result,
        )?;
    }
    // RDKit❗✔️:   inmol.setProp(common_properties::_smilesAtomOutputOrder, atomOrdering, true);
    // RDKit❗✔️:   inmol.setProp(common_properties::_smilesBondOutputOrder, bondOrdering, true);
    result.source_orders_written = true;
    Ok(result)
}

fn validate_writer_carrier_valence_lists(query: &QueryGraph) -> Result<(), SmartsWriteError> {
    // RDKit❗✔️:   for (auto &atom : mol.atoms()) {
    // RDKit❗✔️:     atom->updatePropertyCache(false);
    // RDKit❗✔️:   }
    // RDKit❗✔️: void Atom::updatePropertyCache(bool strict) {
    // RDKit❗✔️:   calcExplicitValence(strict);
    // RDKit❗✔️:   calcImplicitValence(strict);
    // RDKit❗✔️: }
    // RDKit✔️✔️:   const auto &ovalens =
    // RDKit✔️✔️:       PeriodicTable::getTable()->getValenceList(atom.getAtomicNum());
    // RDKit✔️✔️:   const INT_VECT &getValenceList(UINT atomicNumber) const {
    // RDKit✔️✔️:     PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    // RDKit✔️✔️:     return byanum[atomicNumber].ValenceList();
    // RDKit✔️✔️:   }
    // Source boundary: FragmentSmartsConstruct updates all carriers, including
    // atoms outside a fragment selection. calculateExplicitValence requests
    // this list unconditionally, before the atom/query writer runs. Reuse the
    // foundational table owner; a numeric query leaf is not a carrier identity.
    // This helper reproduces that lookup/precondition only; it does not assign
    // valence caches or claim to reproduce the complete property-cache update.
    // Complexity: one indexed lookup per carrier, O(V), no allocation. Repeating
    // the same immutable lookups for later components cannot change the result.
    for atom in query.atoms() {
        cosmolkit_core::required_valence_list(atom.atomic_number())?;
    }
    Ok(())
}

fn renumber_ring_edges(
    ring_edges: &mut [(BondId, AtomId, AtomId, usize)],
    tree_children: &[Vec<(BondId, AtomId)>],
    starts: &[AtomId],
    atom_count: usize,
) {
    // RDKit✔️✔️:     for (auto bIdx : atomRingClosures[atomIdx]) {
    // RDKit✔️✔️:       unsigned int ringIdx = std::numeric_limits<unsigned int>::max();
    // RDKit✔️✔️:       if (bond->getPropIfPresent(common_properties::_TraversalRingClosureBond,
    // RDKit✔️✔️:                                  ringIdx)) {
    // RDKit✔️✔️:         // this is end of the ring closure
    // RDKit✔️✔️:         // we can just pull the ring index from the bond itself:
    // RDKit✔️✔️:         molStack.push_back(MolStackElem(bond, atomIdx));
    // RDKit✔️✔️:         molStack.push_back(MolStackElem(ringIdx));
    // RDKit✔️✔️:         // don't make the ring digit immediately available again: we don't want
    // RDKit✔️✔️:         // to have the same
    // RDKit✔️✔️:         // ring digit opening and closing rings on an atom.
    // RDKit✔️✔️:         ringsClosed.push_back(ringIdx - 1);
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         auto lowestRingIdx = cyclesAvailable.find_first();
    // RDKit✔️✔️:         cyclesAvailable.set(lowestRingIdx, false);
    // RDKit✔️✔️:         ++lowestRingIdx;
    // RDKit✔️✔️:         bond->setProp(common_properties::_TraversalRingClosureBond,
    // RDKit✔️✔️:                       static_cast<unsigned int>(lowestRingIdx));
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     for (auto ringIdx : ringsClosed) {
    // RDKit✔️✔️:       cyclesAvailable.set(ringIdx);
    // RDKit✔️✔️:     }
    if ring_edges.is_empty() {
        return;
    }

    let mut visited = vec![false; atom_count];
    let mut atom_order = Vec::with_capacity(atom_count);
    for start in starts {
        collect_query_atom_order(*start, tree_children, &mut visited, &mut atom_order);
    }
    let mut positions = vec![usize::MAX; atom_count];
    for (position, atom) in atom_order.into_iter().enumerate() {
        positions[atom.index()] = position;
    }

    let mut occurrences = Vec::with_capacity(ring_edges.len() * 2);
    for (edge_index, (_, first, second, _)) in ring_edges.iter().enumerate() {
        occurrences.push((positions[first.index()], edge_index));
        occurrences.push((positions[second.index()], edge_index));
    }
    occurrences.sort_unstable();

    let mut labels = vec![0usize; ring_edges.len()];
    let mut available = BTreeSet::new();
    let mut next_label = 1usize;
    let mut occurrence_index = 0usize;
    while occurrence_index < occurrences.len() {
        let atom_position = occurrences[occurrence_index].0;
        let mut group_end = occurrence_index;
        while group_end < occurrences.len() && occurrences[group_end].0 == atom_position {
            group_end += 1;
        }

        let mut closed_at_atom = Vec::new();
        for (_, edge_index) in &occurrences[occurrence_index..group_end] {
            if labels[*edge_index] == 0 {
                labels[*edge_index] = available.pop_first().unwrap_or_else(|| {
                    let label = next_label;
                    next_label += 1;
                    label
                });
            } else {
                closed_at_atom.push(labels[*edge_index]);
            }
        }
        for label in closed_at_atom {
            available.insert(label);
        }
        occurrence_index = group_end;
    }
    for (edge, label) in ring_edges.iter_mut().zip(labels) {
        edge.3 = label;
    }
}

fn collect_query_atom_order(
    atom: AtomId,
    tree_children: &[Vec<(BondId, AtomId)>],
    visited: &mut [bool],
    order: &mut Vec<AtomId>,
) {
    if visited[atom.index()] {
        return;
    }
    visited[atom.index()] = true;
    order.push(atom);
    for (_, child) in &tree_children[atom.index()] {
        collect_query_atom_order(*child, tree_children, visited, order);
    }
}

fn classify_query_graph(
    query: &QueryGraph,
    atom: AtomId,
    parent_bond: Option<BondId>,
    selected: &BTreeSet<usize>,
    allowed_bonds: &BTreeSet<usize>,
    visited: &mut [bool],
    seen_bonds: &mut BTreeSet<usize>,
    tree_children: &mut [Vec<(BondId, AtomId)>],
    ring_edges: &mut Vec<(BondId, AtomId, AtomId, usize)>,
    next_ring: &mut usize,
) -> Result<(), SmartsWriteError> {
    visited[atom.index()] = true;
    let mut incident = query
        .adjacency()
        .get(atom.index())
        .into_iter()
        .flatten()
        .filter_map(|(other, bond)| {
            (selected.contains(other) && allowed_bonds.contains(bond))
                .then_some((BondId::new(*bond), AtomId::new(*other)))
        })
        .collect::<Vec<_>>();
    incident.sort_by_key(|(bond, other)| (other.index(), bond.index()));
    for (bond, other) in incident {
        if Some(bond) == parent_bond || !seen_bonds.insert(bond.index()) {
            continue;
        }
        if visited[other.index()] {
            ring_edges.push((bond, atom, other, *next_ring));
            *next_ring += 1;
        } else {
            tree_children[atom.index()].push((bond, other));
            classify_query_graph(
                query,
                other,
                Some(bond),
                selected,
                allowed_bonds,
                visited,
                seen_bonds,
                tree_children,
                ring_edges,
                next_ring,
            )?;
        }
    }
    Ok(())
}

fn emit_query_graph(
    query: &QueryGraph,
    atom: AtomId,
    tree_children: &[Vec<(BondId, AtomId)>],
    ring_edges: &[(BondId, AtomId, AtomId, usize)],
    visited: &mut [bool],
    params: &SmartsWriteParams,
    result: &mut SmartsWriteResult,
) -> Result<(), SmartsWriteError> {
    // RDKit✔️✔️:       case Canon::MOL_STACK_ATOM: {
    // RDKit✔️✔️:         auto *atm = msCI.obj.atom;
    // RDKit✔️✔️:         res << SmartsWrite::GetAtomSmarts(atm, params);
    // RDKit✔️✔️:         atomOrdering.push_back(atm->getIdx());
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       case Canon::MOL_STACK_BOND: {
    // RDKit✔️✔️:         auto *bnd = msCI.obj.bond;
    // RDKit✔️✔️:         res << SmartsWrite::GetBondSmarts(bnd, params, msCI.number);
    // RDKit✔️✔️:         bondOrdering.push_back(bnd->getIdx());
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       case Canon::MOL_STACK_BRANCH_OPEN: {
    // RDKit✔️✔️:         res << "(";
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       case Canon::MOL_STACK_BRANCH_CLOSE: {
    // RDKit✔️✔️:         res << ")";
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    visited[atom.index()] = true;
    let query_atom = query
        .atom(atom.index())
        .ok_or(SmartsWriteError::FragmentAtomOutOfRange { atom: atom.index() })?;
    let mut atom_params = *params;
    if query.prop("_doIsoSmiles").is_some() {
        atom_params.isomeric_smiles = true;
    }
    result
        .smarts
        .extend_bytes(query_atom_to_smarts(query_atom, &atom_params)?.as_bytes());
    result.atom_ordering.push(atom);
    for (bond, first, _second, ring_number) in ring_edges
        .iter()
        .filter(|(_, first, second, _)| *first == atom || *second == atom)
    {
        if *first == atom {
            result.smarts.extend_bytes(
                (&query_bond_to_smarts(
                    query
                        .bond(bond.index())
                        .ok_or(SmartsWriteError::FragmentBondOutOfRange { bond: bond.index() })?,
                    params,
                    Some(atom.index()),
                )?)
                    .as_ref(),
            );
            result.bond_ordering.push(*bond);
        }
        if *ring_number < 10 {
            result
                .smarts
                .extend_bytes((&ring_number.to_string()).as_ref());
        } else {
            result.smarts.push_byte(b'%');
            result
                .smarts
                .extend_bytes((&ring_number.to_string()).as_ref());
        }
    }
    let children = &tree_children[atom.index()];
    for (index, (bond, other)) in children.iter().enumerate() {
        if index + 1 != children.len() {
            result.smarts.push_byte(b'(');
        }
        result.smarts.extend_bytes(
            (&query_bond_to_smarts(
                query
                    .bond(bond.index())
                    .ok_or(SmartsWriteError::FragmentBondOutOfRange { bond: bond.index() })?,
                params,
                Some(atom.index()),
            )?)
                .as_ref(),
        );
        result.bond_ordering.push(*bond);
        emit_query_graph(
            query,
            *other,
            tree_children,
            ring_edges,
            visited,
            params,
            result,
        )?;
        if index + 1 != children.len() {
            result.smarts.push_byte(b')');
        }
    }
    Ok(())
}
