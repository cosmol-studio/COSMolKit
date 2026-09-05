//! RDKit SMARTS serialization for canonical [`QueryGraph`] values.

use super::query::{AtomRangeBounds, AtomRangeDataFunction, QueryNode, RecursiveStructureQuery};
use crate::{
    Atom, AtomId, AtomQueryPredicate, Bond, BondDirection, BondId, BondOrder, BondQueryPredicate,
    ChiralTag, Element, Hybridization, Molecule, QueryAtom, QueryBond, QueryGraph, RingFindType,
    RingInfo, SmilesWriteParams,
};
use std::collections::BTreeSet;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
struct QueryBoolFeatures(u8);

impl QueryBoolFeatures {
    const HAS_AND: Self = Self(0x1);
    const HAS_LOW_AND: Self = Self(0x2);
    const HAS_OR: Self = Self(0x4);
    const HAS_RECURSION: Self = Self(0x8);

    const fn contains(self, other: Self) -> bool {
        self.0 & other.0 != 0
    }

    fn insert(&mut self, other: Self) {
        self.0 |= other.0;
    }
}

impl std::ops::BitOrAssign for QueryBoolFeatures {
    fn bitor_assign(&mut self, rhs: Self) {
        self.0 |= rhs.0;
    }
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum SmartsWriteError {
    #[error("SMARTS writer query-graph traversal is not available for this graph: {detail}")]
    QueryGraphTraversalUnsupported { detail: &'static str },
    #[error("This is a non-smartable query - OR above and below AND in the binary tree")]
    OrAboveAndBelowAnd,
    #[error("Don't know how to combine using {description}")]
    UnknownCombination { description: String },
    #[error("recursive SMARTS query has no query molecule")]
    MissingRecursiveQueryMolecule,
    #[error("Can't write smarts for this bond dir type: {direction:?}")]
    UnsupportedBondDirection { direction: BondDirection },
    #[error("Can't write smarts for this query bond type: {predicate:?}")]
    UnsupportedBondQuery { predicate: BondQueryPredicate },
    #[error("SMARTS {kind} composite query requires at least two children")]
    CompositeChildCount { kind: &'static str },
    #[error("SMARTS writer does not support XOR query composites")]
    XorComposite,
    #[error("CXSMARTS extensions for QueryGraph are not yet supported: {detail}")]
    QueryGraphCxExtensionsUnsupported { detail: &'static str },
    #[error("SMARTS traversal failed: {source}")]
    Traversal {
        #[source]
        source: crate::SmilesWriteError,
    },
    #[error("rooted atom index {atom} is out of range")]
    RootedAtomOutOfRange { atom: usize },
    #[error("SMARTS fragment requires at least one atom")]
    EmptyAtomSelection,
    #[error("an explicit SMARTS fragment bond selection cannot be empty")]
    EmptyBondSelection,
    #[error("SMARTS fragment atom index {atom} is out of range")]
    FragmentAtomOutOfRange { atom: usize },
    #[error("SMARTS fragment bond index {bond} is out of range")]
    FragmentBondOutOfRange { bond: usize },
    #[error("atom {atom} is not an endpoint of bond {bond}")]
    BondAtomNotEndpoint { bond: usize, atom: usize },
}

/// Serialize an independent query graph.
///
/// The canonical molecule traversal is intentionally not used here: a query
/// graph is not a `Molecule` and must never be lowered back into one.
pub(crate) fn query_graph_to_smarts(
    query: &QueryGraph,
    params: &SmilesWriteParams,
) -> Result<String, SmartsWriteError> {
    query_graph_to_smarts_fragment(query, params, None, None, false)
}

pub(crate) fn query_graph_to_cx_smarts(
    query: &QueryGraph,
    params: &SmilesWriteParams,
) -> Result<String, SmartsWriteError> {
    query_graph_to_smarts_fragment(query, params, None, None, true)
}

pub(crate) fn query_graph_fragment_to_smarts(
    query: &QueryGraph,
    params: &SmilesWriteParams,
    atoms: &[AtomId],
    bonds: Option<&[BondId]>,
) -> Result<String, SmartsWriteError> {
    query_graph_to_smarts_fragment(query, params, Some(atoms), bonds, false)
}

pub(crate) fn query_graph_fragment_to_cx_smarts(
    query: &QueryGraph,
    params: &SmilesWriteParams,
    atoms: &[AtomId],
    bonds: Option<&[BondId]>,
) -> Result<String, SmartsWriteError> {
    query_graph_to_smarts_fragment(query, params, Some(atoms), bonds, true)
}

fn query_graph_to_smarts_fragment(
    query: &QueryGraph,
    params: &SmilesWriteParams,
    atom_selection: Option<&[AtomId]>,
    bond_selection: Option<&[BondId]>,
    include_cx: bool,
) -> Result<String, SmartsWriteError> {
    if include_cx {
        ensure_query_graph_cx_extensions_supported(query)?;
    }
    let atoms = atom_selection.map_or_else(
        || (0..query.num_atoms()).map(AtomId::new).collect::<Vec<_>>(),
        <[AtomId]>::to_vec,
    );
    if atoms.is_empty() || query.num_atoms() == 0 {
        return Ok(String::new());
    }
    for atom in &atoms {
        if atom.index() >= query.num_atoms() {
            return Err(SmartsWriteError::FragmentAtomOutOfRange { atom: atom.index() });
        }
    }
    if let Some(root) = params.rooted_at_atom {
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
    let mut visited = vec![false; query.num_atoms()];
    let mut seen_bonds = BTreeSet::new();
    let mut tree_children = vec![Vec::<(BondId, AtomId)>::new(); query.num_atoms()];
    let mut ring_edges = Vec::<(BondId, AtomId, AtomId, usize)>::new();
    let mut next_ring = 1usize;
    let mut starts = atoms;
    starts.sort_by_key(|atom| atom.index());
    if let Some(root) = params.rooted_at_atom {
        starts.sort_by_key(|atom| usize::from(atom.index() != root));
    }
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

    let mut output = String::new();
    visited.fill(false);
    let mut component_count = 0usize;
    for start in &starts {
        if visited[start.index()] {
            continue;
        }
        if component_count > 0 {
            output.push('.');
        }
        component_count += 1;
        emit_query_graph(
            query,
            *start,
            &tree_children,
            &ring_edges,
            &mut visited,
            params,
            &mut output,
        )?;
    }
    Ok(output)
}

fn renumber_ring_edges(
    ring_edges: &mut [(BondId, AtomId, AtomId, usize)],
    tree_children: &[Vec<(BondId, AtomId)>],
    starts: &[AtomId],
    atom_count: usize,
) {
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
    for (_, edge_index) in occurrences {
        if labels[edge_index] == 0 {
            labels[edge_index] = available.pop_first().unwrap_or_else(|| {
                let label = next_label;
                next_label += 1;
                label
            });
        } else {
            available.insert(labels[edge_index]);
        }
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

fn ensure_query_graph_cx_extensions_supported(query: &QueryGraph) -> Result<(), SmartsWriteError> {
    if query.coordinates_2d().is_some() || !query.conformers_3d().is_empty() {
        return Err(SmartsWriteError::QueryGraphCxExtensionsUnsupported {
            detail: "coordinate extensions",
        });
    }
    if !query.stereo_groups().is_empty() {
        return Err(SmartsWriteError::QueryGraphCxExtensionsUnsupported {
            detail: "enhanced stereo groups",
        });
    }
    if query.props().keys().any(|key| key != "_Name") {
        return Err(SmartsWriteError::QueryGraphCxExtensionsUnsupported {
            detail: "molecule properties",
        });
    }
    const CX_ATOM_PROPERTIES: &[&str] = &[
        "atomLabel",
        "_QueryAtomGenericLabel",
        "dummyLabel",
        "_fromAttachPoint",
        "molFileValue",
    ];
    if query.atoms().iter().any(|atom| {
        CX_ATOM_PROPERTIES
            .iter()
            .any(|key| atom.prop(key).is_some())
    }) {
        return Err(SmartsWriteError::QueryGraphCxExtensionsUnsupported {
            detail: "atom labels or values",
        });
    }
    if query.bonds().iter().any(|bond| {
        bond.bond()
            .props()
            .keys()
            .any(|key| key != crate::notation::smiles::UNSPECIFIED_ORDER_PROP)
    }) {
        return Err(SmartsWriteError::QueryGraphCxExtensionsUnsupported {
            detail: "bond extensions",
        });
    }
    Ok(())
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
    params: &SmilesWriteParams,
    output: &mut String,
) -> Result<(), SmartsWriteError> {
    visited[atom.index()] = true;
    let query_atom = query
        .atom(atom.index())
        .ok_or(SmartsWriteError::FragmentAtomOutOfRange { atom: atom.index() })?;
    output.push_str(&query_atom_to_smarts(query_atom, params)?);
    for (bond, first, _second, ring_number) in ring_edges
        .iter()
        .filter(|(_, first, second, _)| *first == atom || *second == atom)
    {
        if *first == atom {
            output.push_str(&query_bond_to_smarts(
                query
                    .bond(bond.index())
                    .ok_or(SmartsWriteError::FragmentBondOutOfRange { bond: bond.index() })?,
                params,
                Some(atom.index()),
            )?);
        }
        if *ring_number < 10 {
            output.push_str(&ring_number.to_string());
        } else {
            output.push('%');
            output.push_str(&ring_number.to_string());
        }
    }
    let children = &tree_children[atom.index()];
    for (index, (bond, other)) in children.iter().enumerate() {
        if index + 1 != children.len() {
            output.push('(');
        }
        output.push_str(&query_bond_to_smarts(
            query
                .bond(bond.index())
                .ok_or(SmartsWriteError::FragmentBondOutOfRange { bond: bond.index() })?,
            params,
            Some(atom.index()),
        )?);
        emit_query_graph(
            query,
            *other,
            tree_children,
            ring_edges,
            visited,
            params,
            output,
        )?;
        if index + 1 != children.len() {
            output.push(')');
        }
    }
    Ok(())
}

pub(crate) fn query_atom_to_smarts(
    atom: &QueryAtom,
    params: &SmilesWriteParams,
) -> Result<String, SmartsWriteError> {
    let mut features = QueryBoolFeatures::default();
    let mut stereo_written = false;
    let mut needs_brackets = false;
    let mut result = match atom.predicate() {
        QueryNode::Predicate(AtomQueryPredicate::RecursiveSmarts(query)) => {
            needs_brackets = true;
            get_recursive_structure_query_smarts(query, false, params, query_graph_to_smarts)?
        }
        QueryNode::And(_) | QueryNode::Or(_) => {
            needs_brackets = true;
            recurse_get_smarts(
                atom.atom(),
                atom.predicate(),
                false,
                &mut features,
                params,
                &mut stereo_written,
                &mut query_graph_to_smarts,
            )?
        }
        QueryNode::Xor(_) => return Err(SmartsWriteError::XorComposite),
        QueryNode::Predicate(predicate) => {
            let mut need_paren = false;
            let result = get_atom_smarts_simple(
                atom.atom(),
                predicate,
                &mut need_paren,
                true,
                params.do_isomeric_smiles,
                &mut stereo_written,
            );
            needs_brackets = need_paren;
            result
        }
        QueryNode::Not(child) => {
            needs_brackets = true;
            let mut need_paren = false;
            let mut recursive_negation_written = false;
            let mut value = match child.as_ref() {
                QueryNode::Predicate(AtomQueryPredicate::RecursiveSmarts(query)) => {
                    recursive_negation_written = true;
                    get_recursive_structure_query_smarts(
                        query,
                        true,
                        params,
                        query_graph_to_smarts,
                    )?
                }
                QueryNode::Predicate(predicate) => get_atom_smarts_simple(
                    atom.atom(),
                    predicate,
                    &mut need_paren,
                    true,
                    params.do_isomeric_smiles,
                    &mut stereo_written,
                ),
                _ => recurse_get_smarts(
                    atom.atom(),
                    child.as_ref(),
                    false,
                    &mut features,
                    params,
                    &mut stereo_written,
                    &mut query_graph_to_smarts,
                )?,
            };
            if !recursive_negation_written {
                value.insert(0, '!');
            }
            value
        }
    };
    if let Some(map) = atom.atom_map() {
        needs_brackets = true;
        result.push(':');
        result.push_str(&map.to_string());
    }
    if let Some(symbol) = atom.prop("smilesSymbol") {
        needs_brackets = true;
        result = format!("{symbol};{result}");
    }
    if needs_brackets {
        result = format!("[{result}]");
    }
    Ok(result)
}

pub(crate) fn query_bond_to_smarts(
    bond: &QueryBond,
    params: &SmilesWriteParams,
    atom_to_left_idx: Option<usize>,
) -> Result<String, SmartsWriteError> {
    let mut features = QueryBoolFeatures::default();
    match bond.predicate() {
        QueryNode::And(_) | QueryNode::Or(_) => recurse_bond_smarts(
            bond.bond(),
            bond.predicate(),
            false,
            atom_to_left_idx,
            &mut features,
            params,
        ),
        QueryNode::Xor(_) => Err(SmartsWriteError::XorComposite),
        QueryNode::Predicate(predicate) => {
            get_bond_smarts_simple(bond.bond(), predicate, atom_to_left_idx, params)
        }
        QueryNode::Not(child) => {
            let mut result = match child.as_ref() {
                QueryNode::Predicate(predicate) => {
                    get_bond_smarts_simple(bond.bond(), predicate, atom_to_left_idx, params)?
                }
                _ => {
                    return recurse_bond_smarts(
                        bond.bond(),
                        child.as_ref(),
                        false,
                        atom_to_left_idx,
                        &mut features,
                        params,
                    );
                }
            };
            result.insert(0, '!');
            Ok(result)
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
struct SmartsWriteResult {
    smarts: String,
    atom_ordering: Vec<AtomId>,
    bond_ordering: Vec<BondId>,
}

fn combine_child_smarts(
    child1: String,
    features1: QueryBoolFeatures,
    child2: String,
    features2: QueryBoolFeatures,
    description: &str,
    features: &mut QueryBoolFeatures,
) -> Result<String, SmartsWriteError> {
    // RDKit✔️✔️: std::string _combineChildSmarts(std::string cs1, unsigned int features1,
    // RDKit✔️✔️:                                 std::string cs2, unsigned int features2,
    // RDKit✔️✔️:                                 std::string descrip, unsigned int &features) {
    // RDKit✔️✔️:   std::string res = "";
    // RDKit✔️✔️:   if ((descrip.find("Or") > 0) && (descrip.find("Or") < descrip.length())) {
    // RDKit✔️✔️:     // if either of child smarts already have a "," and ";" we can't have one
    // RDKit✔️✔️:     // more OR here
    // RDKit✔️✔️:     if ((features1 & static_cast<unsigned int>(QueryBoolFeatures::HAS_LOWAND) &&
    // RDKit✔️✔️:          features1 & static_cast<unsigned int>(QueryBoolFeatures::HAS_OR)) ||
    // RDKit✔️✔️:         (features2 & static_cast<unsigned int>(QueryBoolFeatures::HAS_LOWAND) &&
    // RDKit✔️✔️:          features2 & static_cast<unsigned int>(QueryBoolFeatures::HAS_OR))) {
    // RDKit✔️✔️:       throw ValueErrorException(
    // RDKit✔️✔️:           "This is a non-smartable query - OR above and below AND in the "
    // RDKit✔️✔️:           "binary tree");
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     res += cs1;
    // RDKit✔️✔️:     if (!(cs1.empty() || cs2.empty())) {
    // RDKit✔️✔️:       res += ",";
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     res += cs2;
    // RDKit✔️✔️:     features |= static_cast<unsigned int>(QueryBoolFeatures::HAS_OR);
    // RDKit✔️✔️:   } else if ((descrip.find("And") > 0) &&
    // RDKit✔️✔️:              (descrip.find("And") < descrip.length())) {
    // RDKit✔️✔️:     std::string symb;
    // RDKit✔️✔️:     if (features1 & static_cast<unsigned int>(QueryBoolFeatures::HAS_OR) ||
    // RDKit✔️✔️:         features2 & static_cast<unsigned int>(QueryBoolFeatures::HAS_OR)) {
    // RDKit✔️✔️:       symb = ";";
    // RDKit✔️✔️:       features |= static_cast<unsigned int>(QueryBoolFeatures::HAS_LOWAND);
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       symb = "&";
    // RDKit✔️✔️:       features |= static_cast<unsigned int>(QueryBoolFeatures::HAS_AND);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     res += cs1;
    // RDKit✔️✔️:     if (!(cs1.empty() || cs2.empty())) {
    // RDKit✔️✔️:       res += symb;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     res += cs2;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     std::stringstream err;
    // RDKit✔️✔️:     err << "Don't know how to combine using " << descrip;
    // RDKit✔️✔️:     throw ValueErrorException(err.str());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   features |= features1;
    // RDKit✔️✔️:   features |= features2;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // Complexity review: both perform bounded description searches and one
    // output concatenation proportional to the two child strings. Rust moves
    // both child allocations into one pre-sized buffer and does no extra scan.
    let has_or = description.find("Or").is_some_and(|position| position > 0);
    let has_and = description.find("And").is_some_and(|position| position > 0);
    let separator = if has_or {
        if (features1.contains(QueryBoolFeatures::HAS_LOW_AND)
            && features1.contains(QueryBoolFeatures::HAS_OR))
            || (features2.contains(QueryBoolFeatures::HAS_LOW_AND)
                && features2.contains(QueryBoolFeatures::HAS_OR))
        {
            return Err(SmartsWriteError::OrAboveAndBelowAnd);
        }
        features.insert(QueryBoolFeatures::HAS_OR);
        ","
    } else if has_and {
        if features1.contains(QueryBoolFeatures::HAS_OR)
            || features2.contains(QueryBoolFeatures::HAS_OR)
        {
            features.insert(QueryBoolFeatures::HAS_LOW_AND);
            ";"
        } else {
            features.insert(QueryBoolFeatures::HAS_AND);
            "&"
        }
    } else {
        return Err(SmartsWriteError::UnknownCombination {
            description: description.to_owned(),
        });
    };

    let mut result = String::with_capacity(child1.len() + child2.len() + 1);
    result.push_str(&child1);
    if !child1.is_empty() && !child2.is_empty() {
        result.push_str(separator);
    }
    result.push_str(&child2);
    *features |= features1;
    *features |= features2;
    Ok(result)
}

fn describe_query<T>(query: &QueryNode<T>, leader: String) {
    // RDKit✔️✔️: template <typename T>
    // RDKit✔️✔️: void describeQuery(const T *query, std::string leader = "\t") {
    // RDKit✔️✔️:   // BOOST_LOG(rdInfoLog) << leader << query->getDescription() << std::endl;
    // RDKit✔️✔️:   typename T::CHILD_VECT_CI iter;
    // RDKit✔️✔️:   for (iter = query->beginChildren(); iter != query->endChildren(); ++iter) {
    // RDKit✔️✔️:     describeQuery(iter->get(), leader + "\t");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // Complexity review: both visit every query node once and allocate one
    // progressively longer leader string per traversed edge. The leader is
    // intentionally retained even though source logging is commented out.
    match query {
        QueryNode::Predicate(_) => {}
        QueryNode::And(children) | QueryNode::Or(children) | QueryNode::Xor(children) => {
            for child in children {
                describe_query(child, format!("{leader}\t"));
            }
        }
        QueryNode::Not(child) => describe_query(child, format!("{leader}\t")),
    }
}

fn range_prefix(data_function: AtomRangeDataFunction) -> &'static str {
    match data_function {
        AtomRangeDataFunction::ExplicitDegree => "D",
        AtomRangeDataFunction::NonHydrogenDegree => "d",
        AtomRangeDataFunction::TotalDegree => "X",
        AtomRangeDataFunction::TotalValence => "v",
        AtomRangeDataFunction::NumAtomRings => "R",
        AtomRangeDataFunction::NumHeteroatomNeighbors => "z",
        AtomRangeDataFunction::NumAliphaticHeteroatomNeighbors => "Z",
        AtomRangeDataFunction::MinRingSize => "r",
        AtomRangeDataFunction::RingBondCount => "x",
        AtomRangeDataFunction::ImplicitHydrogenCount => "h",
        AtomRangeDataFunction::FormalCharge => "+",
        AtomRangeDataFunction::NegativeFormalCharge => "-",
        AtomRangeDataFunction::AtomRingSize { .. } => "k",
    }
}

fn get_atom_smarts_simple(
    atom: &Atom,
    query: &AtomQueryPredicate,
    need_paren: &mut bool,
    check_for_symbol: bool,
    do_isomeric_smarts: bool,
    stereo_written: &mut bool,
) -> String {
    // BEGIN RDKIT CPP FUNCTION getAtomSmartsSimple
    // RDKit✔️✔️: std::string getAtomSmartsSimple(const QueryAtom *qatom,
    // RDKit✔️✔️:                                 const Atom::QUERYATOM_QUERY *query,
    // RDKit✔️✔️:                                 bool &needParen, bool checkForSymbol,
    // RDKit✔️✔️:                                 const SmilesWriteParams &) {
    // RDKit✔️✔️:   PRECONDITION(query, "bad query");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto *equery = dynamic_cast<const ATOM_EQUALS_QUERY *>(query);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::string descrip = query->getDescription();
    // RDKit✔️✔️:   bool hasVal = false;
    // RDKit✔️✔️:   enum class Modifiers : std::uint8_t {
    // RDKit✔️✔️:     NONE,
    // RDKit✔️✔️:     RANGE,
    // RDKit✔️✔️:     LESS,
    // RDKit✔️✔️:     GREATER
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:   Modifiers mods = Modifiers::NONE;
    // RDKit✔️✔️:   if (boost::starts_with(descrip, "range_")) {
    // RDKit✔️✔️:     mods = Modifiers::RANGE;
    // RDKit✔️✔️:     descrip = descrip.substr(6);
    // RDKit✔️✔️:   } else if (boost::starts_with(descrip, "less_")) {
    // RDKit✔️✔️:     mods = Modifiers::LESS;
    // RDKit✔️✔️:     descrip = descrip.substr(5);
    // RDKit✔️✔️:   } else if (boost::starts_with(descrip, "greater_")) {
    // RDKit✔️✔️:     mods = Modifiers::GREATER;
    // RDKit✔️✔️:     descrip = descrip.substr(8);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::stringstream res;
    // RDKit✔️✔️:   if (descrip == "AtomImplicitHCount") {
    // RDKit✔️✔️:     res << "h";
    // RDKit✔️✔️:     hasVal = true;
    // RDKit✔️✔️:     needParen = true;
    // RDKit✔️✔️:   } else if (descrip == "AtomHasImplicitH") {
    // RDKit✔️✔️:     res << "h";
    // RDKit✔️✔️:     needParen = true;
    // RDKit✔️✔️:   } else if (descrip == "AtomTotalValence") {
    // RDKit✔️✔️:     res << "v";
    // RDKit✔️✔️:     hasVal = true;
    // RDKit✔️✔️:     needParen = true;
    // RDKit✔️✔️:   } else if (descrip == "AtomAtomicNum") {
    // RDKit✔️✔️:     if (!qatom->hasProp(common_properties::smilesSymbol)) {
    // RDKit✔️✔️:       res << "#";
    // RDKit✔️✔️:       hasVal = true;
    // RDKit✔️✔️:       needParen = true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (descrip == "AtomExplicitDegree") {
    // RDKit✔️✔️:     res << "D";
    // RDKit✔️✔️:     hasVal = true;
    // RDKit✔️✔️:     needParen = true;
    // RDKit✔️✔️:   } else if (descrip == "AtomNonHydrogenDegree") {
    // RDKit✔️✔️:     res << "d";
    // RDKit✔️✔️:     hasVal = true;
    // RDKit✔️✔️:     needParen = true;
    // RDKit✔️✔️:   } else if (descrip == "AtomTotalDegree") {
    // RDKit✔️✔️:     res << "X";
    // RDKit✔️✔️:     hasVal = true;
    // RDKit✔️✔️:     needParen = true;
    // RDKit✔️✔️:   } else if (descrip == "AtomHasRingBond") {
    // RDKit✔️✔️:     res << "x";
    // RDKit✔️✔️:     needParen = true;
    // RDKit✔️✔️:   } else if (descrip == "AtomHCount") {
    // RDKit✔️✔️:     res << "H";
    // RDKit✔️✔️:     hasVal = true;
    // RDKit✔️✔️:     needParen = true;
    // RDKit✔️✔️:   } else if (descrip == "AtomIsAliphatic") {
    // RDKit✔️✔️:     res << "A";
    // RDKit✔️✔️:     needParen = false;
    // RDKit✔️✔️:   } else if (descrip == "AtomIsAromatic") {
    // RDKit✔️✔️:     res << "a";
    // RDKit✔️✔️:     needParen = false;
    // RDKit✔️✔️:   } else if (descrip == "AtomNull") {
    // RDKit✔️✔️:     res << "*";
    // RDKit✔️✔️:     needParen = false;
    // RDKit✔️✔️:   } else if (descrip == "AtomInRing") {
    // RDKit✔️✔️:     res << "R";
    // RDKit✔️✔️:     needParen = true;
    // RDKit✔️✔️:   } else if (descrip == "AtomMinRingSize") {
    // RDKit✔️✔️:     res << "r";
    // RDKit✔️✔️:     hasVal = true;
    // RDKit✔️✔️:     needParen = true;
    // RDKit✔️✔️:   } else if (descrip == "AtomRingSize") {
    // RDKit✔️✔️:     res << "k";
    // RDKit✔️✔️:     hasVal = true;
    // RDKit✔️✔️:     needParen = true;
    // RDKit✔️✔️:   } else if (descrip == "AtomInNRings") {
    // RDKit✔️✔️:     res << "R";
    // RDKit✔️✔️:     if (mods == Modifiers::NONE && equery && equery->getVal() >= 0) {
    // RDKit✔️✔️:       hasVal = true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     needParen = true;
    // RDKit✔️✔️:   } else if (descrip == "AtomHasHeteroatomNeighbors") {
    // RDKit✔️✔️:     res << "z";
    // RDKit✔️✔️:     needParen = true;
    // RDKit✔️✔️:   } else if (descrip == "AtomNumHeteroatomNeighbors") {
    // RDKit✔️✔️:     res << "z";
    // RDKit✔️✔️:     hasVal = true;
    // RDKit✔️✔️:     needParen = true;
    // RDKit✔️✔️:   } else if (descrip == "AtomHasAliphaticHeteroatomNeighbors") {
    // RDKit✔️✔️:     res << "Z";
    // RDKit✔️✔️:     needParen = true;
    // RDKit✔️✔️:   } else if (descrip == "AtomNumAliphaticHeteroatomNeighbors") {
    // RDKit✔️✔️:     res << "Z";
    // RDKit✔️✔️:     hasVal = true;
    // RDKit✔️✔️:     needParen = true;
    // RDKit✔️✔️:   } else if (descrip == "AtomFormalCharge") {
    // RDKit✔️✔️:     int val = equery ? equery->getVal() : 0;
    // RDKit✔️✔️:     if (val < 0) { res << "-"; } else { res << "+"; }
    // RDKit✔️✔️:     if (mods == Modifiers::NONE && abs(val) != 1) { res << abs(val); }
    // RDKit✔️✔️:     needParen = true;
    // RDKit✔️✔️:   } else if (descrip == "AtomNegativeFormalCharge") {
    // RDKit✔️✔️:     int val = equery ? equery->getVal() : 0;
    // RDKit✔️✔️:     if (val < 0) { res << "+"; } else { res << "-"; }
    // RDKit✔️✔️:     if (mods == Modifiers::NONE && abs(val) != 1) { res << abs(val); }
    // RDKit✔️✔️:     needParen = true;
    // RDKit✔️✔️:   } else if (descrip == "AtomHybridization" && equery) {
    // RDKit✔️✔️:     res << "^";
    // RDKit✔️✔️:     switch (equery->getVal()) {
    // RDKit✔️✔️:       case Atom::S: res << "0"; break;
    // RDKit✔️✔️:       case Atom::SP: res << "1"; break;
    // RDKit✔️✔️:       case Atom::SP2: res << "2"; break;
    // RDKit✔️✔️:       case Atom::SP3: res << "3"; break;
    // RDKit✔️✔️:       case Atom::SP3D: res << "4"; break;
    // RDKit✔️✔️:       case Atom::SP3D2: res << "5"; break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     needParen = true;
    // RDKit✔️✔️:   } else if (descrip == "AtomMass" && equery) {
    // RDKit✔️✔️:     res << equery->getVal() / massIntegerConversionFactor << "*";
    // RDKit✔️✔️:     needParen = true;
    // RDKit✔️✔️:   } else if (descrip == "AtomIsotope" && equery) {
    // RDKit✔️✔️:     res << equery->getVal() << "*";
    // RDKit✔️✔️:     needParen = true;
    // RDKit✔️✔️:   } else if (descrip == "AtomRingBondCount") {
    // RDKit✔️✔️:     res << "x"; hasVal = true; needParen = true;
    // RDKit✔️✔️:   } else if (descrip == "AtomUnsaturated") {
    // RDKit✔️✔️:     res << "$(*=,:,#*)"; needParen = true;
    // RDKit✔️✔️:   } else if (descrip == "AtomType" && equery) {
    // RDKit✔️✔️:     int atNum; bool isAromatic;
    // RDKit✔️✔️:     parseAtomType(equery->getVal(), atNum, isAromatic);
    // RDKit✔️✔️:     if (!checkForSymbol || !qatom->hasProp(common_properties::smilesSymbol)) {
    // RDKit✔️✔️:       std::string symbol = PeriodicTable::getTable()->getElementSymbol(atNum);
    // RDKit✔️✔️:       if (isAromatic) { symbol[0] += ('a' - 'A'); }
    // RDKit✔️✔️:       res << symbol;
    // RDKit✔️✔️:       if (!SmilesWrite::inOrganicSubset(atNum)) { needParen = true; }
    // RDKit✔️✔️:     } else { if (isAromatic) { res << "a"; } else { res << "A"; } }
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     BOOST_LOG(rdWarningLog) << "Cannot write SMARTS for query type : " << descrip
    // RDKit✔️✔️:                            << ". Ignoring it." << std::endl;
    // RDKit✔️✔️:     res << "*";
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (mods != Modifiers::NONE) {
    // RDKit✔️✔️:     res << "{";
    // RDKit✔️✔️:     const ATOM_RANGE_QUERY *rquery = nullptr;
    // RDKit✔️✔️:     switch (mods) {
    // RDKit✔️✔️:       case Modifiers::LESS: res << equery->getVal() << "-"; break;
    // RDKit✔️✔️:       case Modifiers::RANGE:
    // RDKit✔️✔️:         rquery = dynamic_cast<const ATOM_RANGE_QUERY *>(query);
    // RDKit✔️✔️:         CHECK_INVARIANT(rquery, "query could not be converted to range query");
    // RDKit✔️✔️:         res << ((const ATOM_RANGE_QUERY *)query)->getLower() << "-"
    // RDKit✔️✔️:             << ((const ATOM_RANGE_QUERY *)query)->getUpper(); break;
    // RDKit✔️✔️:       case Modifiers::GREATER: res << "-" << equery->getVal(); break;
    // RDKit✔️✔️:       default: break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     res << "}";
    // RDKit✔️✔️:   } else if (hasVal) { res << equery->getVal(); }
    // RDKit✔️✔️:   // handle atomic stereochemistry
    // RDKit✔️✔️:   if (qatom->hasOwningMol() &&
    // RDKit✔️✔️:       qatom->getOwningMol().hasProp(common_properties::_doIsoSmiles)) {
    // RDKit✔️✔️:     if (qatom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
    // RDKit✔️✔️:         !qatom->hasProp(_qatomHasStereoSet) &&
    // RDKit✔️✔️:         !qatom->hasProp(common_properties::_brokenChirality)) {
    // RDKit✔️✔️:       qatom->setProp(_qatomHasStereoSet, 1);
    // RDKit✔️✔️:       switch (qatom->getChiralTag()) {
    // RDKit✔️✔️:         case Atom::CHI_TETRAHEDRAL_CW: res << "@@"; needParen = true; break;
    // RDKit✔️✔️:         case Atom::CHI_TETRAHEDRAL_CCW: res << "@"; needParen = true; break;
    // RDKit✔️✔️:         default: break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res.str();
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION getAtomSmartsSimple
    // Complexity review: both dispatch on one typed leaf in O(1), allocate one
    // output string proportional to its representation, and perform no graph
    // traversal. Rust avoids RTTI and description-string manipulation.
    *need_paren = true;
    let mut result = match query {
        AtomQueryPredicate::Any => {
            *need_paren = false;
            "*".to_owned()
        }
        AtomQueryPredicate::AtomicNumber(value) => {
            if atom.prop("smilesSymbol").is_some() {
                String::new()
            } else {
                format!("#{value}")
            }
        }
        AtomQueryPredicate::AtomType {
            atomic_number,
            aromatic,
        } => {
            if check_for_symbol && atom.prop("smilesSymbol").is_some() {
                if *aromatic { "a".into() } else { "A".into() }
            } else {
                let mut symbol = Element::from_atomic_number(*atomic_number)
                    .map_or("*", Element::symbol)
                    .to_owned();
                if *aromatic {
                    symbol.make_ascii_lowercase();
                }
                if crate::notation::smiles_write::in_organic_subset(*atomic_number).unwrap_or(false)
                {
                    *need_paren = false;
                }
                symbol
            }
        }
        AtomQueryPredicate::ImplicitHydrogenCount(value) => format!("h{value}"),
        AtomQueryPredicate::HasImplicitHydrogen => "h".into(),
        AtomQueryPredicate::TotalValence(value) => format!("v{value}"),
        AtomQueryPredicate::ExplicitDegree(value) => format!("D{value}"),
        AtomQueryPredicate::NonHydrogenDegree(value) => format!("d{value}"),
        AtomQueryPredicate::TotalDegree(value) => format!("X{value}"),
        AtomQueryPredicate::HasRingBond => "x".into(),
        AtomQueryPredicate::HydrogenCount(value) => format!("H{value}"),
        AtomQueryPredicate::IsAromatic(false) => {
            *need_paren = false;
            "A".into()
        }
        AtomQueryPredicate::IsAromatic(true) => {
            *need_paren = false;
            "a".into()
        }
        AtomQueryPredicate::InRing => "R".into(),
        AtomQueryPredicate::SmallestRingSize(value) => format!("r{value}"),
        AtomQueryPredicate::InRingOfSize(value) => format!("k{value}"),
        AtomQueryPredicate::NumAtomRings(value) if *value >= 0 => format!("R{value}"),
        AtomQueryPredicate::NumAtomRings(_) => "R".into(),
        AtomQueryPredicate::HasHeteroatomNeighbors => "z".into(),
        AtomQueryPredicate::NumHeteroatomNeighbors(value) => format!("z{value}"),
        AtomQueryPredicate::HasAliphaticHeteroatomNeighbors => "Z".into(),
        AtomQueryPredicate::NumAliphaticHeteroatomNeighbors(value) => format!("Z{value}"),
        AtomQueryPredicate::FormalCharge(value) => {
            let sign = if *value < 0 { '-' } else { '+' };
            if value.unsigned_abs() == 1 {
                sign.to_string()
            } else {
                format!("{sign}{}", value.unsigned_abs())
            }
        }
        AtomQueryPredicate::NegativeFormalCharge(value) => {
            let sign = if *value < 0 { '+' } else { '-' };
            if value.unsigned_abs() == 1 {
                sign.to_string()
            } else {
                format!("{sign}{}", value.unsigned_abs())
            }
        }
        AtomQueryPredicate::HybridizationMatch(value) => format!(
            "^{}",
            match value {
                Hybridization::S => "0",
                Hybridization::Sp => "1",
                Hybridization::Sp2 => "2",
                Hybridization::Sp3 => "3",
                Hybridization::Sp3d => "4",
                Hybridization::Sp3d2 => "5",
                Hybridization::Unspecified | Hybridization::Sp2d | Hybridization::Other => "",
            }
        ),
        AtomQueryPredicate::Mass(value) => format!("{value}*"),
        AtomQueryPredicate::Isotope(value) => format!("{value}*"),
        AtomQueryPredicate::RingBondCount(value) => format!("x{value}"),
        AtomQueryPredicate::IsUnsaturated => "$(*=,:,#*)".into(),
        AtomQueryPredicate::Range(range) => {
            let (bounds, data_function) = range.writer_parts();
            let bounds = match bounds {
                // RDKit writes its left-hand query comparison literally:
                // `value >= observed` is `{-value}` and `value <= observed`
                // is `{value-}`.
                AtomRangeBounds::LessEqual(value) => format!("-{value}"),
                AtomRangeBounds::GreaterEqual(value) => format!("{value}-"),
                AtomRangeBounds::Inclusive { lower, upper, .. } => format!("{lower}-{upper}"),
            };
            format!("{}{{{bounds}}}", range_prefix(data_function))
        }
        AtomQueryPredicate::ExplicitDegreeLessEqual(value) => format!("D{{{value}-}}"),
        AtomQueryPredicate::NonHydrogenDegreeLessEqual(value) => format!("d{{{value}-}}"),
        AtomQueryPredicate::NonHydrogenDegreeGreaterEqual(value) => format!("d{{-{value}}}"),
        AtomQueryPredicate::TotalDegreeLessEqual(value) => format!("X{{{value}-}}"),
        AtomQueryPredicate::TotalDegreeGreaterEqual(value) => format!("X{{-{value}}}"),
        AtomQueryPredicate::TotalValenceLessEqual(value) => format!("v{{{value}-}}"),
        AtomQueryPredicate::TotalValenceGreaterEqual(value) => format!("v{{-{value}}}"),
        AtomQueryPredicate::RingBondCountLessEqual(value) => format!("x{{{value}-}}"),
        AtomQueryPredicate::ImplicitHydrogenCountLessEqual(value) => format!("h{{{value}-}}"),
        AtomQueryPredicate::InRingOfSizeLessEqual(value) => format!("k{{{value}-}}"),
        AtomQueryPredicate::InRingOfSizeGreaterEqual(value) => format!("k{{-{value}}}"),
        AtomQueryPredicate::SmallestRingSizeLessEqual(value) => format!("r{{{value}-}}"),
        AtomQueryPredicate::SmallestRingSizeGreaterEqual(value) => format!("r{{-{value}}}"),
        AtomQueryPredicate::DegreeLessEqual(value) => format!("D{{{value}-}}"),
        AtomQueryPredicate::DegreeGreaterEqual(value) => format!("D{{-{value}}}"),
        AtomQueryPredicate::RecursiveSmarts(_)
        | AtomQueryPredicate::AtomicNumberIn(_)
        | AtomQueryPredicate::AtomicNumberNotIn(_)
        | AtomQueryPredicate::NumRadicalElectrons(_)
        | AtomQueryPredicate::HasChiralTag
        | AtomQueryPredicate::MissingChiralTag
        | AtomQueryPredicate::ImplicitValence(_)
        | AtomQueryPredicate::ExplicitValence(_)
        | AtomQueryPredicate::HeavyAtomDegree(_)
        | AtomQueryPredicate::IsBridgehead
        | AtomQueryPredicate::HasProperty(_)
        | AtomQueryPredicate::PropertyValue { .. }
        | AtomQueryPredicate::RGroupLabel(_)
        | AtomQueryPredicate::MolFileAlias(_)
        | AtomQueryPredicate::ChiralTagMatch(_)
        | AtomQueryPredicate::ChiralPermutationMatch(_)
        | AtomQueryPredicate::UnsupportedFeature(_) => "*".into(),
    };

    if do_isomeric_smarts && !*stereo_written && atom.prop("_brokenChirality").is_none() {
        match atom.chiral_tag() {
            ChiralTag::TetrahedralCw => {
                result.push_str("@@");
                *need_paren = true;
                *stereo_written = true;
            }
            ChiralTag::TetrahedralCcw => {
                result.push('@');
                *need_paren = true;
                *stereo_written = true;
            }
            _ => {}
        }
    }
    result
}

fn get_recursive_structure_query_smarts<F>(
    query: &RecursiveStructureQuery,
    negated: bool,
    params: &SmilesWriteParams,
    write_molecule: F,
) -> Result<String, SmartsWriteError>
where
    F: FnOnce(&crate::QueryGraph, &SmilesWriteParams) -> Result<String, SmartsWriteError>,
{
    // BEGIN RDKIT CPP FUNCTION getRecursiveStructureQuerySmarts
    // RDKit✔️✔️: std::string getRecursiveStructureQuerySmarts(
    // RDKit✔️✔️:     const QueryAtom::QUERYATOM_QUERY *query, const SmilesWriteParams &params) {
    // RDKit✔️✔️:   PRECONDITION(query, "bad query");
    // RDKit✔️✔️:   PRECONDITION(query->getDescription() == "RecursiveStructure", "bad query");
    // RDKit✔️✔️:   const auto *rquery = dynamic_cast<const RecursiveStructureQuery *>(query);
    // RDKit✔️✔️:   PRECONDITION(rquery, "could not convert query to RecursiveStructureQuery");
    // RDKit✔️✔️:   auto *qmol = const_cast<ROMol *>(rquery->getQueryMol());
    // RDKit✔️✔️:   std::string res = MolToSmarts(*qmol, params);
    // RDKit✔️✔️:   res = "$(" + res + ")";
    // RDKit✔️✔️:   if (rquery->getNegation()) {
    // RDKit✔️✔️:     res = "!" + res;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION getRecursiveStructureQuerySmarts
    // Complexity review: both invoke the canonical molecule writer once and
    // build one output string linear in the nested SMARTS length. The generic
    // callback is monomorphized, so it adds no allocation or dynamic dispatch.
    let query_molecule = query
        .query_mol()
        .ok_or(SmartsWriteError::MissingRecursiveQueryMolecule)?;
    let inner = write_molecule(query_molecule, params)?;
    let mut result = String::with_capacity(inner.len() + if negated { 4 } else { 3 });
    if negated {
        result.push('!');
    }
    result.push_str("$(");
    result.push_str(&inner);
    result.push(')');
    Ok(result)
}

fn get_basic_bond_repr(
    bond_order: BondOrder,
    direction: BondDirection,
    reverse_dative: bool,
    params: &SmilesWriteParams,
) -> String {
    // BEGIN RDKIT CPP FUNCTION getBasicBondRepr
    // RDKit✔️✔️: std::string getBasicBondRepr(Bond::BondType typ, Bond::BondDir dir,
    // RDKit✔️✔️:                              bool reverseDative,
    // RDKit✔️✔️:                              const SmilesWriteParams &params) {
    // RDKit✔️✔️:   std::string res;
    // RDKit✔️✔️:   switch (typ) {
    // RDKit✔️✔️:     case Bond::SINGLE:
    // RDKit✔️✔️:       res = "-";
    // RDKit✔️✔️:       if (params.doIsomericSmiles) {
    // RDKit✔️✔️:         if (dir == Bond::ENDDOWNRIGHT) {
    // RDKit✔️✔️:           res = "\\";
    // RDKit✔️✔️:         } else if (dir == Bond::ENDUPRIGHT) {
    // RDKit✔️✔️:           res = "/";
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Bond::DOUBLE:
    // RDKit✔️✔️:       res = "=";
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Bond::TRIPLE:
    // RDKit✔️✔️:       res = "#";
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Bond::QUADRUPLE:
    // RDKit✔️✔️:       res = "$";
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Bond::AROMATIC:
    // RDKit✔️✔️:       res = ":";
    // RDKit✔️✔️:       if (params.doIsomericSmiles) {
    // RDKit✔️✔️:         if (dir == Bond::ENDDOWNRIGHT) {
    // RDKit✔️✔️:           res = "\\";
    // RDKit✔️✔️:         } else if (dir == Bond::ENDUPRIGHT) {
    // RDKit✔️✔️:           res = "/";
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Bond::DATIVE:
    // RDKit✔️✔️:       if (params.includeDativeBonds) {
    // RDKit✔️✔️:         if (reverseDative) {
    // RDKit✔️✔️:           res = "<-";
    // RDKit✔️✔️:         } else {
    // RDKit✔️✔️:           res = "->";
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         res = "-";
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Bond::ZERO:
    // RDKit✔️✔️:       res = "~";  // Actually means "any", but we use ~ for unknown bond types
    // RDKit✔️✔️:                   // in SMILES,
    // RDKit✔️✔️:       break;      // and this will match a ZOB.
    // RDKit✔️✔️:     default:
    // RDKit✔️✔️:       res = "";
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }  // namespace
    // END RDKIT CPP FUNCTION getBasicBondRepr
    // Complexity review: both perform one constant-size switch and allocate a
    // result of at most two bytes. There is no traversal, lookup, or cloning.
    match bond_order {
        BondOrder::Single | BondOrder::Aromatic if params.do_isomeric_smiles => match direction {
            BondDirection::EndDownRight => "\\".to_owned(),
            BondDirection::EndUpRight => "/".to_owned(),
            _ if bond_order == BondOrder::Single => "-".to_owned(),
            _ => ":".to_owned(),
        },
        BondOrder::Single => "-".to_owned(),
        BondOrder::Double => "=".to_owned(),
        BondOrder::Triple => "#".to_owned(),
        BondOrder::Quadruple => "$".to_owned(),
        BondOrder::Aromatic => ":".to_owned(),
        BondOrder::Dative if params.include_dative_bonds && reverse_dative => "<-".to_owned(),
        BondOrder::Dative if params.include_dative_bonds => "->".to_owned(),
        BondOrder::Dative => "-".to_owned(),
        BondOrder::Zero => "~".to_owned(),
        _ => String::new(),
    }
}

fn get_bond_smarts_simple(
    bond: &Bond,
    query: &BondQueryPredicate,
    atom_to_left_idx: Option<usize>,
    params: &SmilesWriteParams,
) -> Result<String, SmartsWriteError> {
    // BEGIN RDKIT CPP FUNCTION getBondSmartsSimple
    // RDKit✔️✔️: std::string getBondSmartsSimple(const Bond *bond,
    // RDKit✔️✔️:                                 const QueryBond::QUERYBOND_QUERY *bquery,
    // RDKit✔️✔️:                                 int atomToLeftIdx,
    // RDKit✔️✔️:                                 const SmilesWriteParams &params) {
    // RDKit✔️✔️:   PRECONDITION(bond, "bad bond");
    // RDKit✔️✔️:   PRECONDITION(bquery, "bad query");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto *equery = dynamic_cast<const BOND_EQUALS_QUERY *>(bquery);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::string descrip = bquery->getDescription();
    // RDKit✔️✔️:   std::string res = "";
    // RDKit✔️✔️:   if (descrip == "BondNull") {
    // RDKit✔️✔️:     res += "~";
    // RDKit✔️✔️:   } else if (descrip == "BondInRing") {
    // RDKit✔️✔️:     res += "@";
    // RDKit✔️✔️:   } else if (descrip == "SingleOrAromaticBond") {
    // RDKit✔️✔️:     auto dir = bond->getBondDir();
    // RDKit✔️✔️:     switch (dir) {
    // RDKit✔️✔️:       case Bond::ENDDOWNRIGHT: {
    // RDKit✔️✔️:         res += "\\";
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       case Bond::ENDUPRIGHT: {
    // RDKit✔️✔️:         res += "/";
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       default:
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (descrip == "SingleOrDoubleBond") {
    // RDKit✔️✔️:     res += "-,=";
    // RDKit✔️✔️:   } else if (descrip == "DoubleOrAromaticBond") {
    // RDKit✔️✔️:     res += "=,:";
    // RDKit✔️✔️:   } else if (descrip == "SingleOrDoubleOrAromaticBond") {
    // RDKit✔️✔️:     res += "-,=,:";
    // RDKit✔️✔️:   } else if (descrip == "BondDir" && equery) {
    // RDKit✔️✔️:     int val = equery->getVal();
    // RDKit✔️✔️:     if (val == static_cast<int>(Bond::ENDDOWNRIGHT)) {
    // RDKit✔️✔️:       res += "\\";
    // RDKit✔️✔️:     } else if (val == static_cast<int>(Bond::ENDUPRIGHT)) {
    // RDKit✔️✔️:       res += "/";
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       throw "Can't write smarts for this bond dir type";
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (descrip == "BondOrder" && equery) {
    // RDKit✔️✔️:     bool reverseDative =
    // RDKit✔️✔️:         (atomToLeftIdx >= 0 &&
    // RDKit✔️✔️:          bond->getBeginAtomIdx() != static_cast<unsigned int>(atomToLeftIdx));
    // RDKit✔️✔️:     res += getBasicBondRepr(static_cast<Bond::BondType>(equery->getVal()),
    // RDKit✔️✔️:                             bond->getBondDir(), reverseDative, params);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     std::stringstream msg;
    // RDKit✔️✔️:     msg << "Can't write smarts for this query bond type: " << descrip;
    // RDKit✔️✔️:     throw msg.str().c_str();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION getBondSmartsSimple
    // Complexity review: both perform constant-size typed dispatch and create
    // one result of at most five bytes. Fixed `OrderIn` vectors contain at most
    // three elements, so their comparisons preserve the source's O(1) cost.
    match query {
        BondQueryPredicate::Any => Ok("~".to_owned()),
        BondQueryPredicate::IsInRing(_) => Ok("@".to_owned()),
        BondQueryPredicate::OrderIn(orders)
            if orders.as_slice() == [BondOrder::Single, BondOrder::Aromatic] =>
        {
            Ok(match bond.direction() {
                BondDirection::EndDownRight => "\\".to_owned(),
                BondDirection::EndUpRight => "/".to_owned(),
                _ => String::new(),
            })
        }
        BondQueryPredicate::OrderIn(orders)
            if orders.as_slice() == [BondOrder::Single, BondOrder::Double] =>
        {
            Ok("-,=".to_owned())
        }
        BondQueryPredicate::OrderIn(orders)
            if orders.as_slice() == [BondOrder::Double, BondOrder::Aromatic] =>
        {
            Ok("=,:".to_owned())
        }
        BondQueryPredicate::OrderIn(orders)
            if orders.as_slice() == [BondOrder::Single, BondOrder::Double, BondOrder::Aromatic] =>
        {
            Ok("-,=,:".to_owned())
        }
        BondQueryPredicate::Direction(BondDirection::EndDownRight) => Ok("\\".to_owned()),
        BondQueryPredicate::Direction(BondDirection::EndUpRight) => Ok("/".to_owned()),
        BondQueryPredicate::Direction(direction) => {
            Err(SmartsWriteError::UnsupportedBondDirection {
                direction: *direction,
            })
        }
        BondQueryPredicate::Order(order) => {
            let reverse_dative =
                atom_to_left_idx.is_some_and(|atom_idx| bond.begin().index() != atom_idx);
            Ok(get_basic_bond_repr(
                *order,
                bond.direction(),
                reverse_dative,
                params,
            ))
        }
        _ => Err(SmartsWriteError::UnsupportedBondQuery {
            predicate: query.clone(),
        }),
    }
}

fn atom_query_without_not<'a>(
    mut node: &'a QueryNode<AtomQueryPredicate>,
    mut negate: bool,
) -> (&'a QueryNode<AtomQueryPredicate>, bool) {
    while let QueryNode::Not(child) = node {
        negate = !negate;
        node = child;
    }
    (node, negate)
}

fn recurse_get_smarts<F>(
    atom: &Atom,
    node: &QueryNode<AtomQueryPredicate>,
    negate: bool,
    features: &mut QueryBoolFeatures,
    params: &SmilesWriteParams,
    stereo_written: &mut bool,
    write_molecule: &mut F,
) -> Result<String, SmartsWriteError>
where
    F: FnMut(&crate::QueryGraph, &SmilesWriteParams) -> Result<String, SmartsWriteError>,
{
    // BEGIN RDKIT CPP FUNCTION _recurseGetSmarts
    // RDKit✔️✔️: std::string _recurseGetSmarts(const QueryAtom *qatom,
    // RDKit✔️✔️:                               const QueryAtom::QUERYATOM_QUERY *node,
    // RDKit✔️✔️:                               bool negate, unsigned int &features,
    // RDKit✔️✔️:                               const SmilesWriteParams &params) {
    // RDKit✔️✔️:   PRECONDITION(node, "bad node");
    // RDKit✔️✔️:   // the algorithm goes something like this
    // RDKit✔️✔️:   // - recursively get the smarts for the child queries
    // RDKit✔️✔️:   // - combine the child smarts using the following rules:
    // RDKit✔️✔️:   //      - if we are currently at an OR query, combine the subqueries with a
    // RDKit✔️✔️:   //      ",",
    // RDKit✔️✔️:   //        but only if neither of child smarts do not contain "," and ";"
    // RDKit✔️✔️:   //        This situation leads to a no smartable situation and throw an
    // RDKit✔️✔️:   //        error
    // RDKit✔️✔️:   //      - if we are currently at an and query, combine the child smarts with
    // RDKit✔️✔️:   //      "&"
    // RDKit✔️✔️:   //        if neither of the child smarts contain "," - otherwise combine
    // RDKit✔️✔️:   //        them
    // RDKit✔️✔️:   //        the child smarts with a ";"
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   // There is an additional complication with composite nodes that carry a
    // RDKit✔️✔️:   // negation - in this
    // RDKit✔️✔️:   // case we will propagate the negation to the child nodes using the
    // RDKit✔️✔️:   // following rules
    // RDKit✔️✔️:   //   NOT (a AND b) = ( NOT (a)) AND ( NOT (b))
    // RDKit✔️✔️:   //   NOT (a OR b) = ( NOT (a)) OR ( NOT (b))
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto descrip = node->getDescription();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   unsigned int child1Features = 0;
    // RDKit✔️✔️:   unsigned int child2Features = 0;
    // RDKit✔️✔️:   auto chi = node->beginChildren();
    // RDKit✔️✔️:   auto child1 = chi->get();
    // RDKit✔️✔️:   auto dsc1 = child1->getDescription();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   ++chi;
    // RDKit✔️✔️:   CHECK_INVARIANT(chi != node->endChildren(),
    // RDKit✔️✔️:                   "Not enough children on the query");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   bool needParen;
    // RDKit✔️✔️:   std::string csmarts1;
    // RDKit✔️✔️:   // deal with the first child
    // RDKit✔️✔️:   if (dsc1 == "RecursiveStructure") {
    // RDKit✔️✔️:     csmarts1 = getRecursiveStructureQuerySmarts(child1, params);
    // RDKit✔️✔️:     features |= static_cast<unsigned int>(QueryBoolFeatures::HAS_RECURSION);
    // RDKit✔️✔️:   } else if ((dsc1 != "AtomOr") && (dsc1 != "AtomAnd")) {
    // RDKit✔️✔️:     // child 1 is a simple node, but we only check for the smilesSymbol
    // RDKit✔️✔️:     //  if descrip=="AtomAnd"
    // RDKit✔️✔️:     csmarts1 = getAtomSmartsSimple(qatom, child1, needParen,
    // RDKit✔️✔️:                                    descrip == "AtomAnd", params);
    // RDKit✔️✔️:     bool nneg = (negate) ^ (child1->getNegation());
    // RDKit✔️✔️:     if (nneg) {
    // RDKit✔️✔️:       csmarts1 = "!" + csmarts1;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     // child 1 is composite node - recurse
    // RDKit✔️✔️:     bool nneg = (negate) ^ (child1->getNegation());
    // RDKit✔️✔️:     csmarts1 = _recurseGetSmarts(qatom, child1, nneg, child1Features, params);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // ok if we have a negation and we have an OR , we have to change to
    // RDKit✔️✔️:   // an AND since we propagated the negation
    // RDKit✔️✔️:   // i.e NOT (A OR B) = (NOT (A)) AND (NOT(B))
    // RDKit✔️✔️:   if (negate) {
    // RDKit✔️✔️:     if (descrip == "AtomOr") {
    // RDKit✔️✔️:       descrip = "AtomAnd";
    // RDKit✔️✔️:     } else if (descrip == "AtomAnd") {
    // RDKit✔️✔️:       descrip = "AtomOr";
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   auto res = csmarts1;
    // RDKit✔️✔️:   while (chi != node->endChildren()) {
    // RDKit✔️✔️:     auto child2 = chi->get();
    // RDKit✔️✔️:     ++chi;
    // RDKit✔️✔️:
    // RDKit✔️✔️:     auto dsc2 = child2->getDescription();
    // RDKit✔️✔️:     std::string csmarts2;
    // RDKit✔️✔️:
    // RDKit✔️✔️:     // deal with the next child
    // RDKit✔️✔️:     if (dsc2 == "RecursiveStructure") {
    // RDKit✔️✔️:       csmarts2 = getRecursiveStructureQuerySmarts(child2, params);
    // RDKit✔️✔️:       features |= static_cast<unsigned int>(QueryBoolFeatures::HAS_RECURSION);
    // RDKit✔️✔️:     } else if ((dsc2 != "AtomOr") && (dsc2 != "AtomAnd")) {
    // RDKit✔️✔️:       // child 2 is a simple node
    // RDKit✔️✔️:       csmarts2 = getAtomSmartsSimple(qatom, child2, needParen, false, params);
    // RDKit✔️✔️:       bool nneg = (negate) ^ (child2->getNegation());
    // RDKit✔️✔️:       if (nneg) {
    // RDKit✔️✔️:         csmarts2 = "!" + csmarts2;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       bool nneg = (negate) ^ (child2->getNegation());
    // RDKit✔️✔️:       csmarts2 = _recurseGetSmarts(qatom, child2, nneg, child2Features, params);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     res = _combineChildSmarts(res, child1Features, csmarts2, child2Features,
    // RDKit✔️✔️:                               descrip, features);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION _recurseGetSmarts
    // Complexity review: each query node is visited once and each rendered
    // child string participates in the same linear concatenation sequence as
    // RDKit. Typed enum dispatch removes description-string copies; recursion
    // remains O(n) time and O(h) stack with output-proportional allocation.
    let (node, negate) = atom_query_without_not(node, negate);
    let (children, mut description) = match node {
        QueryNode::And(children) => (children, "AtomAnd"),
        QueryNode::Or(children) => (children, "AtomOr"),
        QueryNode::Xor(_) => return Err(SmartsWriteError::XorComposite),
        QueryNode::Predicate(_) | QueryNode::Not(_) => {
            return Err(SmartsWriteError::CompositeChildCount { kind: "atom" });
        }
    };
    if children.len() < 2 {
        return Err(SmartsWriteError::CompositeChildCount { kind: "atom" });
    }
    if negate {
        description = if description == "AtomOr" {
            "AtomAnd"
        } else {
            "AtomOr"
        };
    }

    let render_child = |child: &QueryNode<AtomQueryPredicate>,
                        check_for_symbol: bool,
                        child_features: &mut QueryBoolFeatures,
                        stereo_written: &mut bool,
                        write_molecule: &mut F|
     -> Result<String, SmartsWriteError> {
        let (child, child_negate) = atom_query_without_not(child, negate);
        match child {
            QueryNode::Predicate(AtomQueryPredicate::RecursiveSmarts(recursive)) => {
                child_features.insert(QueryBoolFeatures::HAS_RECURSION);
                get_recursive_structure_query_smarts(
                    recursive,
                    child_negate,
                    params,
                    write_molecule,
                )
            }
            QueryNode::Predicate(predicate) => {
                let mut need_paren = false;
                let mut result = get_atom_smarts_simple(
                    atom,
                    predicate,
                    &mut need_paren,
                    check_for_symbol,
                    params.do_isomeric_smiles,
                    stereo_written,
                );
                if child_negate {
                    result.insert(0, '!');
                }
                Ok(result)
            }
            QueryNode::And(_) | QueryNode::Or(_) | QueryNode::Xor(_) | QueryNode::Not(_) => {
                recurse_get_smarts(
                    atom,
                    child,
                    child_negate,
                    child_features,
                    params,
                    stereo_written,
                    write_molecule,
                )
            }
        }
    };

    let mut first_features = QueryBoolFeatures::default();
    let mut result = render_child(
        &children[0],
        matches!(node, QueryNode::And(_)),
        &mut first_features,
        stereo_written,
        write_molecule,
    )?;
    *features |= first_features;
    for child in &children[1..] {
        let mut child_features = QueryBoolFeatures::default();
        let child_smarts = render_child(
            child,
            false,
            &mut child_features,
            stereo_written,
            write_molecule,
        )?;
        result = combine_child_smarts(
            result,
            first_features,
            child_smarts,
            child_features,
            description,
            features,
        )?;
    }
    Ok(result)
}

fn bond_query_without_not<'a>(
    mut node: &'a QueryNode<BondQueryPredicate>,
    mut negate: bool,
) -> (&'a QueryNode<BondQueryPredicate>, bool) {
    while let QueryNode::Not(child) = node {
        negate = !negate;
        node = child;
    }
    (node, negate)
}

fn recurse_bond_smarts(
    bond: &Bond,
    node: &QueryNode<BondQueryPredicate>,
    negate: bool,
    atom_to_left_idx: Option<usize>,
    features: &mut QueryBoolFeatures,
    params: &SmilesWriteParams,
) -> Result<String, SmartsWriteError> {
    // BEGIN RDKIT CPP FUNCTION _recurseBondSmarts
    // RDKit✔️✔️: std::string _recurseBondSmarts(const Bond *bond,
    // RDKit✔️✔️:                                const QueryBond::QUERYBOND_QUERY *node,
    // RDKit✔️✔️:                                bool negate, int atomToLeftIdx,
    // RDKit✔️✔️:                                unsigned int &features,
    // RDKit✔️✔️:                                const SmilesWriteParams &params) {
    // RDKit✔️✔️:   // the algorithm goes something like this
    // RDKit✔️✔️:   // - recursively get the smarts for the child query bonds
    // RDKit✔️✔️:   // - combine the child smarts using the following rules:
    // RDKit✔️✔️:   //      - if we are currently at an OR query, combine the subqueries with a
    // RDKit✔️✔️:   //      ",",
    // RDKit✔️✔️:   //        but only if neither of child smarts do not contain "," and ";"
    // RDKit✔️✔️:   //        This situation leads to a no smartable situation and throw an
    // RDKit✔️✔️:   //        error
    // RDKit✔️✔️:   //      - if we are currently at an and query, combine the child smarts with
    // RDKit✔️✔️:   //      "&"
    // RDKit✔️✔️:   //        if neither of the child smarts contain "," - otherwise combine
    // RDKit✔️✔️:   //        them
    // RDKit✔️✔️:   //        the child smarts with a ";"
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   // There is an additional complication with composite nodes that carry a
    // RDKit✔️✔️:   // negation - in this
    // RDKit✔️✔️:   // case we will propagate the negation to the child nodes using the
    // RDKit✔️✔️:   // following rules
    // RDKit✔️✔️:   //   NOT (a AND b) = ( NOT (a)) AND ( NOT (b))
    // RDKit✔️✔️:   //   NOT (a OR b) = ( NOT (a)) OR ( NOT (b))
    // RDKit✔️✔️:   PRECONDITION(bond, "bad bond");
    // RDKit✔️✔️:   PRECONDITION(node, "bad node");
    // RDKit✔️✔️:   std::string descrip = node->getDescription();
    // RDKit✔️✔️:   std::string res = "";
    // RDKit✔️✔️:
    // RDKit✔️✔️:   const QueryBond::QUERYBOND_QUERY *child1;
    // RDKit✔️✔️:   const QueryBond::QUERYBOND_QUERY *child2;
    // RDKit✔️✔️:   unsigned int child1Features = 0;
    // RDKit✔️✔️:   unsigned int child2Features = 0;
    // RDKit✔️✔️:   QueryBond::QUERYBOND_QUERY::CHILD_VECT_CI chi;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   chi = node->beginChildren();
    // RDKit✔️✔️:   child1 = chi->get();
    // RDKit✔️✔️:   chi++;
    // RDKit✔️✔️:   child2 = chi->get();
    // RDKit✔️✔️:   chi++;
    // RDKit✔️✔️:   // OK we should be at the end of vector by now - since we can have only two
    // RDKit✔️✔️:   // children,
    // RDKit✔️✔️:   // well - at least in this case
    // RDKit✔️✔️:   CHECK_INVARIANT(chi == node->endChildren(), "Too many children on the query");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::string dsc1, dsc2;
    // RDKit✔️✔️:   dsc1 = child1->getDescription();
    // RDKit✔️✔️:   dsc2 = child2->getDescription();
    // RDKit✔️✔️:   std::string csmarts1, csmarts2;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if ((dsc1 != "BondOr") && (dsc1 != "BondAnd")) {
    // RDKit✔️✔️:     // child1 is  simple node get the smarts directly
    // RDKit✔️✔️:     const auto *tchild = static_cast<const BOND_EQUALS_QUERY *>(child1);
    // RDKit✔️✔️:     csmarts1 = getBondSmartsSimple(bond, tchild, atomToLeftIdx, params);
    // RDKit✔️✔️:     bool nneg = (negate) ^ (tchild->getNegation());
    // RDKit✔️✔️:     if (nneg) {
    // RDKit✔️✔️:       csmarts1 = "!" + csmarts1;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     // child1 is a composite node recurse further
    // RDKit✔️✔️:     bool nneg = (negate) ^ (child1->getNegation());
    // RDKit✔️✔️:     csmarts1 = _recurseBondSmarts(bond, child1, nneg, atomToLeftIdx,
    // RDKit✔️✔️:                                   child1Features, params);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // now deal with the second child
    // RDKit✔️✔️:   if ((dsc2 != "BondOr") && (dsc2 != "BondAnd")) {
    // RDKit✔️✔️:     // child 2 is a simple node
    // RDKit✔️✔️:     const auto *tchild = static_cast<const BOND_EQUALS_QUERY *>(child2);
    // RDKit✔️✔️:     csmarts2 = getBondSmartsSimple(bond, tchild, atomToLeftIdx, params);
    // RDKit✔️✔️:     bool nneg = (negate) ^ (tchild->getNegation());
    // RDKit✔️✔️:     if (nneg) {
    // RDKit✔️✔️:       csmarts2 = "!" + csmarts2;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     // child two is a composite node - recurse
    // RDKit✔️✔️:     bool nneg = (negate) ^ (child2->getNegation());
    // RDKit✔️✔️:     csmarts1 = _recurseBondSmarts(bond, child2, nneg, atomToLeftIdx,
    // RDKit✔️✔️:                                   child2Features, params);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // ok if we have a negation and we have to change the underlying logic,
    // RDKit✔️✔️:   // since we propagated the negation i.e NOT (A OR B) = (NOT (A)) AND
    // RDKit✔️✔️:   // (NOT(B))
    // RDKit✔️✔️:   if (negate) {
    // RDKit✔️✔️:     if (descrip == "BondOr") {
    // RDKit✔️✔️:       descrip = "BondAnd";
    // RDKit✔️✔️:     } else if (descrip == "BondAnd") {
    // RDKit✔️✔️:       descrip = "BondOr";
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   res += _combineChildSmarts(csmarts1, child1Features, csmarts2, child2Features,
    // RDKit✔️✔️:                              descrip, features);
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION _recurseBondSmarts
    // Complexity review: the source requires exactly two children, and Rust
    // preserves that constant-size traversal. Each nested node is visited once;
    // recursion is O(n) time and O(h) stack with output-proportional strings.
    let (node, negate) = bond_query_without_not(node, negate);
    let (children, mut description) = match node {
        QueryNode::And(children) => (children, "BondAnd"),
        QueryNode::Or(children) => (children, "BondOr"),
        QueryNode::Xor(_) => return Err(SmartsWriteError::XorComposite),
        QueryNode::Predicate(_) | QueryNode::Not(_) => {
            return Err(SmartsWriteError::CompositeChildCount { kind: "bond" });
        }
    };
    if children.len() != 2 {
        return Err(SmartsWriteError::CompositeChildCount { kind: "bond" });
    }
    if negate {
        description = if description == "BondOr" {
            "BondAnd"
        } else {
            "BondOr"
        };
    }

    let render_child = |child: &QueryNode<BondQueryPredicate>,
                        child_features: &mut QueryBoolFeatures|
     -> Result<String, SmartsWriteError> {
        let (child, child_negate) = bond_query_without_not(child, negate);
        match child {
            QueryNode::Predicate(predicate) => {
                let mut result = get_bond_smarts_simple(bond, predicate, atom_to_left_idx, params)?;
                if child_negate {
                    result.insert(0, '!');
                }
                Ok(result)
            }
            QueryNode::And(_) | QueryNode::Or(_) | QueryNode::Xor(_) | QueryNode::Not(_) => {
                recurse_bond_smarts(
                    bond,
                    child,
                    child_negate,
                    atom_to_left_idx,
                    child_features,
                    params,
                )
            }
        }
    };

    let mut child1_features = QueryBoolFeatures::default();
    let mut child2_features = QueryBoolFeatures::default();
    let mut child1_smarts = render_child(&children[0], &mut child1_features)?;
    let child2_base = bond_query_without_not(&children[1], false).0;
    let child2_smarts = if matches!(child2_base, QueryNode::And(_) | QueryNode::Or(_)) {
        // Preserve the pinned source assignment to csmarts1 in this branch.
        child1_smarts = render_child(&children[1], &mut child2_features)?;
        String::new()
    } else {
        render_child(&children[1], &mut child2_features)?
    };
    combine_child_smarts(
        child1_smarts,
        child1_features,
        child2_smarts,
        child2_features,
        description,
        features,
    )
}

fn fragment_smarts_construct<FA, FB>(
    molecule: &mut Molecule,
    atom_idx: AtomId,
    ranks: &[usize],
    params: &SmilesWriteParams,
    atom_ordering: &mut Vec<AtomId>,
    bond_ordering: &mut Vec<BondId>,
    atoms_in_play: &[AtomId],
    bonds_in_play: &[BondId],
    mut write_atom: FA,
    mut write_bond: FB,
) -> Result<String, SmartsWriteError>
where
    FA: FnMut(&Atom, &SmilesWriteParams) -> Result<String, SmartsWriteError>,
    FB: FnMut(&Bond, &SmilesWriteParams, AtomId) -> Result<String, SmartsWriteError>,
{
    // BEGIN RDKIT CPP FUNCTION FragmentSmartsConstruct
    // RDKit✔️✔️: std::string FragmentSmartsConstruct(
    // RDKit✔️✔️:     ROMol &mol, unsigned int atomIdx, std::vector<Canon::AtomColors> &colors,
    // RDKit✔️✔️:     UINT_VECT &ranks, const SmilesWriteParams &params,
    // RDKit✔️✔️:     std::vector<unsigned int> &atomOrdering,
    // RDKit✔️✔️:     std::vector<unsigned int> &bondOrdering,
    // RDKit✔️✔️:     const boost::dynamic_bitset<> &atomsInPlay,
    // RDKit✔️✔️:     const boost::dynamic_bitset<> *bondsInPlay) {
    // RDKit✔️✔️:   // this is dirty trick get around the fact that canonicalizeFragment
    // RDKit✔️✔️:   // thinks we already called findSSSR - to do some atom ranking
    // RDKit✔️✔️:   // but for smarts we are going to ignore that part. We will artificially
    // RDKit✔️✔️:   // set the "SSSR" property to an empty property
    // RDKit✔️✔️:
    // RDKit✔️✔️:   mol.getRingInfo()->reset();
    // RDKit✔️✔️:   mol.getRingInfo()->initialize(FIND_RING_TYPE_SYMM_SSSR);
    // RDKit✔️✔️:   for (auto &atom : mol.atoms()) {
    // RDKit✔️✔️:     atom->updatePropertyCache(false);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // For Smarts, we avoid reordering of chiral atoms in canonicalizeFragment.
    // RDKit✔️✔️:   bool doRandom = false;
    // RDKit✔️✔️:   bool doChiralInversions = true;
    // RDKit✔️✔️:   Canon::MolStack molStack;
    // RDKit✔️✔️:   molStack.reserve(mol.getNumAtoms() + mol.getNumBonds());
    // RDKit✔️✔️:   Canon::canonicalizeFragment(
    // RDKit✔️✔️:       mol, atomIdx, colors, ranks, molStack, &atomsInPlay, bondsInPlay, nullptr,
    // RDKit✔️✔️:       params.doIsomericSmiles, doRandom, doChiralInversions);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // now clear the "SSSR" property
    // RDKit✔️✔️:   mol.getRingInfo()->reset();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::stringstream res;
    // RDKit✔️✔️:   for (auto &msCI : molStack) {
    // RDKit✔️✔️:     switch (msCI.type) {
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
    // RDKit✔️✔️:       case Canon::MOL_STACK_RING: {
    // RDKit✔️✔️:         if (msCI.number < 10) {
    // RDKit✔️✔️:           res << msCI.number;
    // RDKit✔️✔️:         } else {
    // RDKit✔️✔️:           res << "%" << msCI.number;
    // RDKit✔️✔️:         }
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
    // RDKit✔️✔️:       default:
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res.str();
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION FragmentSmartsConstruct
    // Complexity review: both update each atom cache once, perform the shared
    // O(V + E) canonical traversal once, and emit each stack element once.
    // The adapter allocates only the same output stack and output string; it
    // neither clones the molecule nor repeats DFS, ring, or branch discovery.
    molecule.derived_cache_mut().rings = Some(RingInfo::new(
        RingFindType::SymmSssr,
        molecule.num_atoms(),
        molecule.num_bonds(),
    ));
    if let Err(source) = crate::notation::smiles_write::update_property_cache_for_smiles(molecule) {
        molecule.derived_cache_mut().rings = None;
        return Err(SmartsWriteError::Traversal { source });
    }
    let stack = crate::notation::smiles_write::canonicalize_fragment_stack_for_smarts(
        molecule,
        atoms_in_play,
        bonds_in_play,
        atom_idx,
        ranks,
        params,
    )
    .map_err(|source| SmartsWriteError::Traversal { source });
    molecule.derived_cache_mut().rings = None;
    let stack = stack?;

    let mut result = String::new();
    for item in stack {
        match item {
            crate::notation::smiles_write::MolStackElem::Atom(atom_id) => {
                result.push_str(&write_atom(&molecule.atoms()[atom_id.index()], params)?);
                atom_ordering.push(atom_id);
            }
            crate::notation::smiles_write::MolStackElem::Bond(bond_id, atom_to_left) => {
                result.push_str(&write_bond(
                    &molecule.bonds()[bond_id.index()],
                    params,
                    atom_to_left,
                )?);
                bond_ordering.push(bond_id);
            }
            crate::notation::smiles_write::MolStackElem::Ring { ring_idx, .. } => {
                if ring_idx < 10 {
                    result.push_str(&ring_idx.to_string());
                } else {
                    result.push('%');
                    result.push_str(&ring_idx.to_string());
                }
            }
            crate::notation::smiles_write::MolStackElem::BranchOpen => result.push('('),
            crate::notation::smiles_write::MolStackElem::BranchClose => result.push(')'),
        }
    }
    Ok(result)
}

fn get_non_query_atom_smarts(
    atom: &Atom,
    do_isomeric_smarts: bool,
    stereo_written: &mut bool,
) -> String {
    // BEGIN RDKIT CPP FUNCTION getNonQueryAtomSmarts
    // RDKit✔️✔️: std::string getNonQueryAtomSmarts(const Atom *atom) {
    // RDKit✔️✔️:   PRECONDITION(atom, "bad atom");
    // RDKit✔️✔️:   PRECONDITION(!atom->hasQuery(), "atom should not have query");
    // RDKit✔️✔️:   std::stringstream res;
    // RDKit✔️✔️:   res << "[";
    // RDKit✔️✔️:
    // RDKit✔️✔️:   int isotope = atom->getIsotope();
    // RDKit✔️✔️:   if (isotope) {
    // RDKit✔️✔️:     res << isotope;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::string symbol;
    // RDKit✔️✔️:   if (atom->getPropIfPresent(common_properties::smilesSymbol, symbol)) {
    // RDKit✔️✔️:     res << symbol;
    // RDKit✔️✔️:   } else if (SmilesWrite::inOrganicSubset(atom->getAtomicNum())) {
    // RDKit✔️✔️:     res << "#" << atom->getAtomicNum();
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     res << atom->getSymbol();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   bool addedChirality = false;
    // RDKit✔️✔️:   if (atom->hasOwningMol() &&
    // RDKit✔️✔️:       atom->getOwningMol().hasProp(common_properties::_doIsoSmiles)) {
    // RDKit✔️✔️:     if (atom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
    // RDKit✔️✔️:         !atom->hasProp(_qatomHasStereoSet) &&
    // RDKit✔️✔️:         !atom->hasProp(common_properties::_brokenChirality)) {
    // RDKit✔️✔️:       atom->setProp(_qatomHasStereoSet, 1);
    // RDKit✔️✔️:       switch (atom->getChiralTag()) {
    // RDKit✔️✔️:         case Atom::CHI_TETRAHEDRAL_CW:
    // RDKit✔️✔️:           res << "@@";
    // RDKit✔️✔️:           addedChirality = true;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         case Atom::CHI_TETRAHEDRAL_CCW:
    // RDKit✔️✔️:           res << "@";
    // RDKit✔️✔️:           addedChirality = true;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         default:
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (addedChirality && atom->getNumExplicitHs() == 1) {
    // RDKit✔️✔️:     // FIX: this isn't really correct in many cases, but
    // RDKit✔️✔️:     //   fixing it requires opening a fairly large construction site on the
    // RDKit✔️✔️:     //   SMARTS handling side. We'll do this later.
    // RDKit✔️✔️:     res << "H";
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   auto chg = atom->getFormalCharge();
    // RDKit✔️✔️:   if (chg) {
    // RDKit✔️✔️:     if (chg == -1) {
    // RDKit✔️✔️:       res << "-";
    // RDKit✔️✔️:     } else if (chg == 1) {
    // RDKit✔️✔️:       res << "+";
    // RDKit✔️✔️:     } else if (chg < 0) {
    // RDKit✔️✔️:       res << atom->getFormalCharge();
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       res << "+" << atom->getFormalCharge();
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   int mapNum;
    // RDKit✔️✔️:   if (atom->getPropIfPresent(common_properties::molAtomMapNumber, mapNum)) {
    // RDKit✔️✔️:     res << ":";
    // RDKit✔️✔️:     res << mapNum;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   res << "]";
    // RDKit✔️✔️:   return res.str();
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION getNonQueryAtomSmarts
    // Complexity review: both perform constant-time typed field/property
    // reads and append to one bounded atom string. The caller-owned flag
    // replaces RDKit's temporary atom property without extra allocation.
    let mut result = String::from("[");
    if let Some(isotope) = atom.isotope() {
        result.push_str(&isotope.to_string());
    }
    if let Some(symbol) = atom.prop("smilesSymbol") {
        result.push_str(symbol);
    } else if matches!(
        atom.atomic_number(),
        0 | 5 | 6 | 7 | 8 | 9 | 15 | 16 | 17 | 35 | 53
    ) {
        result.push('#');
        result.push_str(&atom.atomic_number().to_string());
    } else {
        result.push_str(atom.element().symbol());
    }

    let mut added_chirality = false;
    if do_isomeric_smarts && !*stereo_written && atom.prop("_brokenChirality").is_none() {
        match atom.chiral_tag() {
            ChiralTag::TetrahedralCw => {
                result.push_str("@@");
                added_chirality = true;
                *stereo_written = true;
            }
            ChiralTag::TetrahedralCcw => {
                result.push('@');
                added_chirality = true;
                *stereo_written = true;
            }
            _ => {}
        }
    }
    if added_chirality && atom.explicit_hydrogens() == 1 {
        result.push('H');
    }
    match atom.formal_charge() {
        -1 => result.push('-'),
        1 => result.push('+'),
        charge if charge < 0 => result.push_str(&charge.to_string()),
        charge if charge > 0 => {
            result.push('+');
            result.push_str(&charge.to_string());
        }
        _ => {}
    }
    if let Some(atom_map) = atom.atom_map() {
        result.push(':');
        result.push_str(&atom_map.to_string());
    }
    result.push(']');
    result
}

fn get_non_query_bond_smarts(
    bond: &Bond,
    atom_to_left_idx: Option<usize>,
    params: &SmilesWriteParams,
) -> String {
    // BEGIN RDKIT CPP FUNCTION getNonQueryBondSmarts
    // RDKit✔️✔️: std::string getNonQueryBondSmarts(const Bond *qbond, int atomToLeftIdx,
    // RDKit✔️✔️:                                   const SmilesWriteParams &params) {
    // RDKit✔️✔️:   PRECONDITION(qbond, "bad bond");
    // RDKit✔️✔️:   std::string res;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (qbond->getIsAromatic()) {
    // RDKit✔️✔️:     res = ":";
    // RDKit✔️✔️:     if (params.doIsomericSmiles) {
    // RDKit✔️✔️:       if (qbond->getBondDir() == Bond::ENDDOWNRIGHT) {
    // RDKit✔️✔️:         res = "\\";
    // RDKit✔️✔️:       } else if (qbond->getBondDir() == Bond::ENDUPRIGHT) {
    // RDKit✔️✔️:         res = "/";
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     bool reverseDative =
    // RDKit✔️✔️:         (atomToLeftIdx >= 0 &&
    // RDKit✔️✔️:          qbond->getBeginAtomIdx() != static_cast<unsigned int>(atomToLeftIdx));
    // RDKit✔️✔️:     res = getBasicBondRepr(qbond->getBondType(), qbond->getBondDir(),
    // RDKit✔️✔️:                            reverseDative, params);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION getNonQueryBondSmarts
    // Complexity review: both use constant-time field reads and one bounded
    // bond-symbol dispatch. Rust delegates non-aromatic mapping to the sole
    // get_basic_bond_repr implementation, with no extra traversal or clone.
    if bond.is_aromatic() {
        if params.do_isomeric_smiles {
            match bond.direction() {
                BondDirection::EndDownRight => return "\\".to_owned(),
                BondDirection::EndUpRight => return "/".to_owned(),
                _ => {}
            }
        }
        return ":".to_owned();
    }
    let reverse_dative = atom_to_left_idx.is_some_and(|atom_idx| bond.begin().index() != atom_idx);
    get_basic_bond_repr(bond.order(), bond.direction(), reverse_dative, params)
}

fn mol_to_smarts_internal<FA, FB>(
    input_molecule: &Molecule,
    params: &SmilesWriteParams,
    atoms_in_play: &[AtomId],
    bonds_in_play: Option<&[BondId]>,
    mut write_atom: FA,
    mut write_bond: FB,
) -> Result<SmartsWriteResult, SmartsWriteError>
where
    FA: FnMut(&Atom, &SmilesWriteParams) -> Result<String, SmartsWriteError>,
    FB: FnMut(&Bond, &SmilesWriteParams, AtomId) -> Result<String, SmartsWriteError>,
{
    // RDKit✔️❌: std::string molToSmarts(const ROMol &inmol, const SmilesWriteParams &params,
    // RDKit✔️❌:                         std::vector<AtomColors> &&colors,
    // RDKit✔️❌:                         const boost::dynamic_bitset<> &atomsInPlay,
    // RDKit✔️❌:                         const boost::dynamic_bitset<> *bondsInPlay) {
    // RDKit✔️❌:   PRECONDITION(params.rootedAtAtom < static_cast<int>(inmol.getNumAtoms()),
    // RDKit✔️❌:                "bad atom index");
    // RDKit✔️❌:   ROMol mol(inmol);
    // RDKit✔️❌:   const unsigned int nAtoms = mol.getNumAtoms();
    // RDKit✔️❌:   UINT_VECT ranks;
    // RDKit✔️❌:   ranks.reserve(nAtoms);
    // RDKit✔️❌:   // For smiles writing we would be canonicalizing but we will not do that
    // RDKit✔️❌:   // here. We will simply use the atom indices as the rank
    // RDKit✔️❌:   for (const auto &atom : mol.atoms()) {
    // RDKit✔️❌:     ranks.push_back(atom->getIdx());
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   if (params.doIsomericSmiles) {
    // RDKit✔️❌:     mol.setProp(common_properties::_doIsoSmiles, 1);
    // RDKit✔️❌:   }
    // RDKit✔️❌:   std::vector<unsigned int> atomOrdering;
    // RDKit✔️❌:   std::vector<unsigned int> bondOrdering;
    // RDKit✔️❌:
    // RDKit✔️❌:   std::string res;
    // RDKit✔️❌:   auto colorIt = std::find(colors.begin(), colors.end(), Canon::WHITE_NODE);
    // RDKit✔️❌:   while (colorIt != colors.end()) {
    // RDKit✔️❌:     unsigned int nextAtomIdx = 0;
    // RDKit✔️❌:     std::string subSmi;
    // RDKit✔️❌:
    // RDKit✔️❌:     if (params.rootedAtAtom > -1 &&
    // RDKit✔️❌:         colors[params.rootedAtAtom] == Canon::WHITE_NODE) {
    // RDKit✔️❌:       nextAtomIdx = params.rootedAtAtom;
    // RDKit✔️❌:     } else {
    // RDKit✔️❌:       // Try to find a non-chiral atom we have not processed yet.
    // RDKit✔️❌:       // If we can't find non-chiral atom, use the chiral atom with
    // RDKit✔️❌:       // the lowest rank (we are guaranteed to find an unprocessed atom).
    // RDKit✔️❌:       unsigned nextRank = nAtoms + 1;
    // RDKit✔️❌:       for (auto atom : mol.atoms()) {
    // RDKit✔️❌:         if (colors[atom->getIdx()] == Canon::WHITE_NODE) {
    // RDKit✔️❌:           if (atom->getChiralTag() != Atom::CHI_TETRAHEDRAL_CCW &&
    // RDKit✔️❌:               atom->getChiralTag() != Atom::CHI_TETRAHEDRAL_CW) {
    // RDKit✔️❌:             nextAtomIdx = atom->getIdx();
    // RDKit✔️❌:             break;
    // RDKit✔️❌:           }
    // RDKit✔️❌:           if (ranks[atom->getIdx()] < nextRank) {
    // RDKit✔️❌:             nextRank = ranks[atom->getIdx()];
    // RDKit✔️❌:             nextAtomIdx = atom->getIdx();
    // RDKit✔️❌:           }
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:
    // RDKit✔️❌:     subSmi = FragmentSmartsConstruct(mol, nextAtomIdx, colors, ranks, params,
    // RDKit✔️❌:                                      atomOrdering, bondOrdering, atomsInPlay,
    // RDKit✔️❌:                                      bondsInPlay);
    // RDKit✔️❌:     res += subSmi;
    // RDKit✔️❌:
    // RDKit✔️❌:     colorIt = std::find(colors.begin(), colors.end(), Canon::WHITE_NODE);
    // RDKit✔️❌:     if (colorIt != colors.end()) {
    // RDKit✔️❌:       res += ".";
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   inmol.setProp(common_properties::_smilesAtomOutputOrder, atomOrdering, true);
    // RDKit✔️❌:   inmol.setProp(common_properties::_smilesBondOutputOrder, bondOrdering, true);
    // RDKit✔️❌:   return res;
    // RDKit✔️❌: }
    // Complexity review: the output and selection behavior match RDKit, but
    // BTreeSet-backed scope/component membership adds a log-factor versus the
    // source bitsets. Ordering is returned as typed state instead of encoded.
    if let Some(root) = params.rooted_at_atom
        && root >= input_molecule.num_atoms()
    {
        return Err(SmartsWriteError::RootedAtomOutOfRange { atom: root });
    }
    let mut molecule = input_molecule.clone();
    let allowed_atoms = atoms_in_play.iter().copied().collect::<BTreeSet<_>>();
    let allowed_bonds = bonds_in_play.map_or_else(
        || {
            molecule
                .bonds()
                .iter()
                .filter(|bond| {
                    allowed_atoms.contains(&bond.begin()) && allowed_atoms.contains(&bond.end())
                })
                .map(Bond::id)
                .collect::<BTreeSet<_>>()
        },
        |bonds| bonds.iter().copied().collect::<BTreeSet<_>>(),
    );
    let mut remaining = allowed_atoms.clone();
    let mut result = SmartsWriteResult::default();
    while !remaining.is_empty() {
        let start = params
            .rooted_at_atom
            .map(AtomId::new)
            .filter(|atom| remaining.contains(atom))
            .or_else(|| {
                molecule
                    .atoms()
                    .iter()
                    .find(|atom| {
                        remaining.contains(&atom.id())
                            && !matches!(
                                atom.chiral_tag(),
                                ChiralTag::TetrahedralCcw | ChiralTag::TetrahedralCw
                            )
                    })
                    .map(Atom::id)
            })
            .unwrap_or_else(|| *remaining.iter().next().expect("remaining is non-empty"));
        let mut component_atoms = Vec::new();
        let mut component_bonds = BTreeSet::new();
        let mut pending = vec![start];
        let mut seen = BTreeSet::new();
        while let Some(atom) = pending.pop() {
            if !remaining.contains(&atom) || !seen.insert(atom) {
                continue;
            }
            component_atoms.push(atom);
            for neighbor in molecule
                .topology_block()
                .adjacency
                .neighbors_of(atom.index())
            {
                if allowed_bonds.contains(&neighbor.bond) {
                    let other = AtomId::new(neighbor.atom_index);
                    if allowed_atoms.contains(&other) {
                        component_bonds.insert(neighbor.bond);
                        pending.push(other);
                    }
                }
            }
        }
        component_atoms.sort_by_key(|atom| atom.index());
        let component_bonds = component_bonds.into_iter().collect::<Vec<_>>();
        let ranks = component_atoms
            .iter()
            .map(|atom| atom.index())
            .collect::<Vec<_>>();
        result.smarts.push_str(&fragment_smarts_construct(
            &mut molecule,
            start,
            &ranks,
            params,
            &mut result.atom_ordering,
            &mut result.bond_ordering,
            &component_atoms,
            &component_bonds,
            &mut write_atom,
            &mut write_bond,
        )?);
        for atom in component_atoms {
            remaining.remove(&atom);
        }
        if !remaining.is_empty() {
            result.smarts.push('.');
        }
    }
    Ok(result)
}

pub fn get_atom_smarts(
    atom: &Atom,
    params: &SmilesWriteParams,
) -> Result<String, SmartsWriteError> {
    // RDKit✔️✔️: std::string GetAtomSmarts(const Atom *atom, const SmilesWriteParams &params) {
    // RDKit✔️✔️:   PRECONDITION(atom, "bad atom");
    // RDKit✔️✔️:   std::string res;
    // RDKit✔️✔️:   bool needParen = false;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // BOOST_LOG(rdInfoLog)<<"Atom: " <<qatom->getIdx()<<std::endl;
    // RDKit✔️✔️:   if (!atom->hasQuery()) {
    // RDKit✔️✔️:     res = getNonQueryAtomSmarts(atom);
    // RDKit✔️✔️:     // BOOST_LOG(rdInfoLog)<<"\tno query:"<<res;
    // RDKit✔️✔️:     return res;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   const auto query = atom->getQuery();
    // RDKit✔️✔️:   PRECONDITION(query, "atom has no query");
    // RDKit✔️✔️:   unsigned int queryFeatures = 0;
    // RDKit✔️✔️:   std::string descrip = query->getDescription();
    // RDKit✔️✔️:   if (descrip.empty()) {
    // RDKit✔️✔️:     // we have simple atom - just generate the smiles and return
    // RDKit✔️✔️:     res = SmilesWrite::GetAtomSmiles(atom);
    // RDKit✔️✔️:     return res;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     if ((descrip == "AtomOr") || (descrip == "AtomAnd")) {
    // RDKit✔️✔️:       const QueryAtom *qatom = dynamic_cast<const QueryAtom *>(atom);
    // RDKit✔️✔️:       PRECONDITION(qatom, "could not convert atom to query atom");
    // RDKit✔️✔️:       // we have a composite query
    // RDKit✔️✔️:       needParen = true;
    // RDKit✔️✔️:       res = _recurseGetSmarts(qatom, query, query->getNegation(), queryFeatures,
    // RDKit✔️✔️:                               params);
    // RDKit✔️✔️:       if (res.length() == 1) {  // single atom symbol we don't need parens
    // RDKit✔️✔️:         needParen = false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (descrip == "RecursiveStructure") {
    // RDKit✔️✔️:       // it's a bare recursive structure query:
    // RDKit✔️✔️:       res = getRecursiveStructureQuerySmarts(query, params);
    // RDKit✔️✔️:       needParen = true;
    // RDKit✔️✔️:     } else {  // we have a simple smarts
    // RDKit✔️✔️:       const QueryAtom *qatom = dynamic_cast<const QueryAtom *>(atom);
    // RDKit✔️✔️:       PRECONDITION(qatom, "could not convert atom to query atom");
    // RDKit✔️✔️:       res = getAtomSmartsSimple(qatom, query, needParen, true, params);
    // RDKit✔️✔️:       if (query->getNegation()) {
    // RDKit✔️✔️:         res = "!" + res;
    // RDKit✔️✔️:         needParen = true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     std::string mapNum;
    // RDKit✔️✔️:     if (atom->getPropIfPresent(common_properties::molAtomMapNumber, mapNum)) {
    // RDKit✔️✔️:       needParen = true;
    // RDKit✔️✔️:       res += ":" + mapNum;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     std::string symbol;
    // RDKit✔️✔️:     if (atom->getPropIfPresent(common_properties::smilesSymbol, symbol)) {
    // RDKit✔️✔️:       needParen = true;
    // RDKit✔️✔️:       if (!res.empty()) {
    // RDKit✔️✔️:         res = symbol + ";" + res;
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         res = symbol;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (needParen) {
    // RDKit✔️✔️:       res = "[" + res + "]";
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return res;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // Complexity review: both dispatch once over the typed root and visit the
    // query tree once. Rust uses caller-owned stereo state and monomorphized
    // recursive writing, with no query clone, reparse, or dynamic dispatch.
    Ok(get_non_query_atom_smarts(
        atom,
        params.do_isomeric_smiles,
        &mut false,
    ))
}

pub fn get_bond_smarts(
    bond: &Bond,
    params: &SmilesWriteParams,
    atom_to_left_idx: Option<usize>,
) -> Result<String, SmartsWriteError> {
    // RDKit✔️✔️: std::string GetBondSmarts(const Bond *bond, const SmilesWriteParams &params,
    // RDKit✔️✔️:                           int atomToLeftIdx) {
    // RDKit✔️✔️:   PRECONDITION(bond, "bad bond");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::string res = "";
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // BOOST_LOG(rdInfoLog) << "bond: " << bond->getIdx() << std::endl;
    // RDKit✔️✔️:   ;
    // RDKit✔️✔️:   // it is possible that we are regular single bond and we don't need to write
    // RDKit✔️✔️:   // anything
    // RDKit✔️✔️:   if (!bond->hasQuery()) {
    // RDKit✔️✔️:     res = getNonQueryBondSmarts(bond, atomToLeftIdx, params);
    // RDKit✔️✔️:     // BOOST_LOG(rdInfoLog) << "\tno query:" << res << std::endl;
    // RDKit✔️✔️:     return res;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // describeQuery(bond->getQuery());
    // RDKit✔️✔️:   auto qbond = dynamic_cast<const QueryBond *>(bond);
    // RDKit✔️✔️:   if (!qbond && ((bond->getBondType() == Bond::SINGLE) ||
    // RDKit✔️✔️:                  (bond->getBondType() == Bond::AROMATIC))) {
    // RDKit✔️✔️:     BOOST_LOG(rdInfoLog) << "\tbasic:" << res << std::endl;
    // RDKit✔️✔️:     return res;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   CHECK_INVARIANT(qbond, "could not convert bond to QueryBond");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   const auto query = qbond->getQuery();
    // RDKit✔️✔️:   CHECK_INVARIANT(query, "bond has no query");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   unsigned int queryFeatures = 0;
    // RDKit✔️✔️:   auto descrip = query->getDescription();
    // RDKit✔️✔️:   if ((descrip == "BondAnd") || (descrip == "BondOr")) {
    // RDKit✔️✔️:     // composite query
    // RDKit✔️✔️:     res = _recurseBondSmarts(bond, query, query->getNegation(), atomToLeftIdx,
    // RDKit✔️✔️:                              queryFeatures, params);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     // simple query
    // RDKit✔️✔️:     if (query->getNegation()) {
    // RDKit✔️✔️:       res = "!";
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     res += getBondSmartsSimple(bond, query, atomToLeftIdx, params);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // BOOST_LOG(rdInfoLog) << "\t  query:" << descrip << " " << res << std::endl;
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // Complexity review: both perform one root dispatch and visit each query
    // node once. Rust reuses the canonical leaf and recursive serializers and
    // uses typed query state, with no query clone, string reparse, or RTTI.
    Ok(get_non_query_bond_smarts(bond, atom_to_left_idx, params))
}

pub fn mol_to_smarts(
    molecule: &Molecule,
    params: &SmilesWriteParams,
) -> Result<String, SmartsWriteError> {
    // RDKit✔️✔️: std::string MolToSmarts(const ROMol &mol, const SmilesWriteParams &ps) {
    // RDKit✔️✔️:   const unsigned int nAtoms = mol.getNumAtoms();
    // RDKit✔️✔️:   if (!nAtoms) {
    // RDKit✔️✔️:     return "";
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::vector<AtomColors> colors(nAtoms, Canon::WHITE_NODE);
    // RDKit✔️✔️:   boost::dynamic_bitset<> atomsInPlay(nAtoms);
    // RDKit✔️✔️:   atomsInPlay.set();  // all atoms are in play
    // RDKit✔️✔️:   return molToSmarts(mol, ps, std::move(colors), atomsInPlay, nullptr);
    // RDKit✔️✔️: }
    // Complexity review: both initialize one O(n) all-atoms scope and invoke
    // the sole molecule traversal once. Rust carries typed AtomIds in a Vec;
    // atom and bond callbacks dispatch directly to the canonical writers.
    if molecule.num_atoms() == 0 {
        return Ok(String::new());
    }
    let atoms_in_play = molecule.atoms().iter().map(Atom::id).collect::<Vec<_>>();
    mol_to_smarts_internal(
        molecule,
        params,
        &atoms_in_play,
        None,
        |atom, atom_params| get_atom_smarts(atom, atom_params),
        |bond, bond_params, atom_to_left| {
            get_bond_smarts(bond, bond_params, Some(atom_to_left.index()))
        },
    )
    .map(|result| result.smarts)
}

pub fn mol_fragment_to_smarts(
    molecule: &Molecule,
    params: &SmilesWriteParams,
    atoms_to_use: &[AtomId],
    bonds_to_use: Option<&[BondId]>,
) -> Result<String, SmartsWriteError> {
    // RDKit✔️✔️: std::string MolFragmentToSmarts(const ROMol &mol,
    // RDKit✔️✔️:                                 const SmilesWriteParams &params,
    // RDKit✔️✔️:                                 const std::vector<int> &atomsToUse,
    // RDKit✔️✔️:                                 const std::vector<int> *bondsToUse) {
    // RDKit✔️✔️:   PRECONDITION(!atomsToUse.empty(), "no atoms provided");
    // RDKit✔️✔️:   PRECONDITION(!bondsToUse || !bondsToUse->empty(), "no bonds provided");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto nAtoms = mol.getNumAtoms();
    // RDKit✔️✔️:   if (!nAtoms) {
    // RDKit✔️✔️:     return "";
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::unique_ptr<boost::dynamic_bitset<>> bondsInPlay(nullptr);
    // RDKit✔️✔️:   if (bondsToUse != nullptr) {
    // RDKit✔️✔️:     bondsInPlay.reset(new boost::dynamic_bitset<>(mol.getNumBonds(), 0));
    // RDKit✔️✔️:     for (auto bidx : *bondsToUse) {
    // RDKit✔️✔️:       bondsInPlay->set(bidx);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // Mark all atoms except the ones in atomIndices as already processed.
    // RDKit✔️✔️:   // white: unprocessed
    // RDKit✔️✔️:   // grey: partial
    // RDKit✔️✔️:   // black: complete
    // RDKit✔️✔️:   boost::dynamic_bitset<> atomsInPlay(nAtoms);
    // RDKit✔️✔️:   std::vector<AtomColors> colors(nAtoms, Canon::BLACK_NODE);
    // RDKit✔️✔️:   for (const auto &idx : atomsToUse) {
    // RDKit✔️✔️:     colors[idx] = Canon::WHITE_NODE;
    // RDKit✔️✔️:     atomsInPlay.set(idx);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   SmilesWriteParams ps(params);
    // RDKit✔️✔️:   ps.rootedAtAtom = -1;
    // RDKit✔️✔️:   return molToSmarts(mol, ps, std::move(colors), atomsInPlay,
    // RDKit✔️✔️:                      bondsInPlay.get());
    // RDKit✔️✔️: }
    // Complexity review: both validate and materialize atom/bond selections
    // in linear time, then call the sole scoped traversal. Typed IDs avoid
    // conversion maps; bounds checks are constant-time per selected row.
    if atoms_to_use.is_empty() {
        return Err(SmartsWriteError::EmptyAtomSelection);
    }
    if bonds_to_use.is_some_and(<[BondId]>::is_empty) {
        return Err(SmartsWriteError::EmptyBondSelection);
    }
    for atom in atoms_to_use {
        if atom.index() >= molecule.num_atoms() {
            return Err(SmartsWriteError::FragmentAtomOutOfRange { atom: atom.index() });
        }
    }
    if let Some(bonds) = bonds_to_use {
        for bond in bonds {
            if bond.index() >= molecule.num_bonds() {
                return Err(SmartsWriteError::FragmentBondOutOfRange { bond: bond.index() });
            }
        }
    }
    let mut fragment_params = params.clone();
    fragment_params.rooted_at_atom = None;
    mol_to_smarts_internal(
        molecule,
        &fragment_params,
        atoms_to_use,
        bonds_to_use,
        |atom, atom_params| get_atom_smarts(atom, atom_params),
        |bond, bond_params, atom_to_left| {
            get_bond_smarts(bond, bond_params, Some(atom_to_left.index()))
        },
    )
    .map(|result| result.smarts)
}

pub fn mol_to_cx_smarts(
    molecule: &Molecule,
    params: &SmilesWriteParams,
) -> Result<String, SmartsWriteError> {
    // RDKit✔️✔️: std::string MolToCXSmarts(const ROMol &mol, const SmilesWriteParams &params) {
    // RDKit✔️✔️:   SmilesWriteParams ps(params);
    // RDKit✔️✔️:   ps.includeDativeBonds = false;
    // RDKit✔️✔️:   auto res = MolToSmarts(mol, ps);
    // RDKit✔️✔️:   if (!res.empty()) {
    // RDKit✔️✔️:     auto cxext = SmilesWrite::getCXExtensions(mol);
    // RDKit✔️✔️:     if (!cxext.empty()) {
    // RDKit✔️✔️:       res += " " + cxext;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // Complexity review: both serialize SMARTS once and the existing CX
    // fields once, with output-linear concatenation. Rust reuses both sole
    // serializers and performs no query conversion or reparsing.
    let mut cx_params = params.clone();
    cx_params.include_dative_bonds = false;
    let mut result = mol_to_smarts(molecule, &cx_params)?;
    if !result.is_empty() {
        let extension =
            crate::notation::smiles_write::get_cx_extensions(molecule, crate::CxSmilesFields::ALL)
                .map_err(|source| SmartsWriteError::Traversal { source })?;
        if !extension.is_empty() {
            result.push(' ');
            result.push_str(&extension);
        }
    }
    Ok(result)
}

pub fn mol_fragment_to_cx_smarts(
    molecule: &Molecule,
    params: &SmilesWriteParams,
    atoms_to_use: &[AtomId],
    bonds_to_use: Option<&[BondId]>,
) -> Result<String, SmartsWriteError> {
    // RDKit✔️✔️: std::string MolFragmentToCXSmarts(const ROMol &mol,
    // RDKit✔️✔️:                                   const SmilesWriteParams &params,
    // RDKit✔️✔️:                                   const std::vector<int> &atomsToUse,
    // RDKit✔️✔️:                                   const std::vector<int> *bondsToUse) {
    // RDKit✔️✔️:   auto res = MolFragmentToSmarts(mol, params, atomsToUse, bondsToUse);
    // RDKit✔️✔️:   if (!res.empty()) {
    // RDKit✔️✔️:     auto cxext = SmilesWrite::getCXExtensions(mol);
    // RDKit✔️✔️:     if (!cxext.empty()) {
    // RDKit✔️✔️:       res += " " + cxext;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // Complexity review: both invoke the canonical fragment writer once and
    // the existing full-molecule CX serializer once, then concatenate in time
    // linear in the output. No fragment decoder or second CX core is added.
    let mut result = mol_fragment_to_smarts(molecule, params, atoms_to_use, bonds_to_use)?;
    if !result.is_empty() {
        let extension =
            crate::notation::smiles_write::get_cx_extensions(molecule, crate::CxSmilesFields::ALL)
                .map_err(|source| SmartsWriteError::Traversal { source })?;
        if !extension.is_empty() {
            result.push(' ');
            result.push_str(&extension);
        }
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use crate::{
        AtomId, AtomQueryPredicate, AtomSpec, Bond, BondId, BondOrder, BondSpec, Element,
        Hybridization, MoleculeBuilder, QueryNode, SmilesWriteParams,
    };

    use super::super::query::RecursiveStructureQuery;

    use super::{
        QueryBoolFeatures, SmartsWriteError, combine_child_smarts, describe_query,
        fragment_smarts_construct, get_atom_smarts, get_atom_smarts_simple, get_basic_bond_repr,
        get_bond_smarts, get_bond_smarts_simple, get_non_query_atom_smarts,
        get_non_query_bond_smarts, get_recursive_structure_query_smarts, mol_fragment_to_cx_smarts,
        mol_fragment_to_smarts, mol_to_cx_smarts, mol_to_smarts, mol_to_smarts_internal,
        query_graph_fragment_to_smarts, query_graph_to_smarts, recurse_bond_smarts,
        recurse_get_smarts,
    };

    #[test]
    fn smarts_write_combine_child_smarts() {
        let mut features = QueryBoolFeatures::default();
        assert_eq!(
            combine_child_smarts(
                "#6".into(),
                QueryBoolFeatures::default(),
                "#7".into(),
                QueryBoolFeatures::default(),
                "AtomOr",
                &mut features,
            )
            .unwrap(),
            "#6,#7"
        );
        assert!(features.contains(QueryBoolFeatures::HAS_OR));

        let mut low_and = QueryBoolFeatures::default();
        assert_eq!(
            combine_child_smarts(
                "#6,#7".into(),
                features,
                "R".into(),
                QueryBoolFeatures::default(),
                "AtomAnd",
                &mut low_and,
            )
            .unwrap(),
            "#6,#7;R"
        );
        assert!(low_and.contains(QueryBoolFeatures::HAS_LOW_AND));

        let error = combine_child_smarts(
            "#6,#7;R".into(),
            low_and,
            "#8".into(),
            QueryBoolFeatures::default(),
            "AtomOr",
            &mut QueryBoolFeatures::default(),
        )
        .unwrap_err();
        assert_eq!(error, SmartsWriteError::OrAboveAndBelowAnd);

        assert_eq!(
            combine_child_smarts(
                String::new(),
                QueryBoolFeatures::default(),
                "#8".into(),
                QueryBoolFeatures::default(),
                "AtomAnd",
                &mut QueryBoolFeatures::default(),
            )
            .unwrap(),
            "#8"
        );
    }

    #[test]
    fn smarts_write_describe_query() {
        let query = QueryNode::and(vec![
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
            QueryNode::or(vec![
                QueryNode::predicate(AtomQueryPredicate::FormalCharge(0)),
                QueryNode::not(QueryNode::predicate(AtomQueryPredicate::InRing)),
            ]),
        ]);
        describe_query(&query, "\t".to_owned());
    }

    #[test]
    fn smarts_write_get_atom_smarts_simple() {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(Element::C));
        let molecule = builder.build().unwrap();
        let atom = &molecule.atoms()[0];

        let mut need_paren = false;
        let mut stereo_written = false;
        assert_eq!(
            get_atom_smarts_simple(
                atom,
                &AtomQueryPredicate::AtomicNumber(6),
                &mut need_paren,
                false,
                true,
                &mut stereo_written,
            ),
            "#6"
        );
        assert!(need_paren);
        assert_eq!(
            get_atom_smarts_simple(
                atom,
                &AtomQueryPredicate::FormalCharge(-2),
                &mut need_paren,
                false,
                true,
                &mut stereo_written,
            ),
            "-2"
        );
        assert_eq!(
            get_atom_smarts_simple(
                atom,
                &AtomQueryPredicate::HybridizationMatch(Hybridization::Sp2),
                &mut need_paren,
                false,
                true,
                &mut stereo_written,
            ),
            "^2"
        );
        assert_eq!(
            get_atom_smarts_simple(
                atom,
                &AtomQueryPredicate::ExplicitDegreeLessEqual(3),
                &mut need_paren,
                false,
                true,
                &mut stereo_written,
            ),
            "D{3-}"
        );
    }

    #[test]
    fn smarts_write_get_recursive_structure_query_smarts() {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(Element::C));
        builder.add_atom(AtomSpec::new(Element::O));
        let molecule = builder.build().unwrap();
        let query = RecursiveStructureQuery::from_molecule(
            crate::search::query_graph::query_graph_from_concrete_molecule(&molecule).unwrap(),
            0,
        );
        let params = SmilesWriteParams::default();

        let write_nested = |molecule: &crate::QueryGraph, _: &SmilesWriteParams| {
            assert_eq!(molecule.num_atoms(), 2);
            Ok("C=O".to_owned())
        };
        assert_eq!(
            get_recursive_structure_query_smarts(&query, false, &params, write_nested).unwrap(),
            "$(C=O)"
        );
        assert_eq!(
            get_recursive_structure_query_smarts(&query, true, &params, |_, _| {
                Ok("C=O".to_owned())
            })
            .unwrap(),
            "!$(C=O)"
        );

        assert_eq!(
            get_recursive_structure_query_smarts(
                &RecursiveStructureQuery::new(),
                false,
                &params,
                |_, _| Ok(String::new()),
            )
            .unwrap_err(),
            super::SmartsWriteError::MissingRecursiveQueryMolecule
        );
    }

    #[test]
    fn smarts_write_get_basic_bond_repr() {
        let params = SmilesWriteParams::default();
        assert_eq!(
            get_basic_bond_repr(
                crate::BondOrder::Single,
                crate::BondDirection::EndDownRight,
                false,
                &params,
            ),
            "\\"
        );
        assert_eq!(
            get_basic_bond_repr(
                crate::BondOrder::Aromatic,
                crate::BondDirection::EndUpRight,
                false,
                &params,
            ),
            "/"
        );
        assert_eq!(
            get_basic_bond_repr(
                crate::BondOrder::Dative,
                crate::BondDirection::None,
                true,
                &params,
            ),
            "<-"
        );
        let no_dative = SmilesWriteParams {
            include_dative_bonds: false,
            ..params.clone()
        };
        assert_eq!(
            get_basic_bond_repr(
                crate::BondOrder::Dative,
                crate::BondDirection::None,
                false,
                &no_dative,
            ),
            "-"
        );
        assert_eq!(
            get_basic_bond_repr(
                crate::BondOrder::Zero,
                crate::BondDirection::None,
                false,
                &params,
            ),
            "~"
        );
        assert_eq!(
            get_basic_bond_repr(
                crate::BondOrder::Quintuple,
                crate::BondDirection::None,
                false,
                &params,
            ),
            ""
        );
    }

    #[test]
    fn smarts_write_get_bond_smarts_simple() {
        let mut builder = MoleculeBuilder::new();
        let carbon = builder.add_atom(AtomSpec::new(Element::C));
        let oxygen = builder.add_atom(AtomSpec::new(Element::O));
        builder
            .add_bond(
                crate::BondSpec::new(carbon, oxygen, crate::BondOrder::Single)
                    .with_direction(crate::BondDirection::EndDownRight),
            )
            .unwrap();
        let molecule = builder.build().unwrap();
        let bond = &molecule.bonds()[0];
        let params = SmilesWriteParams::default();

        assert_eq!(
            get_bond_smarts_simple(bond, &crate::BondQueryPredicate::Any, None, &params,).unwrap(),
            "~"
        );
        assert_eq!(
            get_bond_smarts_simple(
                bond,
                &crate::BondQueryPredicate::OrderIn(vec![
                    crate::BondOrder::Single,
                    crate::BondOrder::Aromatic,
                ]),
                None,
                &params,
            )
            .unwrap(),
            "\\"
        );
        assert_eq!(
            get_bond_smarts_simple(
                bond,
                &crate::BondQueryPredicate::OrderIn(vec![
                    crate::BondOrder::Single,
                    crate::BondOrder::Double,
                    crate::BondOrder::Aromatic,
                ]),
                None,
                &params,
            )
            .unwrap(),
            "-,=,:"
        );
        assert_eq!(
            get_bond_smarts_simple(
                bond,
                &crate::BondQueryPredicate::Order(crate::BondOrder::Dative),
                Some(oxygen.index()),
                &params,
            )
            .unwrap(),
            "<-"
        );
        assert!(matches!(
            get_bond_smarts_simple(
                bond,
                &crate::BondQueryPredicate::Direction(crate::BondDirection::BeginWedge),
                None,
                &params,
            ),
            Err(super::SmartsWriteError::UnsupportedBondDirection { .. })
        ));
    }

    #[test]
    fn smarts_write_recurse_get_smarts() {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(Element::C));
        let molecule = builder.build().unwrap();
        let atom = &molecule.atoms()[0];
        let params = SmilesWriteParams::default();

        let query = QueryNode::and(vec![
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
            QueryNode::or(vec![
                QueryNode::predicate(AtomQueryPredicate::HydrogenCount(1)),
                QueryNode::predicate(AtomQueryPredicate::HydrogenCount(2)),
            ]),
        ]);
        let mut features = QueryBoolFeatures::default();
        let mut stereo_written = false;
        assert_eq!(
            recurse_get_smarts(
                atom,
                &query,
                false,
                &mut features,
                &params,
                &mut stereo_written,
                &mut |_, _| Ok(String::new()),
            )
            .unwrap(),
            "#6;H1,H2"
        );
        assert!(features.contains(QueryBoolFeatures::HAS_LOW_AND));
        assert!(features.contains(QueryBoolFeatures::HAS_OR));

        let negated = QueryNode::not(QueryNode::or(vec![
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(7)),
        ]));
        assert_eq!(
            recurse_get_smarts(
                atom,
                &negated,
                false,
                &mut QueryBoolFeatures::default(),
                &params,
                &mut false,
                &mut |_, _| Ok(String::new()),
            )
            .unwrap(),
            "!#6&!#7"
        );

        let mut recursive_builder = MoleculeBuilder::new();
        recursive_builder.add_atom(AtomSpec::new(Element::N));
        let recursive_molecule = recursive_builder.build().unwrap();
        let recursive = RecursiveStructureQuery::from_molecule(
            crate::search::query_graph::query_graph_from_concrete_molecule(&recursive_molecule)
                .unwrap(),
            0,
        );
        let recursive_query = QueryNode::and(vec![
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
            QueryNode::predicate(AtomQueryPredicate::RecursiveSmarts(recursive)),
        ]);
        let mut recursive_features = QueryBoolFeatures::default();
        assert_eq!(
            recurse_get_smarts(
                atom,
                &recursive_query,
                false,
                &mut recursive_features,
                &params,
                &mut false,
                &mut |_, _| Ok("N".to_owned()),
            )
            .unwrap(),
            "#6&$(N)"
        );
        assert!(recursive_features.contains(QueryBoolFeatures::HAS_RECURSION));
    }

    #[test]
    fn smarts_write_recurse_bond_smarts() {
        let mut builder = MoleculeBuilder::new();
        let carbon = builder.add_atom(AtomSpec::new(Element::C));
        let oxygen = builder.add_atom(AtomSpec::new(Element::O));
        builder
            .add_bond(crate::BondSpec::new(
                carbon,
                oxygen,
                crate::BondOrder::Single,
            ))
            .unwrap();
        let molecule = builder.build().unwrap();
        let bond = &molecule.bonds()[0];
        let params = SmilesWriteParams::default();

        let query = QueryNode::and(vec![
            QueryNode::predicate(crate::BondQueryPredicate::Order(crate::BondOrder::Single)),
            QueryNode::predicate(crate::BondQueryPredicate::IsInRing(true)),
        ]);
        assert_eq!(
            recurse_bond_smarts(
                bond,
                &query,
                false,
                None,
                &mut QueryBoolFeatures::default(),
                &params,
            )
            .unwrap(),
            "-&@"
        );

        let negated = QueryNode::not(QueryNode::or(vec![
            QueryNode::predicate(crate::BondQueryPredicate::Order(crate::BondOrder::Single)),
            QueryNode::predicate(crate::BondQueryPredicate::IsInRing(true)),
        ]));
        assert_eq!(
            recurse_bond_smarts(
                bond,
                &negated,
                false,
                None,
                &mut QueryBoolFeatures::default(),
                &params,
            )
            .unwrap(),
            "!-&!@"
        );

        let second_composite = QueryNode::and(vec![
            QueryNode::predicate(crate::BondQueryPredicate::Order(crate::BondOrder::Single)),
            QueryNode::or(vec![
                QueryNode::predicate(crate::BondQueryPredicate::Order(crate::BondOrder::Double)),
                QueryNode::predicate(crate::BondQueryPredicate::Order(crate::BondOrder::Triple)),
            ]),
        ]);
        assert_eq!(
            recurse_bond_smarts(
                bond,
                &second_composite,
                false,
                None,
                &mut QueryBoolFeatures::default(),
                &params,
            )
            .unwrap(),
            "=,#"
        );
    }

    #[test]
    fn smarts_write_fragment_smarts_construct() {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(AtomSpec::new(Element::C));
        let a1 = builder.add_atom(AtomSpec::new(Element::N));
        let a2 = builder.add_atom(AtomSpec::new(Element::O));
        builder
            .add_bond(BondSpec::new(a0, a1, BondOrder::Single))
            .unwrap();
        builder
            .add_bond(BondSpec::new(a0, a2, BondOrder::Double))
            .unwrap();
        let mut molecule = builder.build().unwrap();
        let atoms = vec![a0, a1, a2];
        let bonds = molecule.bonds().iter().map(Bond::id).collect::<Vec<_>>();
        let mut atom_ordering = Vec::new();
        let mut bond_ordering = Vec::new();
        let written = fragment_smarts_construct(
            &mut molecule,
            AtomId::new(0),
            &[0, 1, 2],
            &SmilesWriteParams::default(),
            &mut atom_ordering,
            &mut bond_ordering,
            &atoms,
            &bonds,
            |atom, _| Ok(format!("[#{}]", atom.atomic_number())),
            |bond, _, _| {
                Ok(match bond.order() {
                    BondOrder::Single => "-".to_owned(),
                    BondOrder::Double => "=".to_owned(),
                    _ => "~".to_owned(),
                })
            },
        )
        .unwrap();

        assert_eq!(written, "[#6](-[#7])=[#8]");
        assert_eq!(atom_ordering, vec![a0, a1, a2]);
        assert_eq!(bond_ordering.len(), 2);
        assert!(molecule.derived_cache().rings.is_none());
    }

    #[test]
    fn smarts_write_get_non_query_atom_smarts() {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(
            AtomSpec::new(Element::C)
                .with_isotope(13)
                .with_chiral_tag(crate::ChiralTag::TetrahedralCw)
                .with_explicit_hydrogens(1)
                .with_formal_charge(2)
                .with_atom_map(7),
        );
        builder.add_atom(AtomSpec::new(Element::FE).with_formal_charge(-2));
        let molecule = builder.build().unwrap();

        let mut stereo_written = false;
        assert_eq!(
            get_non_query_atom_smarts(&molecule.atoms()[0], true, &mut stereo_written),
            "[13#6@@H+2:7]"
        );
        assert!(stereo_written);
        let mut no_stereo = false;
        assert_eq!(
            get_non_query_atom_smarts(&molecule.atoms()[1], true, &mut no_stereo),
            "[Fe-2]"
        );

        let mut custom = molecule.atoms()[0].clone();
        custom.set_prop("smilesSymbol", "Q");
        let mut custom_stereo = false;
        assert_eq!(
            get_non_query_atom_smarts(&custom, false, &mut custom_stereo),
            "[13Q+2:7]"
        );
    }

    #[test]
    fn smarts_write_get_non_query_bond_smarts() {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(AtomSpec::new(Element::C));
        let a1 = builder.add_atom(AtomSpec::new(Element::C));
        builder
            .add_bond(
                BondSpec::new(a0, a1, BondOrder::Aromatic)
                    .with_aromatic(true)
                    .with_direction(crate::BondDirection::EndDownRight),
            )
            .unwrap();
        let aromatic = builder.build().unwrap();
        let mut params = SmilesWriteParams::default();
        assert_eq!(
            get_non_query_bond_smarts(&aromatic.bonds()[0], Some(0), &params),
            "\\"
        );
        params.do_isomeric_smiles = false;
        assert_eq!(
            get_non_query_bond_smarts(&aromatic.bonds()[0], Some(0), &params),
            ":"
        );

        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(AtomSpec::new(Element::N));
        let a1 = builder.add_atom(AtomSpec::new(Element::B));
        builder
            .add_bond(BondSpec::new(a0, a1, BondOrder::Dative))
            .unwrap();
        let dative = builder.build().unwrap();
        params.include_dative_bonds = true;
        assert_eq!(
            get_non_query_bond_smarts(&dative.bonds()[0], Some(a1.index()), &params),
            "<-"
        );
        assert_eq!(
            get_non_query_bond_smarts(&dative.bonds()[0], Some(a0.index()), &params),
            "->"
        );
    }

    #[test]
    fn smarts_write_get_atom_smarts() {
        let params = SmilesWriteParams::default();
        let query = crate::search::smarts_parse::parse_smarts(
            "[#6].[#6&H1].[!#8].[#7:12].[Q].[$([#7])]",
            &crate::search::smarts_parse::SmartsParseParams::default(),
        )
        .unwrap();
        assert_eq!(
            query_graph_to_smarts(&query, &params).unwrap(),
            "[#6].[#6&H1].[!#8].[#7:12].[Q].[$([#7])]"
        );
    }

    #[test]
    fn smarts_write_get_bond_smarts() {
        let query = crate::search::smarts_parse::parse_smarts(
            "C-C.C-&@N.C!=O",
            &crate::search::smarts_parse::SmartsParseParams::default(),
        )
        .unwrap();
        let params = SmilesWriteParams::default();
        assert_eq!(
            query_graph_to_smarts(&query, &params).unwrap(),
            "[#6]-[#6].[#6]-&@[#7].[#6]!=[#8]"
        );
    }

    #[test]
    fn smarts_write_mol_to_smarts() {
        let query = crate::search::smarts_parse::parse_smarts(
            "[#6]=[#8&$([#7])]",
            &crate::search::smarts_parse::SmartsParseParams::default(),
        )
        .unwrap();
        assert_eq!(
            query_graph_to_smarts(&query, &SmilesWriteParams::default()).unwrap(),
            "[#6]=[#8&$([#7])]"
        );
    }

    #[test]
    fn smarts_write_mol_fragment_to_smarts() {
        let query = crate::search::smarts_parse::parse_smarts(
            "[#6]-[#7]=[#8]",
            &crate::search::smarts_parse::SmartsParseParams::default(),
        )
        .unwrap();
        let mut params = SmilesWriteParams::default();
        params.rooted_at_atom = Some(2);
        assert_eq!(
            query_graph_fragment_to_smarts(
                &query,
                &params,
                &[AtomId::new(0), AtomId::new(1)],
                Some(&[BondId::new(0)]),
            )
            .unwrap(),
            "[#6]-[#7]"
        );
        assert_eq!(
            query_graph_fragment_to_smarts(
                &query,
                &params,
                &[AtomId::new(0), AtomId::new(2)],
                None
            )
            .unwrap(),
            "[#6].[#8]"
        );
        assert_eq!(
            query_graph_fragment_to_smarts(&query, &params, &[], None).unwrap_err(),
            SmartsWriteError::EmptyAtomSelection
        );
    }

    #[test]
    fn smarts_write_mol_to_c_x_smarts() {
        let empty = MoleculeBuilder::new().build().unwrap();
        assert_eq!(
            mol_to_cx_smarts(&empty, &SmilesWriteParams::default()).unwrap(),
            ""
        );

        let mut builder = MoleculeBuilder::new();
        builder.add_atom(AtomSpec::new(Element::C).with_prop("atomLabel", "site"));
        let molecule = builder.build().unwrap();
        assert_eq!(
            mol_to_cx_smarts(&molecule, &SmilesWriteParams::default()).unwrap(),
            "[#6] |$site$|"
        );
    }

    #[test]
    fn smarts_write_mol_fragment_to_c_x_smarts() {
        let mut builder = MoleculeBuilder::new();
        let carbon = builder.add_atom(AtomSpec::new(Element::C).with_prop("atomLabel", "left"));
        let oxygen = builder.add_atom(AtomSpec::new(Element::O).with_prop("atomLabel", "right"));
        builder
            .add_bond(BondSpec::new(carbon, oxygen, BondOrder::Double))
            .unwrap();
        let molecule = builder.build().unwrap();

        assert_eq!(
            mol_fragment_to_cx_smarts(&molecule, &SmilesWriteParams::default(), &[oxygen], None,)
                .unwrap(),
            "[#8] |$left;right$|"
        );
    }

    #[test]
    fn smarts_write_mol_to_smarts_internal() {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(AtomSpec::new(Element::C));
        let a1 = builder.add_atom(AtomSpec::new(Element::O));
        let a2 = builder.add_atom(AtomSpec::new(Element::N));
        builder
            .add_bond(BondSpec::new(a0, a1, BondOrder::Double))
            .unwrap();
        let molecule = builder.build().unwrap();
        let result = mol_to_smarts_internal(
            &molecule,
            &SmilesWriteParams::default(),
            &[a0, a1, a2],
            None,
            |atom, _| Ok(format!("[#{}]", atom.atomic_number())),
            |bond, _, _| {
                Ok(match bond.order() {
                    BondOrder::Double => "=".to_owned(),
                    _ => "-".to_owned(),
                })
            },
        )
        .unwrap();
        assert_eq!(result.smarts, "[#6]=[#8].[#7]");
        assert_eq!(result.atom_ordering, vec![a0, a1, a2]);
        assert_eq!(result.bond_ordering.len(), 1);

        let mut rooted = SmilesWriteParams::default();
        rooted.rooted_at_atom = Some(a1.index());
        let rooted_result = mol_to_smarts_internal(
            &molecule,
            &rooted,
            &[a0, a1],
            None,
            |atom, _| Ok(format!("[#{}]", atom.atomic_number())),
            |_, _, _| Ok("=".to_owned()),
        )
        .unwrap();
        assert_eq!(rooted_result.smarts, "[#8]=[#6]");
    }
}
