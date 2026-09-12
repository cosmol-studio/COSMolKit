//! RDKit SMARTS serialization for canonical [`QueryGraph`] values.

use cosmolkit_model::{
    Atom, AtomId, AtomQueryPredicate, AtomRangeBounds, AtomRangeDataFunction, Bond, BondId,
    BondQueryPredicate, QueryAtom, QueryBond, QueryGraph, QueryNode, RecursiveStructureQuery,
};
use cosmolkit_types::{BondDirection, BondOrder, ChiralTag, Element, Hybridization};
use std::collections::BTreeSet;

/// Options for QueryGraph-native SMARTS serialization.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SmartsWriteParams {
    /// Emit atom-map labels stored on query atoms.
    pub include_atom_maps: bool,
    /// Emit atom and bond stereochemistry.
    pub do_isomeric_smiles: bool,
    /// Preserve directional dative bond tokens.
    pub include_dative_bonds: bool,
    /// Start graph traversal at this atom index when present.
    pub rooted_at_atom: Option<usize>,
}

impl Default for SmartsWriteParams {
    fn default() -> Self {
        Self {
            include_atom_maps: true,
            do_isomeric_smiles: true,
            include_dative_bonds: true,
            rooted_at_atom: None,
        }
    }
}

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
    #[error("query graph is invalid: {0}")]
    InvalidGraph(String),
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
pub fn query_graph_to_smarts(
    query: &QueryGraph,
    params: &SmartsWriteParams,
) -> Result<String, SmartsWriteError> {
    query_graph_to_smarts_fragment(query, params, None, None, false)
}

pub fn query_graph_to_cx_smarts(
    query: &QueryGraph,
    params: &SmartsWriteParams,
) -> Result<String, SmartsWriteError> {
    query_graph_to_smarts_fragment(query, params, None, None, true)
}

pub fn query_graph_fragment_to_smarts(
    query: &QueryGraph,
    params: &SmartsWriteParams,
    atoms: &[AtomId],
    bonds: Option<&[BondId]>,
) -> Result<String, SmartsWriteError> {
    query_graph_to_smarts_fragment(query, params, Some(atoms), bonds, false)
}

pub fn query_graph_fragment_to_cx_smarts(
    query: &QueryGraph,
    params: &SmartsWriteParams,
    atoms: &[AtomId],
    bonds: Option<&[BondId]>,
) -> Result<String, SmartsWriteError> {
    query_graph_to_smarts_fragment(query, params, Some(atoms), bonds, true)
}

fn query_graph_to_smarts_fragment(
    query: &QueryGraph,
    params: &SmartsWriteParams,
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
            .any(|key| key != crate::query_graph_behavior::UNSPECIFIED_ORDER_PROP)
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
    params: &SmartsWriteParams,
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

#[doc(hidden)]
pub fn query_atom_to_smarts(
    atom: &QueryAtom,
    params: &SmartsWriteParams,
) -> Result<String, SmartsWriteError> {
    let mut features = QueryBoolFeatures::default();
    let mut stereo_written = false;
    let mut needs_brackets;
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
    if params.include_atom_maps
        && let Some(map) = atom.atom_map()
    {
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

#[doc(hidden)]
pub fn query_bond_to_smarts(
    bond: &QueryBond,
    params: &SmartsWriteParams,
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
    // RDKit✔️✔️:                                 const SmartsWriteParams &) {
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
                if in_organic_subset(*atomic_number).unwrap_or(false) {
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
    params: &SmartsWriteParams,
    write_molecule: F,
) -> Result<String, SmartsWriteError>
where
    F: FnOnce(&crate::QueryGraph, &SmartsWriteParams) -> Result<String, SmartsWriteError>,
{
    // BEGIN RDKIT CPP FUNCTION getRecursiveStructureQuerySmarts
    // RDKit✔️✔️: std::string getRecursiveStructureQuerySmarts(
    // RDKit✔️✔️:     const QueryAtom::QUERYATOM_QUERY *query, const SmartsWriteParams &params) {
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
        .query_graph()
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
    params: &SmartsWriteParams,
) -> String {
    // BEGIN RDKIT CPP FUNCTION getBasicBondRepr
    // RDKit✔️✔️: std::string getBasicBondRepr(Bond::BondType typ, Bond::BondDir dir,
    // RDKit✔️✔️:                              bool reverseDative,
    // RDKit✔️✔️:                              const SmartsWriteParams &params) {
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

fn in_organic_subset(atomic_number: u8) -> Result<bool, std::convert::Infallible> {
    // BEGIN RDKIT CPP FUNCTION inOrganicSubset
    // RDKit✔️✔️: const int atomicSmiles[] = {0, 5, 6, 7, 8, 9, 15, 16, 17, 35, 53, -1};
    // RDKit✔️✔️: bool inOrganicSubset(int atomicNumber) {
    // RDKit✔️✔️:   unsigned int idx = 0;
    // RDKit✔️✔️:   while (atomicSmiles[idx] < atomicNumber && atomicSmiles[idx] != -1) {
    // RDKit✔️✔️:     ++idx;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return atomicSmiles[idx] == atomicNumber;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION inOrganicSubset
    Ok(matches!(
        atomic_number,
        0 | 5 | 6 | 7 | 8 | 9 | 15 | 16 | 17 | 35 | 53
    ))
}

fn get_bond_smarts_simple(
    bond: &Bond,
    query: &BondQueryPredicate,
    atom_to_left_idx: Option<usize>,
    params: &SmartsWriteParams,
) -> Result<String, SmartsWriteError> {
    // BEGIN RDKIT CPP FUNCTION getBondSmartsSimple
    // RDKit✔️✔️: std::string getBondSmartsSimple(const Bond *bond,
    // RDKit✔️✔️:                                 const QueryBond::QUERYBOND_QUERY *bquery,
    // RDKit✔️✔️:                                 int atomToLeftIdx,
    // RDKit✔️✔️:                                 const SmartsWriteParams &params) {
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
    params: &SmartsWriteParams,
    stereo_written: &mut bool,
    write_molecule: &mut F,
) -> Result<String, SmartsWriteError>
where
    F: FnMut(&crate::QueryGraph, &SmartsWriteParams) -> Result<String, SmartsWriteError>,
{
    // BEGIN RDKIT CPP FUNCTION _recurseGetSmarts
    // RDKit✔️✔️: std::string _recurseGetSmarts(const QueryAtom *qatom,
    // RDKit✔️✔️:                               const QueryAtom::QUERYATOM_QUERY *node,
    // RDKit✔️✔️:                               bool negate, unsigned int &features,
    // RDKit✔️✔️:                               const SmartsWriteParams &params) {
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
    params: &SmartsWriteParams,
) -> Result<String, SmartsWriteError> {
    // BEGIN RDKIT CPP FUNCTION _recurseBondSmarts
    // RDKit✔️✔️: std::string _recurseBondSmarts(const Bond *bond,
    // RDKit✔️✔️:                                const QueryBond::QUERYBOND_QUERY *node,
    // RDKit✔️✔️:                                bool negate, int atomToLeftIdx,
    // RDKit✔️✔️:                                unsigned int &features,
    // RDKit✔️✔️:                                const SmartsWriteParams &params) {
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

/// Serialize a complete detached query graph, including cycles and recursive
/// SMARTS predicates.
pub fn write_smarts(
    graph: &QueryGraph,
    params: &SmartsWriteParams,
) -> Result<String, SmartsWriteError> {
    graph
        .validate()
        .map_err(|error| SmartsWriteError::InvalidGraph(error.to_string()))?;
    query_graph_to_smarts(graph, params)
}

/// Serialize one detached query atom.
pub fn atom_to_smarts(
    graph: &QueryGraph,
    atom_id: AtomId,
    params: &SmartsWriteParams,
) -> Result<String, SmartsWriteError> {
    graph
        .validate()
        .map_err(|error| SmartsWriteError::InvalidGraph(error.to_string()))?;
    let atom =
        graph
            .atom(atom_id.index())
            .ok_or_else(|| SmartsWriteError::FragmentAtomOutOfRange {
                atom: atom_id.index(),
            })?;
    query_atom_to_smarts(atom, params)
}

/// Serialize one detached query bond.
pub fn bond_to_smarts(graph: &QueryGraph, bond_id: BondId) -> Result<String, SmartsWriteError> {
    graph
        .validate()
        .map_err(|error| SmartsWriteError::InvalidGraph(error.to_string()))?;
    let bond =
        graph
            .bond(bond_id.index())
            .ok_or_else(|| SmartsWriteError::FragmentBondOutOfRange {
                bond: bond_id.index(),
            })?;
    query_bond_to_smarts(bond, &SmartsWriteParams::default(), None)
}
