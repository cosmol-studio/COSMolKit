//! RDKit SMARTS serialization for canonical [`QueryGraph`] values.

use cosmolkit_model::{
    AdjacencyList, AtomId, AtomQueryPredicate, AtomRangeBounds, AtomRangeDataFunction, Bond,
    BondId, BondQueryPredicate, QueryAtom, QueryBond, QueryGraph, QueryNode,
    RecursiveStructureQuery, SGroupConnection, StereoGroup, StereoGroupKind, SubstanceGroup,
    SubstanceGroupKind, query_substance_groups,
};
use cosmolkit_types::{BondDirection, BondOrder, BondStereo, ChiralTag, Element, Hybridization};
use std::collections::{BTreeMap, BTreeSet};

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
    #[error("Can't write smarts for this query atom type: {predicate:?}")]
    UnsupportedAtomQuery { predicate: AtomQueryPredicate },
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
    query_graph_to_smarts_fragment_result(query, params, atom_selection, bond_selection, include_cx)
        .map(|result| result.smarts)
}

fn query_graph_to_smarts_fragment_result(
    query: &QueryGraph,
    params: &SmartsWriteParams,
    atom_selection: Option<&[AtomId]>,
    bond_selection: Option<&[BondId]>,
    include_cx: bool,
) -> Result<SmartsWriteResult, SmartsWriteError> {
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
            result.smarts.push('.');
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
    if include_cx && !result.smarts.is_empty() {
        // RDKit✔️✔️:   if (!res.empty()) {
        // RDKit✔️✔️:     auto cxext = SmilesWrite::getCXExtensions(mol);
        // RDKit✔️✔️:     if (!cxext.empty()) {
        // RDKit✔️✔️:       res += " " + cxext;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        let extension =
            write_query_cx_extensions(query, &result.atom_ordering, &result.bond_ordering)?;
        if !extension.is_empty() {
            result.smarts.push(' ');
            result.smarts.push_str(&extension);
        }
    }
    Ok(result)
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

fn append_query_cx_extension(addition: String, output: &mut String) {
    // RDKit✔️✔️: void appendToCXExtension(const std::string &addition, std::string &base) {
    // RDKit✔️✔️:   if (!addition.empty()) {
    // RDKit✔️✔️:     if (base.size() > 1) {
    // RDKit✔️✔️:       base += ",";
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     base += addition;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    if !addition.is_empty() {
        if output.len() > 1 {
            output.push(',');
        }
        output.push_str(&addition);
    }
}

fn query_cx_atom_positions(atom_order: &[AtomId], atom_count: usize) -> Vec<Option<usize>> {
    let mut positions = vec![None; atom_count];
    for (position, atom) in atom_order.iter().copied().enumerate() {
        positions[atom.index()] = Some(position);
    }
    positions
}

fn query_cx_source_reverse_atom_order(atom_order: &[AtomId], atom_count: usize) -> Vec<usize> {
    // RDKit value-initializes the reverse vector. For fragment CXSMARTS this
    // intentionally leaves every unselected source atom mapped to output 0.
    let mut positions = vec![0; atom_count];
    for (position, atom) in atom_order.iter().copied().enumerate() {
        positions[atom.index()] = position;
    }
    positions
}

fn query_cx_zero_small(value: f64) -> f64 {
    if value.abs() < 1e-4 { 0.0 } else { value }
}

// `boost::format("%g")` uses six significant digits and switches to
// scientific notation outside the source's fixed-format exponent window.
fn query_cx_format_general(value: f64) -> String {
    if value == 0.0 {
        return "0".to_owned();
    }
    if !value.is_finite() {
        return value.to_string();
    }
    let exponent = value.abs().log10().floor() as i32;
    let scale = 10_f64.powi(5 - exponent);
    let rounded = (value * scale).round() / scale;
    let exponent = rounded.abs().log10().floor() as i32;
    if !(-4..6).contains(&exponent) {
        let mut text = format!("{rounded:.5e}");
        let exponent_at = text.find('e').expect("scientific format contains exponent");
        let exponent_text = text.split_off(exponent_at);
        while text.ends_with('0') {
            text.pop();
        }
        if text.ends_with('.') {
            text.pop();
        }
        let exponent_value = exponent_text[1..].parse::<i32>().unwrap_or(exponent);
        format!("{text}e{exponent_value:+03}")
    } else {
        let decimals = usize::try_from((5 - exponent).max(0)).unwrap_or(0);
        let mut text = format!("{rounded:.decimals$}");
        if text.contains('.') {
            while text.ends_with('0') {
                text.pop();
            }
            if text.ends_with('.') {
                text.pop();
            }
        }
        text
    }
}

fn write_query_cx_coordinates(query: &QueryGraph, atom_order: &[AtomId]) -> Option<String> {
    // RDKit✔️✔️: std::string get_coords_block(const ROMol &mol,
    // RDKit✔️✔️:                              const std::vector<unsigned int> &atomOrder) {
    // RDKit✔️✔️:   const auto &conf = mol.getConformer();
    // RDKit✔️✔️:   for (auto idx : atomOrder) {
    // RDKit✔️✔️:     const auto &pt = conf.getAtomPos(idx);
    // RDKit✔️✔️:     res += boost::str(boost::format("%g,%g,") % zero_small_vals(pt.x) %
    // RDKit✔️✔️:                       zero_small_vals(pt.y));
    // RDKit✔️✔️:     if (conf.is3D()) {
    // RDKit✔️✔️:       auto zc = boost::str(boost::format("%g") % zero_small_vals(pt.z));
    // RDKit✔️✔️:       if (zc != "0") {
    // RDKit✔️✔️:         res += zc;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    if let Some(conformer) = query.conformers_3d().first() {
        return Some(
            atom_order
                .iter()
                .map(|atom| {
                    let point = conformer.coordinates()[atom.index()];
                    let mut text = format!(
                        "{},{},",
                        query_cx_format_general(query_cx_zero_small(point[0])),
                        query_cx_format_general(query_cx_zero_small(point[1]))
                    );
                    if conformer.is_3d() {
                        let z = query_cx_format_general(query_cx_zero_small(point[2]));
                        if z != "0" {
                            text.push_str(&z);
                        }
                    }
                    text
                })
                .collect::<Vec<_>>()
                .join(";"),
        );
    }
    let coordinates = query.coordinates_2d()?;
    Some(
        atom_order
            .iter()
            .map(|atom| {
                let point = coordinates[atom.index()];
                format!(
                    "{},{},",
                    query_cx_format_general(query_cx_zero_small(point[0])),
                    query_cx_format_general(query_cx_zero_small(point[1]))
                )
            })
            .collect::<Vec<_>>()
            .join(";"),
    )
}

fn write_query_cx_atom_labels(query: &QueryGraph, atom_order: &[AtomId]) -> String {
    // RDKit✔️✔️: if (atom->getPropIfPresent(common_properties::_QueryAtomGenericLabel,
    // RDKit✔️✔️:                            lbl)) {
    // RDKit✔️✔️:   res += quote_string(lbl + "_p");
    // RDKit✔️✔️: } else if (!atom->getAtomicNum() &&
    // RDKit✔️✔️:            atom->getPropIfPresent(common_properties::dummyLabel, lbl) &&
    // RDKit✔️✔️:            std::find(SmilesParseOps::pseudoatoms.begin(),
    // RDKit✔️✔️:                      SmilesParseOps::pseudoatoms.end(), lbl) !=
    // RDKit✔️✔️:                SmilesParseOps::pseudoatoms.end()) {
    // RDKit✔️✔️:   res += quote_string(lbl + "_p");
    // RDKit✔️✔️: } else if (!atom->getAtomicNum() &&
    // RDKit✔️✔️:            atom->getPropIfPresent(common_properties::_fromAttachPoint,
    // RDKit✔️✔️:                                       val) &&
    // RDKit✔️✔️:            (val == 1 || val == 2)) {
    // RDKit✔️✔️:   res += quote_string("_AP" + std::to_string(val));
    // RDKit✔️✔️: } else if (atom->getPropIfPresent(common_properties::atomLabel, lbl)) {
    // RDKit✔️✔️:   res += quote_string(lbl);
    // RDKit✔️✔️: }
    const PSEUDOATOMS: [&str; 2] = ["Pol", "Mod"];
    let labels = atom_order
        .iter()
        .map(|atom_id| {
            let atom = &query.atoms()[atom_id.index()];
            if let Some(label) = atom.prop("_QueryAtomGenericLabel") {
                format!("{label}_p")
            } else if atom.atomic_number() == 0
                && atom
                    .prop("dummyLabel")
                    .is_some_and(|label| PSEUDOATOMS.contains(&label))
            {
                format!("{}_p", atom.prop("dummyLabel").unwrap_or_default())
            } else if atom.atomic_number() == 0
                && atom
                    .prop("_fromAttachPoint")
                    .is_some_and(|value| matches!(value, "1" | "2"))
            {
                format!("_AP{}", atom.prop("_fromAttachPoint").unwrap_or_default())
            } else {
                atom.prop("atomLabel").unwrap_or_default().to_owned()
            }
        })
        .collect::<Vec<_>>();
    if labels.iter().all(String::is_empty) {
        String::new()
    } else {
        labels.join(";")
    }
}

fn write_query_cx_atom_values(query: &QueryGraph, atom_order: &[AtomId]) -> String {
    // RDKit✔️✔️: if (mol.getAtomWithIdx(idx)->getPropIfPresent(prop, lbl)) {
    // RDKit✔️✔️:   res += quote_string(lbl);
    // RDKit✔️✔️: }
    let values = atom_order
        .iter()
        .map(|atom| {
            query.atoms()[atom.index()]
                .prop("molFileValue")
                .unwrap_or_default()
                .to_owned()
        })
        .collect::<Vec<_>>();
    if values.iter().all(String::is_empty) {
        String::new()
    } else {
        values.join(";")
    }
}

fn write_query_cx_radicals(query: &QueryGraph, atom_order: &[AtomId]) -> String {
    // RDKit✔️✔️: std::map<unsigned int, std::vector<unsigned int>> rads;
    // RDKit✔️✔️: for (unsigned int i = 0; i < atomOrder.size(); ++i) {
    // RDKit✔️✔️:   auto nrad = mol.getAtomWithIdx(atomOrder[i])->getNumRadicalElectrons();
    // RDKit✔️✔️:   if (nrad) {
    // RDKit✔️✔️:     rads[nrad].push_back(i);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: switch (pr.first) {
    // RDKit✔️✔️:   case 1: res += "^1:"; break;
    // RDKit✔️✔️:   case 2: res += "^2:"; break;
    // RDKit✔️✔️:   case 3: res += "^5:"; break;
    // RDKit✔️✔️:   default: BOOST_LOG(rdWarningLog) << "unsupported number of radical electrons ";
    // RDKit✔️✔️: }
    let mut radicals = BTreeMap::<u8, Vec<usize>>::new();
    for (position, atom) in atom_order.iter().copied().enumerate() {
        let count = query.atoms()[atom.index()].radical_electrons();
        if count != 0 {
            radicals.entry(count).or_default().push(position);
        }
    }
    let mut output = String::new();
    for (count, atoms) in radicals {
        match count {
            1 => output.push_str("^1:"),
            2 => output.push_str("^2:"),
            3 => output.push_str("^5:"),
            _ => {}
        }
        for atom in atoms {
            output.push_str(&atom.to_string());
            output.push(',');
        }
    }
    while output.ends_with(',') {
        output.pop();
    }
    output
}

fn quote_query_cx_atom_property(text: &str) -> String {
    // RDKit✔️✔️: for (auto c : txt) {
    // RDKit✔️✔️:   if (c == '.') {
    // RDKit✔️✔️:     res += "&#46;";
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     res += c;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    text.replace('.', "&#46;")
}

fn write_query_cx_atom_properties(query: &QueryGraph, atom_order: &[AtomId]) -> String {
    // RDKit✔️✔️: constexpr std::array<std::string_view, 7> skip = {
    // RDKit✔️✔️:     common_properties::atomLabel, common_properties::molFileValue,
    // RDKit✔️✔️:     common_properties::molParity, common_properties::molAtomMapNumber,
    // RDKit✔️✔️:     common_properties::molStereoCare, common_properties::molRxnExactChange,
    // RDKit✔️✔️:     common_properties::molInversionFlag};
    // RDKit✔️✔️: for (const auto &pn : atom->getPropList(includePrivate, includeComputed)) {
    // RDKit✔️✔️:   if (std::find(skip.begin(), skip.end(), pn) == skip.end()) {
    // RDKit✔️✔️:     res += boost::str(boost::format(":%d.%s.%s") % which %
    // RDKit✔️✔️:                       quote_atomprop_string(pn) % quote_atomprop_string(pv));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    const SKIP: [&str; 7] = [
        "atomLabel",
        "molFileValue",
        "molParity",
        "molAtomMapNumber",
        "molStereoCare",
        "molRxnExactChange",
        "molInversionFlag",
    ];
    const PSEUDOATOMS: [&str; 2] = ["Pol", "Mod"];
    let mut entries = Vec::new();
    for (position, atom_id) in atom_order.iter().copied().enumerate() {
        let atom = &query.atoms()[atom_id.index()];
        let attachment = atom.atomic_number() == 0 && atom.prop("_fromAttachPoint").is_some();
        for (name, value) in atom.props() {
            if name.starts_with('_') || atom.is_prop_computed(name) || SKIP.contains(&name.as_str())
            {
                continue;
            }
            if name == "dummyLabel"
                && (attachment || value == "*" || PSEUDOATOMS.contains(&value.as_str()))
            {
                continue;
            }
            entries.push(format!(
                "{position}.{}.{}",
                quote_query_cx_atom_property(name),
                quote_query_cx_atom_property(value)
            ));
        }
    }
    if entries.is_empty() {
        String::new()
    } else {
        format!("atomProp:{}", entries.join(":"))
    }
}

fn write_query_cx_bond_config(
    query: &QueryGraph,
    atom_order: &[AtomId],
    bond_order: &[BondId],
    coordinates_included: bool,
) -> String {
    // RDKit✔️✔️: if (!canHaveDirection(*bond)) {
    // RDKit✔️✔️:   continue;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (bd == Bond::BondDir::UNKNOWN) {
    // RDKit✔️✔️:   wType = "w";
    // RDKit✔️✔️: } else if (coordsIncluded || isAnAtropisomer) {
    // RDKit✔️✔️:   if (bd == Bond::BondDir::BEGINWEDGE) {
    // RDKit✔️✔️:     wType = "wU";
    // RDKit✔️✔️:   } else if (bd == Bond::BondDir::BEGINDASH) {
    // RDKit✔️✔️:     wType = "wD";
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    let positions = query_cx_atom_positions(atom_order, query.num_atoms());
    let mut parts = BTreeMap::<&'static str, Vec<String>>::new();
    for (bond_position, bond_id) in bond_order.iter().copied().enumerate() {
        let bond = query.bonds()[bond_id.index()].bond();
        if !matches!(bond.order(), BondOrder::Single | BondOrder::Aromatic) {
            continue;
        }
        let mut direction = match bond.direction() {
            BondDirection::BeginWedge | BondDirection::BeginDash | BondDirection::Unknown => {
                bond.direction()
            }
            _ => BondDirection::None,
        };
        if direction == BondDirection::None {
            direction = match bond
                .prop("_MolFileBondCfg")
                .and_then(|value| value.parse::<u8>().ok())
            {
                Some(1) => BondDirection::BeginWedge,
                Some(2) => BondDirection::Unknown,
                Some(3) => BondDirection::BeginDash,
                _ => BondDirection::None,
            };
        }
        let kind = match direction {
            BondDirection::Unknown => Some("w"),
            BondDirection::BeginWedge if coordinates_included => Some("wU"),
            BondDirection::BeginDash if coordinates_included => Some("wD"),
            _ => None,
        };
        let (Some(kind), Some(begin)) = (kind, positions[bond.begin().index()]) else {
            continue;
        };
        parts
            .entry(kind)
            .or_default()
            .push(format!("{begin}.{bond_position}"));
    }
    parts
        .into_iter()
        .map(|(kind, entries)| format!("{kind}:{}", entries.join(",")))
        .collect::<Vec<_>>()
        .join(",")
}

fn write_query_cx_typed_bonds(
    query: &QueryGraph,
    atom_order: &[AtomId],
    bond_order: &[BondId],
    order: BondOrder,
    symbol: &str,
) -> String {
    // RDKit✔️✔️: for (unsigned int i = 0; i < bondOrder.size(); ++i) {
    // RDKit✔️✔️:   auto idx = bondOrder[i];
    // RDKit✔️✔️:   const auto bond = mol.getBondWithIdx(idx);
    // RDKit✔️✔️:   if (bond->getBondType() != bondType) {
    // RDKit✔️✔️:     continue;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   auto begAtomOrder = std::find(atomOrder.begin(), atomOrder.end(),
    // RDKit✔️✔️:                                  bond->getBeginAtomIdx()) - atomOrder.begin();
    // RDKit✔️✔️:   res += boost::str(boost::format("%d.%d") % begAtomOrder % i);
    // RDKit✔️✔️: }
    let positions = query_cx_atom_positions(atom_order, query.num_atoms());
    let entries = bond_order
        .iter()
        .copied()
        .enumerate()
        .filter_map(|(position, bond_id)| {
            let bond = query.bonds()[bond_id.index()].bond();
            let matches = if order == BondOrder::Dative {
                matches!(bond.order(), BondOrder::Dative | BondOrder::DativeOne)
            } else {
                bond.order() == order
            };
            matches.then(|| {
                positions[bond.begin().index()].map(|begin| format!("{begin}.{position}"))
            })?
        })
        .collect::<Vec<_>>();
    if entries.is_empty() {
        String::new()
    } else {
        format!("{symbol}:{}", entries.join(","))
    }
}

fn write_query_cx_zero_bonds(query: &QueryGraph, bond_order: &[BondId]) -> String {
    // RDKit✔️✔️: if (bond->getBondType() != Bond::BondType::ZERO) {
    // RDKit✔️✔️:   continue;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: res += boost::str(boost::format("%d") % i);
    let entries = bond_order
        .iter()
        .copied()
        .enumerate()
        .filter_map(|(position, bond)| {
            (query.bonds()[bond.index()].bond().order() == BondOrder::Zero)
                .then(|| position.to_string())
        })
        .collect::<Vec<_>>();
    if entries.is_empty() {
        String::new()
    } else {
        format!("Z:{}", entries.join(","))
    }
}

fn query_cx_stereo_kind_order(kind: StereoGroupKind) -> u8 {
    match kind {
        StereoGroupKind::Absolute => 0,
        StereoGroupKind::Or => 1,
        StereoGroupKind::And => 2,
    }
}

fn assign_query_cx_stereo_group_ids(groups: &[StereoGroup]) -> Vec<Option<u32>> {
    // RDKit✔️✔️: if (sg.getWriteId() == 0) {
    // RDKit✔️✔️:   ++nextId;
    // RDKit✔️✔️:   while (nextId < ids.size() && ids[nextId]) {
    // RDKit✔️✔️:     ++nextId;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   sg.setWriteId(nextId);
    // RDKit✔️✔️: }
    // Query parsing retains source read IDs. The current detached model has no
    // separate write-ID carrier, matching parsed source groups whose write ID
    // is zero; assignment is therefore independent within OR and AND kinds.
    let mut ids = vec![None; groups.len()];
    let mut next_or = 0;
    let mut next_and = 0;
    for (index, group) in groups.iter().enumerate() {
        let next = match group.kind() {
            StereoGroupKind::Absolute => continue,
            StereoGroupKind::Or => &mut next_or,
            StereoGroupKind::And => &mut next_and,
        };
        *next += 1;
        ids[index] = Some(*next);
    }
    ids
}

fn write_query_cx_enhanced_stereo(query: &QueryGraph, atom_order: &[AtomId]) -> String {
    // RDKit✔️✔️: const auto newAtomIndexes = getSortedMappedIndexes(atomIds, revOrder);
    // RDKit✔️✔️: if (!newAtomIndexes.empty()) {
    // RDKit✔️✔️:   sortingGroups.emplace_back(sg, newAtomIndexes);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: // sort by 1) StereoGroup type; 2) StereoGroup id; 3) atom indexes
    // RDKit✔️✔️: assignStereoGroupIds(groups);
    // RDKit✔️✔️: case StereoGroupType::STEREO_ABSOLUTE: res << "a:"; break;
    // RDKit✔️✔️: case StereoGroupType::STEREO_OR: res << "o" << sgItr->getWriteId() << ":"; break;
    // RDKit✔️✔️: case StereoGroupType::STEREO_AND: res << "&" << sgItr->getWriteId() << ":"; break;
    let positions = query_cx_source_reverse_atom_order(atom_order, query.num_atoms());
    let mut groups = query
        .stereo_groups()
        .iter()
        .filter_map(|group| {
            let mut atoms = group
                .atoms()
                .iter()
                .map(|atom| positions[atom.index()])
                .collect::<Vec<_>>();
            atoms.sort_unstable();
            (!atoms.is_empty()).then(|| (group.clone(), atoms))
        })
        .collect::<Vec<_>>();
    groups.sort_by(|(left_group, left_atoms), (right_group, right_atoms)| {
        query_cx_stereo_kind_order(left_group.kind())
            .cmp(&query_cx_stereo_kind_order(right_group.kind()))
            .then_with(|| left_atoms.cmp(right_atoms))
    });
    let sorted_groups = groups
        .iter()
        .map(|(group, _)| group.clone())
        .collect::<Vec<_>>();
    let ids = assign_query_cx_stereo_group_ids(&sorted_groups);
    groups
        .into_iter()
        .zip(ids)
        .map(|((group, atoms), id)| {
            let prefix = match group.kind() {
                StereoGroupKind::Absolute => "a".to_owned(),
                StereoGroupKind::Or => format!("o{}", id.expect("OR group has assigned ID")),
                StereoGroupKind::And => format!("&{}", id.expect("AND group has assigned ID")),
            };
            format!(
                "{prefix}:{}",
                atoms
                    .iter()
                    .map(usize::to_string)
                    .collect::<Vec<_>>()
                    .join(",")
            )
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn query_cx_other_atom(bond: &Bond, atom: AtomId) -> Option<AtomId> {
    if bond.begin() == atom {
        Some(bond.end())
    } else if bond.end() == atom {
        Some(bond.begin())
    } else {
        None
    }
}

fn write_query_cx_ring_bond_stereo(
    query: &QueryGraph,
    atom_order: &[AtomId],
    bond_order: &[BondId],
) -> Result<String, SmartsWriteError> {
    // RDKit✔️❌: if (!mol.getRingInfo()->isInitialized()) {
    // RDKit✔️❌:   return "";
    // RDKit✔️❌: }
    // RDKit✔️❌: if (!rinfo->numBondRings(idx) ||
    // RDKit✔️❌:     rinfo->minBondRingSize(idx) <
    // RDKit✔️❌:         Chirality::minRingSizeForDoubleBondStereo) {
    // RDKit✔️❌:   continue;
    // RDKit✔️❌: }
    // RDKit✔️❌: if (bstereo == Bond::BondStereo::STEREOANY) {
    // RDKit✔️❌:   ctu += label;
    // RDKit✔️❌: } else if (bstereo == Bond::BondStereo::STEREOCIS || needSwap) {
    // RDKit✔️❌:   c += label;
    // RDKit✔️❌: } else {
    // RDKit✔️❌:   t += label;
    // RDKit✔️❌: }
    // The shared ring owner accepts ordinary Bond rows, so this detached
    // QueryGraph adapter performs one linear carrier clone before the same
    // ring algorithm. That extra allocation is a known source performance gap.
    let bonds = query
        .bonds()
        .iter()
        .map(|bond| bond.bond().clone())
        .collect::<Vec<_>>();
    let adjacency = AdjacencyList::try_from_topology(query.num_atoms(), &bonds)
        .map_err(|error| SmartsWriteError::InvalidGraph(error.to_string()))?;
    let rings = cosmolkit_core::fast_find_rings_from_parts(query.num_atoms(), &bonds, &adjacency)
        .map_err(|error| SmartsWriteError::InvalidGraph(error.to_string()))?;
    let positions = query_cx_atom_positions(atom_order, query.num_atoms());
    let mut cis = Vec::new();
    let mut trans = Vec::new();
    let mut unknown = Vec::new();
    for (bond_position, bond_id) in bond_order.iter().copied().enumerate() {
        if rings.num_bond_rings(bond_id) == 0 || rings.min_bond_ring_size(bond_id) < 8 {
            continue;
        }
        let bond = query.bonds()[bond_id.index()].bond();
        if !matches!(bond.order(), BondOrder::Double | BondOrder::Aromatic)
            || !matches!(
                bond.stereo(),
                BondStereo::Any | BondStereo::Cis | BondStereo::Trans
            )
        {
            continue;
        }
        if bond.stereo() == BondStereo::Any {
            unknown.push(bond_position.to_string());
            continue;
        }
        let Some([begin_reference, end_reference]) = bond.stereo_atoms() else {
            continue;
        };
        let mut swap = false;
        for (center, opposite, reference) in [
            (bond.begin(), bond.end(), begin_reference),
            (bond.end(), bond.begin(), end_reference),
        ] {
            let neighbors = &query.adjacency()[center.index()];
            if neighbors.len() <= 2 {
                continue;
            }
            let Some(reference_position) = positions[reference.index()] else {
                continue;
            };
            for &(_, incident) in neighbors {
                let neighbor = query_cx_other_atom(query.bonds()[incident].bond(), center)
                    .expect("query adjacency endpoint");
                if neighbor != opposite
                    && neighbor != reference
                    && positions[neighbor.index()]
                        .is_some_and(|position| position < reference_position)
                {
                    swap = !swap;
                }
            }
        }
        if bond.stereo() == BondStereo::Cis || swap {
            cis.push(bond_position.to_string());
        } else {
            trans.push(bond_position.to_string());
        }
    }
    Ok([
        (!cis.is_empty()).then(|| format!("c:{}", cis.join(","))),
        (!trans.is_empty()).then(|| format!("t:{}", trans.join(","))),
        (!unknown.is_empty()).then(|| format!("ctu:{}", unknown.join(","))),
    ]
    .into_iter()
    .flatten()
    .collect::<Vec<_>>()
    .join(","))
}

fn write_query_cx_link_nodes(query: &QueryGraph, atom_order: &[AtomId]) -> String {
    // RDKit✔️✔️: auto linkNodes = MolEnumerator::utils::getMolLinkNodes(mol, strict);
    // RDKit✔️✔️: if (linkNodes.empty()) {
    // RDKit✔️✔️:   return "";
    // RDKit✔️✔️: }
    // RDKit✔️✔️: for (const auto &ln : linkNodes) {
    // RDKit✔️✔️:   unsigned int atomIdx = atomOrder[ln.bondAtoms[0].first];
    // RDKit✔️✔️:   res << atomIdx << ":" << ln.minRep << "." << ln.maxRep;
    // RDKit✔️✔️: }
    let Some(raw) = query.prop("_MolFileLinkNodes") else {
        return String::new();
    };
    let mut entries = Vec::new();
    for item in raw.split('|').filter(|item| !item.trim().is_empty()) {
        let Ok(values) = item
            .split_whitespace()
            .map(str::parse::<usize>)
            .collect::<Result<Vec<_>, _>>()
        else {
            continue;
        };
        if values.len() < 7 || values[2] < 2 {
            continue;
        }
        let Some(center) = values[3].checked_sub(1) else {
            continue;
        };
        let Some(center_output) = atom_order.get(center).map(|atom| atom.index()) else {
            continue;
        };
        let mut entry = format!("{center_output}:{}.{}", values[0], values[1]);
        if query.adjacency().get(center).map_or(0, Vec::len) > 2 {
            let (Some(first), Some(second)) = (values[4].checked_sub(1), values[6].checked_sub(1))
            else {
                continue;
            };
            let (Some(first_output), Some(second_output)) = (
                atom_order.get(first).map(|atom| atom.index()),
                atom_order.get(second).map(|atom| atom.index()),
            ) else {
                continue;
            };
            entry.push_str(&format!(".{first_output}.{second_output}"));
        }
        entries.push(entry);
    }
    if entries.is_empty() {
        String::new()
    } else {
        format!("LN:{}", entries.join(","))
    }
}

fn query_cx_is_data_sgroup(group: &SubstanceGroup) -> bool {
    matches!(group.kind(), SubstanceGroupKind::Data)
        || group
            .props()
            .get("TYPE")
            .is_some_and(|value| value == "DAT")
}

fn query_cx_sgroup_value(group: &SubstanceGroup, key: &str) -> String {
    group.props().get(key).cloned().unwrap_or_default()
}

fn write_query_cx_data_sgroups(query: &QueryGraph, atom_order: &[AtomId]) -> String {
    // RDKit✔️✔️: if (sg.hasProp("TYPE") && sg.getProp<std::string>("TYPE") == "DAT") {
    // RDKit✔️✔️:   res << "SgD:";
    // RDKit✔️✔️:   for (const auto oaid : sg.getAtoms()) {
    // RDKit✔️✔️:     res << revOrder[oaid] << ",";
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   res << ":" << FIELDNAME << ":" << DATAFIELDS << ":" << QUERYOP
    // RDKit✔️✔️:       << ":" << FIELDINFO << ":" << FIELDTAG << ":";
    // RDKit✔️✔️: }
    let positions = query_cx_source_reverse_atom_order(atom_order, query.num_atoms());
    query_substance_groups(query)
        .iter()
        .filter(|group| query_cx_is_data_sgroup(group))
        .filter_map(|group| {
            let atoms = group
                .atoms()
                .iter()
                .map(|atom| positions[atom.index()])
                .map(|position| position.to_string())
                .collect::<Vec<_>>();
            if atoms.is_empty() {
                return None;
            }
            let data = if group.data_fields().is_empty() {
                query_cx_sgroup_value(group, "DATAFIELDS")
            } else {
                group.data_fields().join(",")
            };
            Some(format!(
                "SgD:{}:{}:{}:{}:{}:{}:",
                atoms.join(","),
                query_cx_sgroup_value(group, "FIELDNAME"),
                data,
                query_cx_sgroup_value(group, "QUERYOP"),
                query_cx_sgroup_value(group, "FIELDINFO"),
                query_cx_sgroup_value(group, "FIELDTAG")
            ))
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn query_cx_polymer_type(group: &SubstanceGroup) -> Option<&'static str> {
    match group.kind() {
        SubstanceGroupKind::StructuralRepeatUnit => Some("n"),
        SubstanceGroupKind::Monomer => Some("mon"),
        SubstanceGroupKind::Mer => Some("mer"),
        SubstanceGroupKind::Copolymer => match group
            .subtype()
            .or_else(|| group.props().get("SUBTYPE").map(String::as_str))
        {
            Some("ALT") => Some("alt"),
            Some("RAN") => Some("ran"),
            Some("BLO") => Some("blk"),
            _ => Some("co"),
        },
        SubstanceGroupKind::Crosslink => Some("xl"),
        SubstanceGroupKind::Modification => Some("mod"),
        SubstanceGroupKind::MixtureComponent => Some("mix"),
        SubstanceGroupKind::Formulation => Some("f"),
        SubstanceGroupKind::AnyPolymer => Some("any"),
        SubstanceGroupKind::Graft => Some("grf"),
        SubstanceGroupKind::Generic(value) if value == "GEN" => Some("gen"),
        SubstanceGroupKind::Generic(value) if value == "COM" => Some("c"),
        _ => None,
    }
}

fn query_cx_connection_text(group: &SubstanceGroup) -> String {
    group
        .connection()
        .map(|connection| match connection {
            SGroupConnection::HeadToHead => "hh".to_owned(),
            SGroupConnection::HeadToTail => "ht".to_owned(),
            SGroupConnection::Either => "eu".to_owned(),
            SGroupConnection::Unknown(value) => value.to_ascii_lowercase(),
        })
        .or_else(|| {
            group
                .props()
                .get("CONNECT")
                .map(|value| value.to_ascii_lowercase())
        })
        .unwrap_or_default()
}

fn write_query_cx_polymer_sgroups(
    query: &QueryGraph,
    atom_order: &[AtomId],
    bond_order: &[BondId],
) -> Result<String, SmartsWriteError> {
    // RDKit✔️✔️: if (sg.getPropIfPresent("TYPE", typ) &&
    // RDKit✔️✔️:     reverseTypemap.find(typ) != reverseTypemap.end()) {
    // RDKit✔️✔️:   res << "Sg:" << reverse type << ":";
    // RDKit✔️✔️:   for (const auto oaid : sg.getAtoms()) {
    // RDKit✔️✔️:     res << revAtomOrder[oaid] << ",";
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   res << ":" << LABEL << ":" << lowercase CONNECT << ":";
    // RDKit✔️✔️: }
    let atom_positions = query_cx_source_reverse_atom_order(atom_order, query.num_atoms());
    let mut blocks = Vec::new();
    for group in query_substance_groups(query) {
        let Some(kind) = query_cx_polymer_type(group) else {
            continue;
        };
        let atoms = group
            .atoms()
            .iter()
            .map(|atom| atom_positions[atom.index()])
            .map(|position| position.to_string())
            .collect::<Vec<_>>();
        if atoms.is_empty() {
            continue;
        }
        let crossing_position = |bond: BondId| {
            bond_order
                .get(bond.index())
                .map(|output_bond| output_bond.index())
                .ok_or_else(|| {
                    SmartsWriteError::InvalidGraph(format!(
                        "SGroup {} crossing bond {} is absent from CX bond order",
                        group.id().index(),
                        bond.index()
                    ))
                })
        };
        let head = if group.head_crossing_bonds().len() > 1 {
            group
                .head_crossing_bonds()
                .iter()
                .copied()
                .map(crossing_position)
                .collect::<Result<Vec<_>, _>>()?
                .into_iter()
                .map(|position| position.to_string())
                .collect::<Vec<_>>()
                .join(",")
        } else {
            String::new()
        };
        let tail = if group.crossing_bond_correspondence().len() > 2 {
            group
                .crossing_bond_correspondence()
                .iter()
                .skip(1)
                .step_by(2)
                .copied()
                .map(crossing_position)
                .collect::<Result<Vec<_>, _>>()?
                .into_iter()
                .map(|position| position.to_string())
                .collect::<Vec<_>>()
                .join(",")
        } else {
            String::new()
        };
        blocks.push(format!(
            "Sg:{kind}:{}:{}:{}:{head}:{tail}:",
            atoms.join(","),
            group
                .label()
                .or_else(|| group.props().get("LABEL").map(String::as_str))
                .unwrap_or_default(),
            query_cx_connection_text(group)
        ));
    }
    Ok(blocks.join(","))
}

fn query_cx_sgroup_index(group: &SubstanceGroup) -> usize {
    group
        .props()
        .get("index")
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or_else(|| group.id().index())
}

fn write_query_cx_sgroup_hierarchy(query: &QueryGraph) -> String {
    // RDKit✔️✔️: if (sg.hasProp("_cxsmilesOutputIndex")) {
    // RDKit✔️✔️:   unsigned int sgidx = sg.getIndexInMol();
    // RDKit✔️✔️:   sg.getPropIfPresent("index", sgidx);
    // RDKit✔️✔️:   sgroupOrder[sgidx] = sg.getProp<unsigned int>("_cxsmilesOutputIndex");
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (sg.getPropIfPresent("PARENT", pidx) &&
    // RDKit✔️✔️:     sgroupOrder.find(pidx) != sgroupOrder.end()) {
    // RDKit✔️✔️:   accum[sgroupOrder[pidx]].push_back(sgroupOrder[sgidx]);
    // RDKit✔️✔️: }
    let groups = query_substance_groups(query);
    let mut output_indices = BTreeMap::new();
    let mut next = 0usize;
    for group in groups.iter().filter(|group| query_cx_is_data_sgroup(group)) {
        output_indices.insert(query_cx_sgroup_index(group), next);
        next += 1;
    }
    for group in groups
        .iter()
        .filter(|group| query_cx_polymer_type(group).is_some())
    {
        output_indices.insert(query_cx_sgroup_index(group), next);
        next += 1;
    }
    let source_indices = groups
        .iter()
        .map(|group| (group.id(), query_cx_sgroup_index(group)))
        .collect::<BTreeMap<_, _>>();
    let mut hierarchy = BTreeMap::<usize, Vec<usize>>::new();
    for group in groups {
        let Some(child) = output_indices.get(&query_cx_sgroup_index(group)).copied() else {
            continue;
        };
        let parent_key = group
            .parent()
            .and_then(|parent| source_indices.get(&parent).copied())
            .or_else(|| {
                group
                    .props()
                    .get("PARENT")
                    .and_then(|value| value.parse::<usize>().ok())
            });
        let Some(parent) = parent_key.and_then(|parent| output_indices.get(&parent).copied())
        else {
            continue;
        };
        hierarchy.entry(parent).or_default().push(child);
    }
    if hierarchy.is_empty() {
        String::new()
    } else {
        format!(
            "SgH:{}",
            hierarchy
                .into_iter()
                .map(|(parent, children)| format!(
                    "{parent}:{}",
                    children
                        .iter()
                        .map(usize::to_string)
                        .collect::<Vec<_>>()
                        .join(".")
                ))
                .collect::<Vec<_>>()
                .join(",")
        )
    }
}

fn write_query_cx_extensions(
    query: &QueryGraph,
    atom_order: &[AtomId],
    bond_order: &[BondId],
) -> Result<String, SmartsWriteError> {
    // RDKit✔️✔️: std::string getCXExtensions(const ROMol &mol, std::uint32_t flags) {
    // RDKit✔️✔️:   std::string res = "|";
    // RDKit✔️✔️:   if ((flags & SmilesWrite::CXSmilesFields::CX_COORDS) &&
    // RDKit✔️✔️:       mol.getNumConformers()) {
    // RDKit✔️✔️:     res += "(" + get_coords_block(mol, atomOrder) + ")";
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // labels, values, radicals, atom properties, bond configuration,
    // RDKit✔️✔️:   // coordinate/hydrogen/zero bonds, link nodes, enhanced stereo,
    // RDKit✔️✔️:   // data/polymer SGroups and hierarchy follow in this exact order.
    // RDKit✔️✔️:   if (res.size() > 1) {
    // RDKit✔️✔️:     res += "|";
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     res = "";
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    let mut result = String::from("|");
    let coordinates = write_query_cx_coordinates(query, atom_order);
    if let Some(coordinates) = &coordinates {
        result.push('(');
        result.push_str(coordinates);
        result.push(')');
    }
    let labels = write_query_cx_atom_labels(query, atom_order);
    if !labels.is_empty() {
        append_query_cx_extension(format!("${labels}$"), &mut result);
    }
    let values = write_query_cx_atom_values(query, atom_order);
    if !values.is_empty() {
        append_query_cx_extension(format!("$_AV:{values}$"), &mut result);
    }
    append_query_cx_extension(write_query_cx_radicals(query, atom_order), &mut result);
    append_query_cx_extension(
        write_query_cx_atom_properties(query, atom_order),
        &mut result,
    );
    append_query_cx_extension(
        write_query_cx_bond_config(query, atom_order, bond_order, coordinates.is_some()),
        &mut result,
    );
    append_query_cx_extension(
        write_query_cx_ring_bond_stereo(query, atom_order, bond_order)?,
        &mut result,
    );
    append_query_cx_extension(
        write_query_cx_typed_bonds(query, atom_order, bond_order, BondOrder::Dative, "C"),
        &mut result,
    );
    append_query_cx_extension(
        write_query_cx_typed_bonds(query, atom_order, bond_order, BondOrder::Hydrogen, "H"),
        &mut result,
    );
    append_query_cx_extension(write_query_cx_zero_bonds(query, bond_order), &mut result);
    append_query_cx_extension(write_query_cx_link_nodes(query, atom_order), &mut result);
    append_query_cx_extension(
        write_query_cx_enhanced_stereo(query, atom_order),
        &mut result,
    );
    append_query_cx_extension(write_query_cx_data_sgroups(query, atom_order), &mut result);
    append_query_cx_extension(
        write_query_cx_polymer_sgroups(query, atom_order, bond_order)?,
        &mut result,
    );
    append_query_cx_extension(write_query_cx_sgroup_hierarchy(query), &mut result);
    if result.len() == 1 {
        Ok(String::new())
    } else {
        result.push('|');
        Ok(result)
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
    result
        .smarts
        .push_str(&query_atom_to_smarts(query_atom, params)?);
    result.atom_ordering.push(atom);
    for (bond, first, _second, ring_number) in ring_edges
        .iter()
        .filter(|(_, first, second, _)| *first == atom || *second == atom)
    {
        if *first == atom {
            result.smarts.push_str(&query_bond_to_smarts(
                query
                    .bond(bond.index())
                    .ok_or(SmartsWriteError::FragmentBondOutOfRange { bond: bond.index() })?,
                params,
                Some(atom.index()),
            )?);
            result.bond_ordering.push(*bond);
        }
        if *ring_number < 10 {
            result.smarts.push_str(&ring_number.to_string());
        } else {
            result.smarts.push('%');
            result.smarts.push_str(&ring_number.to_string());
        }
    }
    let children = &tree_children[atom.index()];
    for (index, (bond, other)) in children.iter().enumerate() {
        if index + 1 != children.len() {
            result.smarts.push('(');
        }
        result.smarts.push_str(&query_bond_to_smarts(
            query
                .bond(bond.index())
                .ok_or(SmartsWriteError::FragmentBondOutOfRange { bond: bond.index() })?,
            params,
            Some(atom.index()),
        )?);
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
            result.smarts.push(')');
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
                atom,
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
                atom,
                predicate,
                &mut need_paren,
                true,
                params.do_isomeric_smiles,
                &mut stereo_written,
            )?;
            needs_brackets = need_paren;
            result
        }
        QueryNode::Not(child) => {
            needs_brackets = true;
            let mut need_paren = false;
            let mut propagated_negation_written = false;
            let mut value = match child.as_ref() {
                QueryNode::Predicate(AtomQueryPredicate::RecursiveSmarts(query)) => {
                    propagated_negation_written = true;
                    get_recursive_structure_query_smarts(
                        query,
                        true,
                        params,
                        query_graph_to_smarts,
                    )?
                }
                QueryNode::Predicate(predicate) => get_atom_smarts_simple(
                    atom,
                    predicate,
                    &mut need_paren,
                    true,
                    params.do_isomeric_smiles,
                    &mut stereo_written,
                )?,
                _ => {
                    propagated_negation_written = true;
                    recurse_get_smarts(
                        atom,
                        child.as_ref(),
                        true,
                        &mut features,
                        params,
                        &mut stereo_written,
                        &mut query_graph_to_smarts,
                    )?
                }
            };
            if !propagated_negation_written {
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
        result = if result.is_empty() {
            symbol.to_owned()
        } else {
            format!("{symbol};{result}")
        };
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
    // RDKit✔️✔️: if ((descrip == "BondAnd") || (descrip == "BondOr")) {
    // RDKit✔️✔️:   // composite query
    // RDKit✔️✔️:   res = _recurseBondSmarts(bond, query, query->getNegation(), atomToLeftIdx,
    // RDKit✔️✔️:                            queryFeatures, params);
    // RDKit✔️✔️: } else {
    // RDKit✔️✔️:   // simple query
    // RDKit✔️✔️:   if (query->getNegation()) {
    // RDKit✔️✔️:     res = "!";
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   res += getBondSmartsSimple(bond, query, atomToLeftIdx, params);
    // RDKit✔️✔️: }
    // Complexity review: this dispatch is constant-time. Composite rendering
    // retains the source's O(n) recursive traversal and output allocation.
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
                        true,
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
    atom: &QueryAtom,
    query: &AtomQueryPredicate,
    need_paren: &mut bool,
    check_for_symbol: bool,
    do_isomeric_smarts: bool,
    stereo_written: &mut bool,
) -> Result<String, SmartsWriteError> {
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
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     BOOST_LOG(rdWarningLog) << "Cannot write SMARTS for query type : " << descrip
    // RDKit❗✔️:                            << ". Ignoring it." << std::endl;
    // RDKit❗✔️:     res << "*";
    // RDKit❗✔️:   }
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
    // COSMolKit intentionally returns a typed error for modeled predicates
    // without a completed writer branch; the project contract forbids the
    // source warning-and-wildcard information loss.
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
                let Some(element) = Element::from_atomic_number(*atomic_number) else {
                    return Err(SmartsWriteError::UnsupportedAtomQuery {
                        predicate: query.clone(),
                    });
                };
                let mut symbol = element.symbol().to_owned();
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
                AtomRangeBounds::LessEqual(value) => format!("{value}-"),
                AtomRangeBounds::GreaterEqual(value) => format!("-{value}"),
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
        | AtomQueryPredicate::UnsupportedFeature(_) => {
            return Err(SmartsWriteError::UnsupportedAtomQuery {
                predicate: query.clone(),
            });
        }
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
    Ok(result)
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
    atom: &QueryAtom,
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
                )?;
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
