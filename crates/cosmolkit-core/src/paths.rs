//! Source-backed graph paths and detached topology subsets.
//!
//! This module owns graph algorithms over validated detached model values. It
//! never accepts a live molecule or runtime capability.

use std::collections::{BTreeMap, VecDeque};

use cosmolkit_model::{
    AtomId, AtomMapping, AtomQueryPredicate, BondId, BondMapping, BondQueryPredicate, BondStereo,
    MappingValidationError, QueryAtom, QueryBond, QueryGraph, QueryGraphError, QueryNode,
    StereoGroup, SubstanceGroup, SubstanceGroupId, TopologyBlock, TopologyMapping,
    TopologyValidationError,
};

use crate::{PeriodicTableError, atomic_mass};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PathRepresentation {
    Bonds,
    Atoms,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GraphPath {
    Bonds(Vec<BondId>),
    Atoms(Vec<AtomId>),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PathSearchParams {
    pub representation: PathRepresentation,
    pub use_hydrogens: bool,
    pub rooted_at_atom: Option<AtomId>,
    pub only_shortest_paths: bool,
}

impl Default for PathSearchParams {
    fn default() -> Self {
        Self {
            representation: PathRepresentation::Bonds,
            use_hydrogens: false,
            rooted_at_atom: None,
            only_shortest_paths: false,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct SubgraphSearchParams {
    pub use_hydrogens: bool,
    pub rooted_at_atom: Option<AtomId>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UniqueSubgraphParams {
    pub use_hydrogens: bool,
    pub use_bond_orders: bool,
    pub rooted_at_atom: Option<AtomId>,
    pub extra_atom_invariants: Option<Vec<u32>>,
}

impl Default for UniqueSubgraphParams {
    fn default() -> Self {
        Self {
            use_hydrogens: false,
            use_bond_orders: true,
            rooted_at_atom: None,
            extra_atom_invariants: None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AtomEnvironmentParams {
    pub use_hydrogens: bool,
    pub enforce_radius: bool,
}

impl Default for AtomEnvironmentParams {
    fn default() -> Self {
        Self {
            use_hydrogens: false,
            enforce_radius: true,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConnectedComponents {
    pub atom_to_component: Vec<usize>,
    pub components: Vec<Vec<AtomId>>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AtomEnvironment {
    pub bonds: Vec<BondId>,
    pub atom_distances: Vec<Option<usize>>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct SubtopologyParams {
    pub copy_as_query: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub enum DetachedPathSubgraph {
    Concrete(TopologyBlock),
    Query {
        graph: QueryGraph,
        substance_groups: Vec<SubstanceGroup>,
    },
}

#[derive(Debug, Clone, PartialEq)]
pub struct SubtopologyResult {
    pub subgraph: DetachedPathSubgraph,
    pub mapping: TopologyMapping,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum PathError {
    #[error(transparent)]
    InvalidTopology(TopologyValidationError),
    #[error("{role} atom {atom} is out of range for {atom_count} atoms")]
    AtomOutOfRange {
        role: &'static str,
        atom: AtomId,
        atom_count: usize,
    },
    #[error("shortest-path endpoints must be distinct (both were atom {atom})")]
    EqualShortestPathEndpoints { atom: AtomId },
    #[error("path length range {lower}..={upper} is invalid")]
    InvalidLengthRange { lower: usize, upper: usize },
    #[error("path length {length} cannot be represented by the source unsigned range")]
    LengthOverflow { length: usize },
    #[error("atom path row {position} refers to atom {atom}, outside {atom_count} atoms")]
    AtomPathOutOfRange {
        position: usize,
        atom: AtomId,
        atom_count: usize,
    },
    #[error("extra atom invariants have {actual} rows, expected {expected}")]
    ExtraInvariantLength { actual: usize, expected: usize },
    #[error(transparent)]
    PeriodicTable(PeriodicTableError),
    #[error("detached subset topology is invalid: {0}")]
    InvalidSubsetTopology(TopologyValidationError),
    #[error(transparent)]
    QueryGraph(QueryGraphError),
    #[error(transparent)]
    Mapping(MappingValidationError),
    #[error(transparent)]
    Matrix(crate::MatrixError),
}

pub(crate) trait NeighborSource {
    fn atom_count(&self) -> usize;
    fn visit_neighbors(&self, atom: usize, visitor: &mut dyn FnMut(usize));
}

impl NeighborSource for TopologyBlock {
    fn atom_count(&self) -> usize {
        self.atoms.len()
    }

    fn visit_neighbors(&self, atom: usize, visitor: &mut dyn FnMut(usize)) {
        for neighbor in self.adjacency.neighbors_of(atom) {
            visitor(neighbor.atom_index);
        }
    }
}

pub(crate) fn connected_components_from_source(
    source: &impl NeighborSource,
) -> ConnectedComponents {
    // BEGIN RDKIT CPP FUNCTION MolOps::getMolFrags
    // RDKit✔️✔️: unsigned int getMolFrags(const ROMol &mol, INT_VECT &mapping) {
    // RDKit✔️✔️:   unsigned int natms = mol.getNumAtoms();
    // RDKit✔️✔️:   mapping.resize(natms);
    // RDKit✔️✔️:   return natms ? boost::connected_components(mol.getTopology(), &mapping[0])
    // RDKit✔️✔️:                : 0;
    // RDKit✔️✔️: };
    let mut atom_to_component = vec![usize::MAX; source.atom_count()];
    let mut components = Vec::new();
    for start in 0..source.atom_count() {
        if atom_to_component[start] != usize::MAX {
            continue;
        }
        let component = components.len();
        let mut queue = VecDeque::from([start]);
        atom_to_component[start] = component;
        while let Some(atom) = queue.pop_front() {
            source.visit_neighbors(atom, &mut |neighbor| {
                if atom_to_component[neighbor] == usize::MAX {
                    atom_to_component[neighbor] = component;
                    queue.push_back(neighbor);
                }
            });
        }
        components.push(Vec::new());
    }
    // RDKit✔️✔️:   INT_INT_VECT_MAP comMap;
    // RDKit✔️✔️:   for (unsigned int i = 0; i < mol.getNumAtoms(); i++) {
    // RDKit✔️✔️:     int mi = mapping[i];
    // RDKit✔️✔️:     if (comMap.find(mi) == comMap.end()) {
    // RDKit✔️✔️:       INT_VECT comp;
    // RDKit✔️✔️:       comMap[mi] = comp;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     comMap[mi].push_back(i);
    // RDKit✔️✔️:   }
    for (atom, component) in atom_to_component.iter().copied().enumerate() {
        components[component].push(AtomId::new(atom));
    }
    // RDKit✔️✔️:   for (INT_INT_VECT_MAP_CI mci = comMap.begin(); mci != comMap.end(); mci++) {
    // RDKit✔️✔️:     frags.push_back((*mci).second);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return rdcast<unsigned int>(frags.size());
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION MolOps::getMolFrags
    ConnectedComponents {
        atom_to_component,
        components,
    }
}

pub fn connected_components(topology: &TopologyBlock) -> Result<ConnectedComponents, PathError> {
    topology.validate().map_err(PathError::InvalidTopology)?;
    Ok(connected_components_from_source(topology))
}

pub fn shortest_path(
    topology: &TopologyBlock,
    begin: AtomId,
    end: AtomId,
) -> Result<Vec<AtomId>, PathError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::getShortestPath
    // RDKit✔️✔️: INT_LIST getShortestPath(const ROMol &mol, int aid1, int aid2) {
    // RDKit✔️✔️:   int nats = mol.getNumAtoms();
    // RDKit✔️✔️:   RANGE_CHECK(0, aid1, nats - 1);
    // RDKit✔️✔️:   RANGE_CHECK(0, aid2, nats - 1);
    // RDKit✔️✔️:   CHECK_INVARIANT(aid1 != aid2, "");
    topology.validate().map_err(PathError::InvalidTopology)?;
    validate_required_atom(topology, begin, "begin")?;
    validate_required_atom(topology, end, "end")?;
    if begin == end {
        return Err(PathError::EqualShortestPathEndpoints { atom: begin });
    }
    // RDKit✔️✔️:   INT_VECT pred(nats, -1);  // set all atoms to unprocessed state
    // RDKit✔️✔️:   pred[aid1] = -2;          // marks begin
    // RDKit✔️✔️:   pred[aid2] = -3;          // marks end
    // RDKit✔️✔️:   std::deque<int> bfsQ;
    // RDKit✔️✔️:   bfsQ.push_back(aid1);
    let mut predecessor = vec![None; topology.atoms.len()];
    let mut visited = vec![false; topology.atoms.len()];
    let mut queue = VecDeque::from([begin.index()]);
    visited[begin.index()] = true;
    let mut found = false;
    // RDKit✔️✔️:   while ((!done) && (bfsQ.size() > 0)) {
    // RDKit✔️✔️:     int curAid = bfsQ.front();
    // RDKit✔️✔️:     boost::tie(nbrIdx, endNbrs) =
    // RDKit✔️✔️:         mol.getAtomNeighbors(mol.getAtomWithIdx(curAid));
    while !found {
        let Some(current) = queue.pop_front() else {
            break;
        };
        for neighbor in topology.adjacency.neighbors_of(current) {
            let next = neighbor.atom_index;
            if next == end.index() {
                predecessor[next] = Some(current);
                found = true;
                break;
            }
            if !visited[next] {
                visited[next] = true;
                predecessor[next] = Some(current);
                queue.push_back(next);
            }
        }
    }
    // RDKit✔️✔️:   INT_LIST res;
    // RDKit✔️✔️:   if (done) {
    // RDKit✔️✔️:     int prev = aid2;
    // RDKit✔️✔️:     res.push_back(aid2);
    // RDKit✔️✔️:     while (!done) {
    // RDKit✔️✔️:       prev = pred[prev];
    // RDKit✔️✔️:       if (prev != aid1) { res.push_front(prev); }
    // RDKit✔️✔️:       else { done = true; }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     res.push_front(aid1);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION MolOps::getShortestPath
    if !found {
        return Ok(Vec::new());
    }
    let mut reversed = vec![end];
    let mut current = end.index();
    while current != begin.index() {
        current = predecessor[current].expect("a found BFS path has predecessors");
        reversed.push(AtomId::new(current));
    }
    reversed.reverse();
    Ok(reversed)
}

pub fn all_paths_of_length(
    topology: &TopologyBlock,
    target_length: usize,
    params: &PathSearchParams,
) -> Result<Vec<GraphPath>, PathError> {
    Ok(
        all_paths_in_range(topology, target_length, target_length, params)?
            .remove(&target_length)
            .unwrap_or_default(),
    )
}

pub fn all_paths_in_range(
    topology: &TopologyBlock,
    lower_length: usize,
    upper_length: usize,
    params: &PathSearchParams,
) -> Result<BTreeMap<usize, Vec<GraphPath>>, PathError> {
    // BEGIN RDKIT CPP FUNCTION findAllPathsOfLengthsMtoN
    // RDKit✔️✔️: PRECONDITION(lowerLen <= upperLen, "");
    // RDKit✔️✔️: double *distMat = onlyShortestPaths ? MolOps::getDistanceMat(mol) : nullptr;
    topology.validate().map_err(PathError::InvalidTopology)?;
    validate_range(lower_length, upper_length)?;

    let distances = params
        .only_shortest_paths
        .then(|| crate::matrices::unweighted_distance_steps(topology))
        .transpose()
        .map_err(PathError::Matrix)?;
    let adjacency =
        atom_adjacency_matrix(topology, params.only_shortest_paths || params.use_hydrogens);
    // RDKit✔️✔️:   if (useBonds) {
    // RDKit✔️✔️:     ++lowerLen;
    // RDKit✔️✔️:     ++upperLen;
    // RDKit✔️✔️:   }
    let (atom_lower, atom_upper) = match params.representation {
        PathRepresentation::Bonds => (
            lower_length
                .checked_add(1)
                .ok_or(PathError::LengthOverflow {
                    length: lower_length,
                })?,
            upper_length
                .checked_add(1)
                .ok_or(PathError::LengthOverflow {
                    length: upper_length,
                })?,
        ),
        PathRepresentation::Atoms => (lower_length, upper_length),
    };
    // RDKit✔️✔️:   INT_PATH_LIST_MAP atomPaths = Subgraphs::pathFinderHelper(
    // RDKit✔️✔️:       adjMat, dim, lowerLen, upperLen, rootedAtAtom, distMat);
    let atom_paths = path_finder_helper(
        &adjacency,
        topology.atoms.len(),
        atom_lower,
        atom_upper,
        params.rooted_at_atom,
        distances.as_deref(),
    );

    let mut result = BTreeMap::new();
    // RDKit✔️✔️:   if (!useBonds && lowerLen >= 1) {
    // RDKit✔️✔️:     res[1] = atomPaths[1];
    // RDKit✔️✔️:   }
    if params.representation == PathRepresentation::Atoms && atom_lower >= 1 {
        result.insert(
            1,
            atom_paths
                .get(&1)
                .into_iter()
                .flatten()
                .map(|row| GraphPath::Atoms(row.iter().copied().map(AtomId::new).collect()))
                .collect(),
        );
    }
    // RDKit✔️✔️:   for (unsigned int i = lowerLen; i <= upperLen; ++i) {
    // RDKit✔️✔️:     if (i <= 1) { continue; }
    // RDKit✔️✔️:     std::vector<boost::dynamic_bitset<>> invars;
    if params.representation == PathRepresentation::Bonds || atom_upper > 1 {
        for length in atom_lower..=atom_upper {
            if length <= 1 {
                continue;
            }
            let mut seen_bond_sets: Vec<Vec<bool>> = Vec::new();
            for atom_path in atom_paths.get(&length).into_iter().flatten() {
                let mut bond_set = vec![false; topology.bonds.len()];
                let mut bond_path = Vec::with_capacity(length - 1);
                for pair in atom_path.windows(2) {
                    let bond = bond_between(topology, pair[0], pair[1])
                        .expect("validated adjacency path must have a bond");
                    bond_set[bond.index()] = true;
                    bond_path.push(bond);
                }
                // RDKit✔️✔️:       if (std::find(invars.begin(), invars.end(), invar) == invars.end()) {
                // RDKit✔️✔️:         invars.push_back(invar);
                if seen_bond_sets.contains(&bond_set) {
                    continue;
                }
                seen_bond_sets.push(bond_set);
                match params.representation {
                    PathRepresentation::Bonds => result
                        .entry(length - 1)
                        .or_insert_with(Vec::new)
                        .push(GraphPath::Bonds(bond_path)),
                    PathRepresentation::Atoms => result
                        .entry(length)
                        .or_insert_with(Vec::new)
                        .push(GraphPath::Atoms(
                            atom_path.iter().copied().map(AtomId::new).collect(),
                        )),
                }
            }
        }
    }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION findAllPathsOfLengthsMtoN
    Ok(result)
}

pub fn all_subgraphs_of_length(
    topology: &TopologyBlock,
    target_length: usize,
    params: &SubgraphSearchParams,
) -> Result<Vec<Vec<BondId>>, PathError> {
    topology.validate().map_err(PathError::InvalidTopology)?;
    if target_length == 0 {
        return Ok(Vec::new());
    }
    let neighbors = bond_neighbor_map(topology, params.use_hydrogens);
    Ok(all_subgraphs_of_length_from_neighbors(
        topology,
        &neighbors,
        target_length,
        params.rooted_at_atom,
    ))
}

pub fn all_subgraphs_in_range(
    topology: &TopologyBlock,
    lower_length: usize,
    upper_length: usize,
    params: &SubgraphSearchParams,
) -> Result<BTreeMap<usize, Vec<Vec<BondId>>>, PathError> {
    // BEGIN RDKIT CPP FUNCTION findAllSubgraphsOfLengthsMtoN
    // RDKit✔️✔️: PRECONDITION(lowerLen <= upperLen, "");
    // RDKit✔️✔️: boost::dynamic_bitset<> forbidden(mol.getNumBonds());
    // RDKit✔️✔️: INT_INT_VECT_MAP nbrs;
    // RDKit✔️✔️: Subgraphs::getNbrsList(mol, useHs, nbrs);
    topology.validate().map_err(PathError::InvalidTopology)?;
    validate_range(lower_length, upper_length)?;
    let neighbors = bond_neighbor_map(topology, params.use_hydrogens);
    let mut result = (lower_length..=upper_length)
        .map(|length| (length, Vec::new()))
        .collect::<BTreeMap<_, _>>();
    if upper_length == 0 {
        return Ok(result);
    }
    let mut forbidden = vec![false; topology.bonds.len()];
    // RDKit✔️✔️:   for (auto nbi = nbrs.begin(); nbi != nbrs.end(); nbi++) {
    for (&start, adjacent) in &neighbors {
        if !root_allows_bond(topology, params.rooted_at_atom, start) || forbidden[start] {
            continue;
        }
        forbidden[start] = true;
        recurse_walk_range(
            &neighbors,
            vec![start],
            adjacent.clone(),
            lower_length,
            upper_length,
            forbidden.clone(),
            &mut result,
        );
    }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION findAllSubgraphsOfLengthsMtoN
    Ok(result)
}

pub fn unique_subgraphs_of_length(
    topology: &TopologyBlock,
    target_length: usize,
    params: &UniqueSubgraphParams,
) -> Result<Vec<Vec<BondId>>, PathError> {
    topology.validate().map_err(PathError::InvalidTopology)?;
    if let Some(extra) = &params.extra_atom_invariants
        && extra.len() != topology.atoms.len()
    {
        return Err(PathError::ExtraInvariantLength {
            actual: extra.len(),
            expected: topology.atoms.len(),
        });
    }
    // RDKit✔️✔️: PATH_LIST allSubgraphs =
    // RDKit✔️✔️:     findAllSubgraphsOfLengthN(mol, targetLen, useHs, rootedAtAtom);
    // RDKit✔️✔️: PATH_LIST res = Subgraphs::uniquifyPaths(mol, allSubgraphs, useBO);
    let all = all_subgraphs_of_length(
        topology,
        target_length,
        &SubgraphSearchParams {
            use_hydrogens: params.use_hydrogens,
            rooted_at_atom: params.rooted_at_atom,
        },
    )?;
    let mut result = Vec::new();
    let mut seen = Vec::new();
    for path in all {
        let discriminator = path_discriminator(
            topology,
            &path,
            params.use_bond_orders,
            params.extra_atom_invariants.as_deref(),
        )?;
        if !seen.contains(&discriminator) {
            seen.push(discriminator);
            result.push(path);
        }
    }
    Ok(result)
}

pub fn atom_environment(
    topology: &TopologyBlock,
    radius: usize,
    root: AtomId,
    params: &AtomEnvironmentParams,
) -> Result<AtomEnvironment, PathError> {
    // BEGIN RDKIT CPP FUNCTION findAtomEnvironmentOfRadiusN
    // RDKit✔️✔️: if (rootedAtAtom >= mol.getNumAtoms()) {
    // RDKit✔️✔️:   throw ValueErrorException("bad atom index");
    // RDKit✔️✔️: }
    topology.validate().map_err(PathError::InvalidTopology)?;
    validate_required_atom(topology, root, "root")?;
    let mut atom_distances = vec![None; topology.atoms.len()];
    atom_distances[root.index()] = Some(0);
    if radius == 0 {
        return Ok(AtomEnvironment {
            bonds: Vec::new(),
            atom_distances,
        });
    }
    // RDKit✔️✔️: std::list<std::pair<int, int>> nbrStack;
    let mut layer = VecDeque::new();
    for neighbor in topology.adjacency.neighbors_of(root.index()) {
        if params.use_hydrogens || topology.atoms[neighbor.atom_index].atomic_number() != 1 {
            layer.push_back((root.index(), neighbor.bond));
        }
    }
    let mut bonds_in = vec![false; topology.bonds.len()];
    let mut bonds = Vec::new();
    let mut completed_layers = 0usize;
    for depth in 0..radius {
        if layer.is_empty() {
            break;
        }
        let mut next_layer = VecDeque::new();
        while let Some((start_atom, bond_id)) = layer.pop_front() {
            if bonds_in[bond_id.index()] {
                continue;
            }
            bonds_in[bond_id.index()] = true;
            bonds.push(bond_id);
            let bond = &topology.bonds[bond_id.index()];
            let other = if bond.begin().index() == start_atom {
                bond.end().index()
            } else {
                bond.begin().index()
            };
            let distance = depth + 1;
            atom_distances[other] =
                Some(atom_distances[other].map_or(distance, |old| old.min(distance)));
            if depth < radius - 1 {
                for neighbor in topology.adjacency.neighbors_of(other) {
                    if !bonds_in[neighbor.bond.index()]
                        && (params.use_hydrogens
                            || topology.atoms[neighbor.atom_index].atomic_number() != 1)
                    {
                        next_layer.push_back((other, neighbor.bond));
                    }
                }
            }
        }
        layer = next_layer;
        completed_layers += 1;
    }
    // RDKit✔️✔️: if (i != radius && enforceSize) {
    // RDKit✔️✔️:   res.clear();
    // RDKit✔️✔️:   if (atomMap) { atomMap->clear(); }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION findAtomEnvironmentOfRadiusN
    if completed_layers != radius && params.enforce_radius {
        bonds.clear();
        atom_distances.fill(None);
    }
    Ok(AtomEnvironment {
        bonds,
        atom_distances,
    })
}

pub fn bond_ids_from_atom_path(
    topology: &TopologyBlock,
    atoms: &[AtomId],
) -> Result<Vec<BondId>, PathError> {
    // BEGIN RDKIT CPP FUNCTION bondListFromAtomList
    // RDKit✔️✔️: PATH_TYPE bondListFromAtomList(const ROMol &mol, const PATH_TYPE &atomIds) {
    // RDKit✔️✔️:   PATH_TYPE bids;
    // RDKit✔️✔️:   unsigned int natms = atomIds.size();
    // RDKit✔️✔️:   if (natms <= 1) { return bids; }
    topology.validate().map_err(PathError::InvalidTopology)?;
    for (position, atom) in atoms.iter().copied().enumerate() {
        if atom.index() >= topology.atoms.len() {
            return Err(PathError::AtomPathOutOfRange {
                position,
                atom,
                atom_count: topology.atoms.len(),
            });
        }
    }
    let mut result = Vec::new();
    // RDKit✔️✔️:   for (unsigned int i = 0; i < natms; i++) {
    // RDKit✔️✔️:     for (unsigned int j = i + 1; j < natms; j++) {
    // RDKit✔️✔️:       const Bond *bnd = mol.getBondBetweenAtoms(atomIds[i], atomIds[j]);
    // RDKit✔️✔️:       if (bnd) { bids.push_back(bnd->getIdx()); }
    for i in 0..atoms.len() {
        for j in i + 1..atoms.len() {
            if let Some(bond) = bond_between(topology, atoms[i].index(), atoms[j].index()) {
                result.push(bond);
            }
        }
    }
    // RDKit✔️✔️:   return bids;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION bondListFromAtomList
    Ok(result)
}

pub fn subtopology_from_path(
    topology: &TopologyBlock,
    bonds: &[BondId],
    params: &SubtopologyParams,
) -> Result<SubtopologyResult, PathError> {
    topology.validate().map_err(PathError::InvalidTopology)?;
    // BEGIN RDKIT CPP FUNCTION getSubsetInfo/copySelectedAtomsAndBonds
    // RDKit✔️✔️: } else if (options.method == SubsetMethod::BONDS) {
    // RDKit✔️✔️:   for (const auto &bond_idx : path) {
    // RDKit✔️✔️:     if (bond_idx < num_bonds) {
    // RDKit✔️✔️:       selectedBonds.set(bond_idx);
    // RDKit✔️✔️:       const auto &bnd = mol.getBondWithIdx(bond_idx);
    // RDKit✔️✔️:       selectedAtoms.set(bnd->getBeginAtomIdx());
    // RDKit✔️✔️:       selectedAtoms.set(bnd->getEndAtomIdx());
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    let mut selected_bonds = vec![false; topology.bonds.len()];
    let mut selected_atoms = vec![false; topology.atoms.len()];
    for bond_id in bonds.iter().copied() {
        let Some(bond) = topology.bonds.get(bond_id.index()) else {
            continue;
        };
        selected_bonds[bond_id.index()] = true;
        selected_atoms[bond.begin().index()] = true;
        selected_atoms[bond.end().index()] = true;
    }

    let mut atom_old_to_new = vec![None; topology.atoms.len()];
    let mut atom_new_to_old = Vec::new();
    let mut copied_atoms = Vec::new();
    // RDKit✔️✔️: for (const auto &ref_atom : reference_mol.atoms()) {
    // RDKit✔️✔️:   if (!selectedAtoms[ref_atom->getIdx()]) { continue; }
    // RDKit✔️✔️:   extracted_atom->clearComputedProps();
    for atom in &topology.atoms {
        if !selected_atoms[atom.id().index()] {
            continue;
        }
        let new_id = AtomId::new(copied_atoms.len());
        atom_old_to_new[atom.id().index()] = Some(new_id);
        atom_new_to_old.push(Some(atom.id()));
        let mut copied = atom.clone().with_id(new_id);
        copied.clear_computed_props();
        copied_atoms.push(copied);
    }

    let mut bond_old_to_new = vec![None; topology.bonds.len()];
    let mut bond_new_to_old = Vec::new();
    let mut copied_bonds = Vec::new();
    // RDKit✔️✔️: for (const auto &ref_bond : reference_mol.bonds()) {
    // RDKit✔️✔️:   if (!selectedBonds[ref_bond->getIdx()]) { continue; }
    for bond in &topology.bonds {
        if !selected_bonds[bond.id().index()] {
            continue;
        }
        let new_id = BondId::new(copied_bonds.len());
        let begin =
            atom_old_to_new[bond.begin().index()].expect("selected bond endpoints are selected");
        let end =
            atom_old_to_new[bond.end().index()].expect("selected bond endpoints are selected");
        // RDKit✔️✔️: if (atoms.size() == 2) {
        // RDKit✔️✔️:   if (map1 != atomMapping.end() && map2 != atomMapping.end()) {
        // RDKit✔️✔️:     atoms[0] = map1->second; atoms[1] = map2->second;
        // RDKit✔️✔️:   } else { atoms.clear(); }
        let stereo_atoms = bond.stereo_atoms().and_then(|[left, right]| {
            Some([
                atom_old_to_new
                    .get(left.index())
                    .and_then(|mapped| *mapped)?,
                atom_old_to_new
                    .get(right.index())
                    .and_then(|mapped| *mapped)?,
            ])
        });
        let mut copied = bond.clone();
        copied.clear_computed_props();
        if stereo_atoms.is_none() && matches!(copied.stereo(), BondStereo::Cis | BondStereo::Trans)
        {
            copied.set_stereo_atoms(None);
            copied
                .set_stereo(BondStereo::None)
                .expect("clearing stereo cannot fail");
        }
        copied = copied.remapped(new_id, begin, end, stereo_atoms);
        bond_old_to_new[bond.id().index()] = Some(new_id);
        bond_new_to_old.push(Some(bond.id()));
        copied_bonds.push(copied);
    }

    let mapping = TopologyMapping {
        atoms: AtomMapping {
            old_to_new: atom_old_to_new,
            new_to_old: atom_new_to_old,
        },
        bonds: BondMapping {
            old_to_new: bond_old_to_new,
            new_to_old: bond_new_to_old,
        },
    };
    mapping
        .validate_for_counts(
            topology.atoms.len(),
            copied_atoms.len(),
            topology.bonds.len(),
            copied_bonds.len(),
        )
        .map_err(PathError::Mapping)?;

    let substance_groups = remap_selected_substance_groups(topology, &mapping);
    let stereo_groups = remap_selected_stereo_groups(topology, &mapping);

    let subgraph = if params.copy_as_query {
        // RDKit✔️✔️: std::unique_ptr<Atom> extracted_atom{
        // RDKit✔️✔️:     options.copyAsQuery ? new QueryAtom(*ref_atom) : ref_atom->copy()};
        // RDKit✔️✔️: std::unique_ptr<Bond> extracted_bond{
        // RDKit✔️✔️:     options.copyAsQuery ? new QueryBond(*ref_bond) : ref_bond->copy()};
        let query_atoms = copied_atoms
            .into_iter()
            .map(|atom| {
                let predicate =
                    QueryNode::predicate(AtomQueryPredicate::AtomicNumber(atom.atomic_number()));
                QueryAtom::from_parts(atom, predicate)
            })
            .collect();
        let query_bonds = copied_bonds
            .into_iter()
            .map(|bond| {
                let predicate = if bond.order() == cosmolkit_model::BondOrder::Unspecified {
                    QueryNode::predicate(BondQueryPredicate::Any)
                } else {
                    QueryNode::predicate(BondQueryPredicate::Order(bond.order()))
                };
                QueryBond::from_parts(bond, predicate)
            })
            .collect();
        let graph = QueryGraph::from_parts(
            query_atoms,
            query_bonds,
            BTreeMap::new(),
            Vec::new(),
            Vec::new(),
            stereo_groups,
        )
        .map_err(PathError::QueryGraph)?;
        DetachedPathSubgraph::Query {
            graph,
            substance_groups,
        }
    } else {
        let concrete = TopologyBlock::try_from_parts(
            copied_atoms,
            copied_bonds,
            substance_groups,
            stereo_groups,
        )
        .map_err(PathError::InvalidSubsetTopology)?;
        DetachedPathSubgraph::Concrete(concrete)
    };
    // END RDKIT CPP FUNCTION getSubsetInfo/copySelectedAtomsAndBonds
    Ok(SubtopologyResult { subgraph, mapping })
}

fn validate_required_atom(
    topology: &TopologyBlock,
    atom: AtomId,
    role: &'static str,
) -> Result<(), PathError> {
    if atom.index() >= topology.atoms.len() {
        Err(PathError::AtomOutOfRange {
            role,
            atom,
            atom_count: topology.atoms.len(),
        })
    } else {
        Ok(())
    }
}

fn validate_range(lower: usize, upper: usize) -> Result<(), PathError> {
    if lower > upper {
        Err(PathError::InvalidLengthRange { lower, upper })
    } else {
        Ok(())
    }
}

fn atom_adjacency_matrix(topology: &TopologyBlock, include_hydrogens: bool) -> Vec<bool> {
    // RDKit✔️✔️: for (bondIt = mol.beginBonds(); bondIt != mol.endBonds(); bondIt++) {
    // RDKit✔️✔️:   Atom *beg = (*bondIt)->getBeginAtom();
    // RDKit✔️✔️:   Atom *end = (*bondIt)->getEndAtom();
    // RDKit✔️✔️:   if (useHs || (beg->getAtomicNum() != 1 && end->getAtomicNum() != 1)) {
    // RDKit✔️✔️:     adjMat[beg->getIdx() * dim + end->getIdx()] = 1;
    // RDKit✔️✔️:     adjMat[end->getIdx() * dim + beg->getIdx()] = 1;
    // RDKit✔️✔️:   }
    let dimension = topology.atoms.len();
    let mut adjacency = vec![false; dimension.saturating_mul(dimension)];
    for bond in &topology.bonds {
        let begin = bond.begin().index();
        let end = bond.end().index();
        if include_hydrogens
            || (topology.atoms[begin].atomic_number() != 1
                && topology.atoms[end].atomic_number() != 1)
        {
            adjacency[begin * dimension + end] = true;
            adjacency[end * dimension + begin] = true;
        }
    }
    adjacency
}

fn path_finder_helper(
    adjacency: &[bool],
    dimension: usize,
    minimum_length: usize,
    maximum_length: usize,
    root: Option<AtomId>,
    distances: Option<&[usize]>,
) -> BTreeMap<usize, Vec<Vec<usize>>> {
    // BEGIN RDKIT CPP FUNCTION pathFinderHelper
    // RDKit✔️✔️: if (rootedAtAtom < 0) {
    // RDKit✔️✔️:   for (unsigned int i = 0; i < dim; i++) {
    // RDKit✔️✔️:     PATH_TYPE tPath; tPath.push_back(i); paths.push_back(tPath);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: } else if (rootedAtAtom < static_cast<int>(dim)) {
    // RDKit✔️✔️:   PATH_TYPE tPath; tPath.push_back(rootedAtAtom); paths.push_back(tPath);
    // RDKit✔️✔️: } else { return res; }
    let mut paths = match root {
        None => (0..dimension).map(|atom| vec![atom]).collect(),
        Some(root) if root.index() < dimension => vec![vec![root.index()]],
        Some(_) => return BTreeMap::new(),
    };
    let mut result = BTreeMap::new();
    // RDKit✔️✔️: for (unsigned int length = 1; length < maxLen; length++) {
    // RDKit✔️✔️:   if (length >= minLen) { res[length] = paths; }
    // RDKit✔️✔️:   paths = extendPaths(adjMat, dim, paths, maxLen, distMat);
    // RDKit✔️✔️: }
    for length in 1..maximum_length {
        if length >= minimum_length {
            result.insert(length, paths.clone());
        }
        paths = extend_paths(adjacency, dimension, &paths, maximum_length, distances);
    }
    // RDKit✔️✔️: res[maxLen] = paths;
    // END RDKIT CPP FUNCTION pathFinderHelper
    result.insert(maximum_length, paths);
    result
}

fn extend_paths(
    adjacency: &[bool],
    dimension: usize,
    paths: &[Vec<usize>],
    allow_ring_closures: usize,
    distances: Option<&[usize]>,
) -> Vec<Vec<usize>> {
    // BEGIN RDKIT CPP FUNCTION extendPaths
    // RDKit✔️✔️: for (path = paths.begin(); path != paths.end(); ++path) {
    // RDKit✔️✔️:   unsigned int endIdx = (*path)[path->size() - 1];
    // RDKit✔️✔️:   for (unsigned int otherIdx = 0; otherIdx < dim; otherIdx++) {
    let mut result = Vec::new();
    for path in paths {
        let end = *path.last().expect("path finder never stores an empty path");
        for other in 0..dimension {
            if !adjacency[end * dimension + other] {
                continue;
            }
            // RDKit✔️✔️: if (distMat &&
            // RDKit✔️✔️:     distMat[path->front() * dim + otherIdx] - path->size() < -0.001) {
            // RDKit✔️✔️:   continue;
            // RDKit✔️✔️: }
            if distances.is_some_and(|matrix| matrix[path[0] * dimension + other] < path.len()) {
                continue;
            }
            if !path.contains(&other) {
                let mut extended = path.clone();
                extended.push(other);
                result.push(extended);
            } else if allow_ring_closures > 2 && path.len() == allow_ring_closures - 1 {
                let penultimate = path[path.len() - 2];
                if penultimate != other {
                    let mut extended = path.clone();
                    extended.push(other);
                    result.push(extended);
                }
            }
        }
    }
    // END RDKIT CPP FUNCTION extendPaths
    result
}

fn bond_neighbor_map(topology: &TopologyBlock, use_hydrogens: bool) -> BTreeMap<usize, Vec<usize>> {
    // BEGIN RDKIT CPP FUNCTION getNbrsList
    // RDKit✔️✔️: for (int i = 0; i < nAtoms; i++) {
    // RDKit✔️✔️:   const Atom *atom = mol.getAtomWithIdx(i);
    // RDKit✔️✔️:   if (useHs || atom->getAtomicNum() != 1) {
    // RDKit✔️✔️:     while (bIt1 != end) {
    // RDKit✔️✔️:       const Bond *bond1 = mol[*bIt1];
    let mut result = BTreeMap::<usize, Vec<usize>>::new();
    for atom_index in 0..topology.atoms.len() {
        if !use_hydrogens && topology.atoms[atom_index].atomic_number() == 1 {
            continue;
        }
        let atom_bonds = topology.adjacency.neighbors_of(atom_index);
        for bond1 in atom_bonds {
            if !use_hydrogens && topology.atoms[bond1.atom_index].atomic_number() == 1 {
                continue;
            }
            result.entry(bond1.bond.index()).or_default();
            for bond2 in atom_bonds {
                if bond1.bond != bond2.bond
                    && (use_hydrogens || topology.atoms[bond2.atom_index].atomic_number() != 1)
                {
                    result
                        .get_mut(&bond1.bond.index())
                        .expect("entry was inserted")
                        .push(bond2.bond.index());
                }
            }
        }
    }
    // END RDKIT CPP FUNCTION getNbrsList
    result
}

fn all_subgraphs_of_length_from_neighbors(
    topology: &TopologyBlock,
    neighbors: &BTreeMap<usize, Vec<usize>>,
    target_length: usize,
    root: Option<AtomId>,
) -> Vec<Vec<BondId>> {
    // BEGIN RDKIT CPP FUNCTION findAllSubgraphsOfLengthN
    // RDKit✔️✔️: boost::dynamic_bitset<> forbidden(mol.getNumBonds());
    // RDKit✔️✔️: for (auto nbi = nbrs.begin(); nbi != nbrs.end(); ++nbi) {
    let mut forbidden = vec![false; topology.bonds.len()];
    let mut raw_result = Vec::new();
    for (&start, adjacent) in neighbors {
        if !root_allows_bond(topology, root, start) || forbidden[start] {
            continue;
        }
        forbidden[start] = true;
        recurse_walk(
            neighbors,
            vec![start],
            adjacent.clone(),
            target_length,
            forbidden.clone(),
            &mut raw_result,
        );
    }
    // RDKit✔️✔️: return res;
    // END RDKIT CPP FUNCTION findAllSubgraphsOfLengthN
    raw_result
        .into_iter()
        .map(|row| row.into_iter().map(BondId::new).collect())
        .collect()
}

fn recurse_walk(
    neighbors: &BTreeMap<usize, Vec<usize>>,
    path: Vec<usize>,
    mut candidates: Vec<usize>,
    target_length: usize,
    mut forbidden: Vec<bool>,
    result: &mut Vec<Vec<usize>>,
) {
    // BEGIN RDKIT CPP FUNCTION recurseWalk
    // RDKit✔️✔️: if (spath.size() == targetLen) { res.push_back(spath); return; }
    // RDKit✔️✔️: if (spath.size() > targetLen) { return; }
    if path.len() == target_length {
        result.push(path);
        return;
    }
    if path.len() > target_length {
        return;
    }
    // RDKit✔️✔️: while (cands.size() != 0) {
    // RDKit✔️✔️:   int next = cands.back(); cands.pop_back();
    while let Some(next) = candidates.pop() {
        if forbidden[next] {
            continue;
        }
        forbidden[next] = true;
        let mut stack = candidates.clone();
        if let Some(next_neighbors) = neighbors.get(&next) {
            for &bond in next_neighbors {
                if !forbidden[bond] {
                    stack.push(bond);
                }
            }
        }
        let mut next_path = path.clone();
        next_path.push(next);
        recurse_walk(
            neighbors,
            next_path,
            stack,
            target_length,
            forbidden.clone(),
            result,
        );
    }
    // END RDKIT CPP FUNCTION recurseWalk
}

fn recurse_walk_range(
    neighbors: &BTreeMap<usize, Vec<usize>>,
    path: Vec<usize>,
    mut candidates: Vec<usize>,
    lower_length: usize,
    upper_length: usize,
    mut forbidden: Vec<bool>,
    result: &mut BTreeMap<usize, Vec<Vec<BondId>>>,
) {
    // BEGIN RDKIT CPP FUNCTION recurseWalkRange
    // RDKit✔️✔️: unsigned int nsize = spath.size();
    // RDKit✔️✔️: if ((nsize >= lowerLen) && (nsize <= upperLen)) {
    // RDKit✔️✔️:   res[nsize].push_back(spath);
    // RDKit✔️✔️: }
    let length = path.len();
    if length >= lower_length && length <= upper_length {
        result
            .entry(length)
            .or_default()
            .push(path.iter().copied().map(BondId::new).collect());
    }
    if length >= upper_length {
        return;
    }
    while let Some(next) = candidates.pop() {
        if forbidden[next] {
            continue;
        }
        forbidden[next] = true;
        let mut stack = candidates.clone();
        if let Some(next_neighbors) = neighbors.get(&next) {
            for &bond in next_neighbors {
                if !forbidden[bond] {
                    stack.push(bond);
                }
            }
        }
        let mut next_path = path.clone();
        next_path.push(next);
        recurse_walk_range(
            neighbors,
            next_path,
            stack,
            lower_length,
            upper_length,
            forbidden.clone(),
            result,
        );
    }
    // END RDKIT CPP FUNCTION recurseWalkRange
}

fn root_allows_bond(topology: &TopologyBlock, root: Option<AtomId>, bond: usize) -> bool {
    let Some(root) = root else {
        return true;
    };
    let bond = &topology.bonds[bond];
    bond.begin() == root || bond.end() == root
}

fn bond_between(topology: &TopologyBlock, begin: usize, end: usize) -> Option<BondId> {
    topology
        .adjacency
        .neighbors_of(begin)
        .iter()
        .find(|neighbor| neighbor.atom_index == end)
        .map(|neighbor| neighbor.bond)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct PathDiscriminator(u32, usize, usize);

fn path_discriminator(
    topology: &TopologyBlock,
    path: &[BondId],
    use_bond_orders: bool,
    extra_atom_invariants: Option<&[u32]>,
) -> Result<PathDiscriminator, PathError> {
    // BEGIN RDKIT CPP FUNCTION calcPathDiscriminators
    // RDKit✔️✔️: std::vector<int32_t> atomsUsed(mol.getNumAtoms(), -1);
    // RDKit✔️✔️: std::vector<const Atom *> atoms;
    // RDKit✔️✔️: std::vector<uint32_t> pathDegrees;
    let mut atom_positions = vec![None; topology.atoms.len()];
    let mut atoms = Vec::new();
    let mut degrees = Vec::<u32>::new();
    for bond_id in path {
        let bond = &topology.bonds[bond_id.index()];
        for atom_id in [bond.begin(), bond.end()] {
            match atom_positions[atom_id.index()] {
                Some(position) => degrees[position] += 1,
                None => {
                    atom_positions[atom_id.index()] = Some(atoms.len());
                    atoms.push(atom_id);
                    degrees.push(1);
                }
            }
        }
    }
    // RDKit✔️✔️: uint32_t invar = atom->getAtomicNum();
    // RDKit✔️✔️: hash_combine(invar, pathDegrees[i]);
    // RDKit✔️✔️: hash_combine(invar, atom->getFormalCharge());
    // RDKit✔️✔️: int deltaMass = static_cast<int>(
    // RDKit✔️✔️:     atom->getMass() - PeriodicTable::getTable()->getAtomicWeight(atom->getAtomicNum()));
    // RDKit✔️✔️: hash_combine(invar, deltaMass);
    let mut invariants = Vec::with_capacity(atoms.len());
    for (position, atom_id) in atoms.iter().copied().enumerate() {
        let atom = &topology.atoms[atom_id.index()];
        let mut invariant = u32::from(atom.atomic_number());
        hash_combine(&mut invariant, degrees[position]);
        hash_combine(&mut invariant, i32::from(atom.formal_charge()) as u32);
        let isotope_mass =
            atomic_mass(atom.element(), atom.isotope()).map_err(PathError::PeriodicTable)?;
        let average_mass = atomic_mass(atom.element(), None).map_err(PathError::PeriodicTable)?;
        hash_combine(&mut invariant, (isotope_mass - average_mass) as i32 as u32);
        if atom.is_aromatic() {
            hash_combine(&mut invariant, 1);
        }
        if let Some(extra) = extra_atom_invariants {
            hash_combine(&mut invariant, extra[atom_id.index()]);
        }
        invariants.push(invariant);
    }
    // RDKit✔️✔️: unsigned int nCycles = path.size() / 2 + 1;
    // RDKit✔️✔️: for (unsigned int cycle = 0; cycle < nCycles; ++cycle) {
    for _ in 0..path.len() / 2 + 1 {
        let mut local = vec![Vec::<u32>::new(); atoms.len()];
        for bond_id in path {
            let bond = &topology.bonds[bond_id.index()];
            let begin = atom_positions[bond.begin().index()].expect("path atom was recorded");
            let end = atom_positions[bond.end().index()].expect("path atom was recorded");
            let mut begin_value = invariants[begin];
            let mut end_value = invariants[end];
            if use_bond_orders {
                hash_combine(&mut begin_value, bond.order().rdkit_code() as u32);
                hash_combine(&mut end_value, bond.order().rdkit_code() as u32);
            }
            local[begin].push(end_value);
            local[end].push(begin_value);
        }
        for (invariant, neighbors) in invariants.iter_mut().zip(&mut local) {
            neighbors.sort_unstable();
            *invariant = hash_range(neighbors);
        }
    }
    invariants.sort_unstable();
    Ok(PathDiscriminator(
        hash_range(&invariants),
        path.len(),
        atoms.len(),
    ))
    // END RDKIT CPP FUNCTION calcPathDiscriminators
}

fn hash_combine(seed: &mut u32, value: u32) {
    // RDKit✔️✔️: seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    *seed ^= value
        .wrapping_add(0x9e37_79b9)
        .wrapping_add((*seed).wrapping_shl(6))
        .wrapping_add(*seed >> 2);
}

fn hash_range(values: &[u32]) -> u32 {
    // RDKit✔️✔️: std::hash_result_t seed = 0;
    // RDKit✔️✔️: for (; first != last; ++first) { hash_combine(seed, *first); }
    let mut seed = 0;
    for &value in values {
        hash_combine(&mut seed, value);
    }
    seed
}

fn remap_selected_substance_groups(
    topology: &TopologyBlock,
    mapping: &TopologyMapping,
) -> Vec<SubstanceGroup> {
    let mut selected = topology
        .substance_groups
        .iter()
        .map(|group| {
            group.can_remap_without_parent(&mapping.atoms.old_to_new, &mapping.bonds.old_to_new)
        })
        .collect::<Vec<_>>();
    loop {
        let mut changed = false;
        for (index, group) in topology.substance_groups.iter().enumerate() {
            if selected[index]
                && group
                    .parent()
                    .is_some_and(|parent| !selected.get(parent.index()).copied().unwrap_or(false))
            {
                selected[index] = false;
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
    let mut sgroup_map = vec![None; topology.substance_groups.len()];
    let mut next = 0;
    for (old, keep) in selected.iter().copied().enumerate() {
        if keep {
            sgroup_map[old] = Some(SubstanceGroupId::new(next));
            next += 1;
        }
    }
    topology
        .substance_groups
        .iter()
        .enumerate()
        .filter(|(index, _)| selected[*index])
        .filter_map(|(index, group)| {
            group.remapped(
                sgroup_map[index].expect("selected SGroup has a new id"),
                &mapping.atoms.old_to_new,
                &mapping.bonds.old_to_new,
                &sgroup_map,
            )
        })
        .collect()
}

fn remap_selected_stereo_groups(
    topology: &TopologyBlock,
    mapping: &TopologyMapping,
) -> Vec<StereoGroup> {
    topology
        .stereo_groups
        .iter()
        .filter_map(|group| {
            // RDKit✔️✔️: return objects.empty() ||
            // RDKit✔️✔️:        std::any_of(objects.begin(), objects.end(), [&](auto &object) {
            // RDKit✔️✔️:          return selected_indices[object->getIdx()];
            // RDKit✔️✔️:        });
            let atom_side_selected = group.atoms().is_empty()
                || group
                    .atoms()
                    .iter()
                    .any(|atom| mapping.atoms.old_to_new[atom.index()].is_some());
            let bond_side_selected = group.bonds().is_empty()
                || group
                    .bonds()
                    .iter()
                    .any(|bond| mapping.bonds.old_to_new[bond.index()].is_some());
            if !atom_side_selected || !bond_side_selected {
                return None;
            }
            let atoms = group
                .atoms()
                .iter()
                .filter_map(|atom| mapping.atoms.old_to_new[atom.index()])
                .collect();
            let bonds = group
                .bonds()
                .iter()
                .filter_map(|bond| mapping.bonds.old_to_new[bond.index()])
                .collect();
            let remapped = StereoGroup::new(group.kind(), atoms, bonds);
            Some(match group.id() {
                Some(id) => remapped.with_id(id),
                None => remapped,
            })
        })
        .collect()
}
