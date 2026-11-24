use std::collections::BTreeMap;

use crate::{AtomId, BondId, GraphPath, Molecule, PathError};
use cosmolkit_core as core;

/// Selects the length unit and index kind for linear path enumeration.
#[derive(Clone, Debug, Default, PartialEq, Eq, cosmolkit_macros::Configuration)]
#[configuration(scalar)]
pub struct PathEnumerationParams {
    /// Count and return atom indices instead of bond indices.
    pub use_atom_indices: bool,
    /// Include stored hydrogen atoms when enumerating paths.
    pub use_hydrogens: bool,
    /// Start at this source atom, or search the whole graph when absent.
    pub rooted_at_atom: Option<usize>,
    /// Keep only paths that are shortest between their endpoints.
    pub only_shortest_paths: bool,
}

/// Selection for enumerating connected sets of bonds, including branches.
#[derive(Clone, Debug, Default, PartialEq, Eq, cosmolkit_macros::Configuration)]
#[configuration(scalar)]
pub struct SubgraphEnumerationParams {
    /// Include bonds to stored hydrogen atoms.
    pub use_hydrogens: bool,
    /// Require the connected bond set to contain this source atom when specified.
    pub rooted_at_atom: Option<usize>,
}

/// Selection and equivalence rules for unique connected bond sets.
#[derive(Clone, Debug, PartialEq, Eq, cosmolkit_macros::Configuration)]
#[configuration(scalar)]
pub struct UniqueSubgraphSearchParams {
    /// Include bonds to stored hydrogen atoms.
    pub use_hydrogens: bool,
    /// Restrict enumeration to this source atom when specified.
    pub rooted_at_atom: Option<usize>,
    /// Distinguish bond orders when choosing unique representatives.
    pub use_bond_orders: bool,
}

impl Default for UniqueSubgraphSearchParams {
    fn default() -> Self {
        Self {
            use_hydrogens: false,
            rooted_at_atom: None,
            use_bond_orders: true,
        }
    }
}

/// Hydrogen and radius rules for a rooted atom environment.
#[derive(Clone, Debug, PartialEq, Eq, cosmolkit_macros::Configuration)]
#[configuration(scalar)]
pub struct AtomEnvironmentParams {
    /// Include stored hydrogen neighbors in the environment.
    pub use_hydrogens: bool,
    /// Return an empty environment if traversal does not reach the requested radius.
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

/// Source bonds and the shortest distance of each reached atom from the root.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AtomEnvironment {
    bonds: Vec<BondId>,
    atom_distances: BTreeMap<AtomId, usize>,
}

impl AtomEnvironment {
    pub fn bonds(&self) -> &[BondId] {
        &self.bonds
    }
    pub fn atom_distances(&self) -> &BTreeMap<AtomId, usize> {
        &self.atom_distances
    }
}

impl Molecule {
    fn traversal_atom(&self, atom: AtomId, role: &'static str) -> Result<AtomId, PathError> {
        if atom.index() < self.num_atoms() {
            Ok(atom)
        } else {
            Err(PathError::AtomOutOfRange {
                role,
                atom,
                atom_count: self.num_atoms(),
            })
        }
    }

    fn traversal_root(&self, index: Option<usize>) -> Result<Option<AtomId>, PathError> {
        index
            .map(|i| self.traversal_atom(AtomId::new(i), "root"))
            .transpose()
    }

    pub fn bond_between_atoms(
        &self,
        begin: AtomId,
        end: AtomId,
    ) -> Result<Option<BondId>, PathError> {
        self.traversal_atom(begin, "begin")?;
        self.traversal_atom(end, "end")?;
        Ok(self
            .topology()
            .adjacency
            .neighbors_of(begin.index())
            .iter()
            .find(|n| n.atom_index == end.index())
            .map(|n| n.bond))
    }

    /// Atom indices from begin to end, inclusive; empty if they are disconnected.
    pub fn shortest_path(&self, begin: AtomId, end: AtomId) -> Result<Vec<AtomId>, PathError> {
        core::shortest_path(self.topology(), begin, end)
    }

    pub fn paths_of_length(&self, length: usize) -> Result<Vec<GraphPath>, PathError> {
        self.paths_of_length_with_params(length, &Default::default())
    }

    pub fn paths_of_length_with_params(
        &self,
        length: usize,
        params: &PathEnumerationParams,
    ) -> Result<Vec<GraphPath>, PathError> {
        u32::try_from(length).map_err(|_| PathError::LengthOverflow { length })?;
        let representation = if params.use_atom_indices {
            core::PathRepresentation::Atoms
        } else {
            core::PathRepresentation::Bonds
        };
        core::all_paths_of_length(
            self.topology(),
            length,
            &core::PathSearchParams {
                representation,
                use_hydrogens: params.use_hydrogens,
                rooted_at_atom: self.traversal_root(params.rooted_at_atom)?,
                only_shortest_paths: params.only_shortest_paths,
                ..Default::default()
            },
        )
    }

    pub fn subgraphs_of_length(&self, length: usize) -> Result<Vec<Vec<BondId>>, PathError> {
        self.subgraphs_of_length_with_params(length, &Default::default())
    }

    pub fn subgraphs_of_length_with_params(
        &self,
        length: usize,
        params: &SubgraphEnumerationParams,
    ) -> Result<Vec<Vec<BondId>>, PathError> {
        u32::try_from(length).map_err(|_| PathError::LengthOverflow { length })?;
        core::all_subgraphs_of_length(self.topology(), length, &self.subgraph_search(params)?)
    }

    /// Element `i` contains connected sets of `min + i` bonds, including empty groups.
    pub fn subgraphs_in_range(
        &self,
        min: usize,
        max: usize,
    ) -> Result<Vec<Vec<Vec<BondId>>>, PathError> {
        self.subgraphs_in_range_with_params(min, max, &Default::default())
    }

    pub fn subgraphs_in_range_with_params(
        &self,
        min: usize,
        max: usize,
        params: &SubgraphEnumerationParams,
    ) -> Result<Vec<Vec<Vec<BondId>>>, PathError> {
        u32::try_from(max).map_err(|_| PathError::LengthOverflow { length: max })?;
        Ok(
            core::all_subgraphs_in_range(
                self.topology(),
                min,
                max,
                &self.subgraph_search(params)?,
            )?
            .into_values()
            .collect(),
        )
    }

    fn subgraph_search(
        &self,
        params: &SubgraphEnumerationParams,
    ) -> Result<core::SubgraphSearchParams<'static>, PathError> {
        Ok(core::SubgraphSearchParams {
            use_hydrogens: params.use_hydrogens,
            rooted_at_atom: self.traversal_root(params.rooted_at_atom)?,
            ..Default::default()
        })
    }

    pub fn unique_subgraphs_of_length(&self, length: usize) -> Result<Vec<Vec<BondId>>, PathError> {
        self.unique_subgraphs_of_length_with_params(length, &Default::default())
    }

    pub fn unique_subgraphs_of_length_with_params(
        &self,
        length: usize,
        params: &UniqueSubgraphSearchParams,
    ) -> Result<Vec<Vec<BondId>>, PathError> {
        u32::try_from(length).map_err(|_| PathError::LengthOverflow { length })?;
        core::unique_subgraphs_of_length(
            self.topology(),
            length,
            &core::UniqueSubgraphParams {
                use_hydrogens: params.use_hydrogens,
                rooted_at_atom: self.traversal_root(params.rooted_at_atom)?,
                use_bond_orders: params.use_bond_orders,
                ..Default::default()
            },
        )
    }

    pub fn atom_environment(
        &self,
        radius: usize,
        root: AtomId,
    ) -> Result<AtomEnvironment, PathError> {
        self.atom_environment_with_params(radius, root, &Default::default())
    }

    pub fn atom_environment_with_params(
        &self,
        radius: usize,
        root: AtomId,
        params: &AtomEnvironmentParams,
    ) -> Result<AtomEnvironment, PathError> {
        u32::try_from(radius).map_err(|_| PathError::LengthOverflow { length: radius })?;
        let env = core::atom_environment(
            self.topology(),
            radius,
            root,
            &core::AtomEnvironmentParams {
                use_hydrogens: params.use_hydrogens,
                enforce_radius: params.enforce_radius,
            },
        )?;
        let atom_distances = env
            .atom_distances
            .into_iter()
            .enumerate()
            .filter_map(|(i, d)| d.map(|d| (AtomId::new(i), d)))
            .collect();
        Ok(AtomEnvironment {
            bonds: env.bonds,
            atom_distances,
        })
    }
}
