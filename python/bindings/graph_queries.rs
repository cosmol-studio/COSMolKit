use crate::{
    arguments::{Argument, configuration_arguments},
    canonical_errors::CanonicalRankError,
    matrix_errors::MatrixError,
    native::Molecule,
    native_numerical_values::DenseMatrix,
    params::{
        AtomEnvironmentParams, AtomRankParams, PathEnumerationParams, SubgraphEnumerationParams,
        UniqueSubgraphSearchParams,
    },
    path_errors::PathError,
};
use cosmolkit as ck;
use pyo3::prelude::*;
use std::collections::BTreeMap;
crate::configuration::type_projection! { value(eq)
    /// Bonds reached around a source atom and the shortest distance to each reached atom.
    pub struct AtomEnvironment { pub(crate) inner:ck::AtomEnvironment }
}
#[pymethods]
impl AtomEnvironment {
    /// Bonds reached in traversal order.
    #[getter]
    fn bonds(&self) -> Vec<usize> {
        self.inner.bonds().iter().map(|i| i.index()).collect()
    }
    /// Shortest distance from the root for each reached atom index.
    #[getter]
    fn atom_distances(&self) -> BTreeMap<usize, usize> {
        self.inner
            .atom_distances()
            .iter()
            .map(|(i, d)| (i.index(), *d))
            .collect()
    }
}
fn path_indices(path: ck::GraphPath) -> Vec<usize> {
    match path {
        ck::GraphPath::Atoms(v) => v.into_iter().map(|i| i.index()).collect(),
        ck::GraphPath::Bonds(v) => bond_indices(v),
    }
}
fn bond_indices(v: Vec<ck::BondId>) -> Vec<usize> {
    v.into_iter().map(|i| i.index()).collect()
}
#[pymethods]
impl Molecule {
    /// Canonical ranks in source atom order. Sort indices by rank to obtain a permutation.
    #[pyo3(signature = (params=Argument::Missing, *, break_ties=Argument::Missing, include_chirality=Argument::Missing, include_isotopes=Argument::Missing, include_atom_maps=Argument::Missing, include_chiral_presence=Argument::Missing))]
    fn canonical_atom_ranks(
        &self,
        py: Python<'_>,
        params: Argument<Option<PyRef<'_, AtomRankParams>>>,
        break_ties: Argument<Option<bool>>,
        include_chirality: Argument<Option<bool>>,
        include_isotopes: Argument<Option<bool>>,
        include_atom_maps: Argument<Option<bool>>,
        include_chiral_presence: Argument<Option<bool>>,
    ) -> PyResult<Vec<usize>> {
        let has_keywords = break_ties.is_supplied()
            || include_chirality.is_supplied()
            || include_isotopes.is_supplied()
            || include_atom_maps.is_supplied()
            || include_chiral_presence.is_supplied();
        configuration_arguments!(params; break_ties, include_chirality, include_isotopes, include_atom_maps, include_chiral_presence);
        if has_keywords {
            self._canonical_atom_ranks_keywords(
                break_ties,
                include_chirality,
                include_isotopes,
                include_atom_maps,
                include_chiral_presence,
            )
        } else {
            self._canonical_atom_ranks_params(py, params.as_deref())
        }
    }

    #[doc = r#"{"role":"overload","target":"Molecule.canonical_atom_ranks"}"#]
    #[pyo3(signature=(params=None))]
    fn _canonical_atom_ranks_params(
        &self,
        py: Python<'_>,
        params: Option<&AtomRankParams>,
    ) -> PyResult<Vec<usize>> {
        let p = params
            .map(|p| p.to_rust(py))
            .transpose()?
            .unwrap_or_default();
        self.inner
            .canonical_atom_ranks_with_params(&p)
            .map_err(CanonicalRankError::from_rust)
    }

    #[doc = r#"{"role":"overload","target":"Molecule.canonical_atom_ranks"}"#]
    #[pyo3(signature = (*, break_ties=None, include_chirality=None, include_isotopes=None, include_atom_maps=None, include_chiral_presence=None))]
    fn _canonical_atom_ranks_keywords(
        &self,
        break_ties: Option<bool>,
        include_chirality: Option<bool>,
        include_isotopes: Option<bool>,
        include_atom_maps: Option<bool>,
        include_chiral_presence: Option<bool>,
    ) -> PyResult<Vec<usize>> {
        let defaults = ck::AtomRankParams::default();
        let p = ck::AtomRankParams {
            break_ties: break_ties.unwrap_or(defaults.break_ties),
            include_chirality: include_chirality.unwrap_or(defaults.include_chirality),
            include_isotopes: include_isotopes.unwrap_or(defaults.include_isotopes),
            include_atom_maps: include_atom_maps.unwrap_or(defaults.include_atom_maps),
            include_chiral_presence: include_chiral_presence
                .unwrap_or(defaults.include_chiral_presence),
        };
        self.inner
            .canonical_atom_ranks_with_params(&p)
            .map_err(CanonicalRankError::from_rust)
    }

    /// Canonical ranks in source atom order. Sort indices by rank to obtain a permutation. Uses the supplied configuration.
    #[pyo3(signature = (params))]
    fn canonical_atom_ranks_with_params(
        &self,
        py: Python<'_>,
        params: &AtomRankParams,
    ) -> PyResult<Vec<usize>> {
        self._canonical_atom_ranks_params(py, Some(params))
    }

    /// Rank selected atoms using only the specified bonds. Unselected atoms have rank -1.
    #[pyo3(signature = (atoms, bonds=None, atom_symbols=None, params=Argument::Missing, *, break_ties=Argument::Missing, include_chirality=Argument::Missing, include_isotopes=Argument::Missing, include_atom_maps=Argument::Missing, include_chiral_presence=Argument::Missing))]
    fn canonical_atom_ranks_in_fragment(
        &self,
        py: Python<'_>,
        atoms: Vec<usize>,
        bonds: Option<Vec<usize>>,
        atom_symbols: Option<Vec<String>>,
        params: Argument<Option<PyRef<'_, AtomRankParams>>>,
        break_ties: Argument<Option<bool>>,
        include_chirality: Argument<Option<bool>>,
        include_isotopes: Argument<Option<bool>>,
        include_atom_maps: Argument<Option<bool>>,
        include_chiral_presence: Argument<Option<bool>>,
    ) -> PyResult<Vec<i64>> {
        let has_keywords = break_ties.is_supplied()
            || include_chirality.is_supplied()
            || include_isotopes.is_supplied()
            || include_atom_maps.is_supplied()
            || include_chiral_presence.is_supplied();
        configuration_arguments!(params; break_ties, include_chirality, include_isotopes, include_atom_maps, include_chiral_presence);
        if has_keywords {
            self._canonical_atom_ranks_in_fragment_keywords(
                atoms,
                bonds,
                atom_symbols,
                break_ties,
                include_chirality,
                include_isotopes,
                include_atom_maps,
                include_chiral_presence,
            )
        } else {
            self._canonical_atom_ranks_in_fragment_params(
                py,
                atoms,
                bonds,
                atom_symbols,
                params.as_deref(),
            )
        }
    }

    #[doc = r#"{"role":"overload","target":"Molecule.canonical_atom_ranks_in_fragment"}"#]
    #[pyo3(signature=(atoms, bonds=None, atom_symbols=None, params=None))]
    fn _canonical_atom_ranks_in_fragment_params(
        &self,
        py: Python<'_>,
        atoms: Vec<usize>,
        bonds: Option<Vec<usize>>,
        atom_symbols: Option<Vec<String>>,
        params: Option<&AtomRankParams>,
    ) -> PyResult<Vec<i64>> {
        let p = params
            .map(|p| p.to_rust(py))
            .transpose()?
            .unwrap_or_default();
        self.inner
            .canonical_atom_ranks_in_fragment_with_params(
                &atoms.into_iter().map(ck::AtomId::new).collect::<Vec<_>>(),
                &bonds
                    .unwrap_or_default()
                    .into_iter()
                    .map(ck::BondId::new)
                    .collect::<Vec<_>>(),
                atom_symbols.as_deref(),
                &p,
            )
            .map_err(CanonicalRankError::from_rust)
    }

    #[doc = r#"{"role":"overload","target":"Molecule.canonical_atom_ranks_in_fragment"}"#]
    #[pyo3(signature = (atoms, bonds=None, atom_symbols=None, *, break_ties=None, include_chirality=None, include_isotopes=None, include_atom_maps=None, include_chiral_presence=None))]
    fn _canonical_atom_ranks_in_fragment_keywords(
        &self,
        atoms: Vec<usize>,
        bonds: Option<Vec<usize>>,
        atom_symbols: Option<Vec<String>>,
        break_ties: Option<bool>,
        include_chirality: Option<bool>,
        include_isotopes: Option<bool>,
        include_atom_maps: Option<bool>,
        include_chiral_presence: Option<bool>,
    ) -> PyResult<Vec<i64>> {
        let defaults = ck::AtomRankParams::default();
        let p = ck::AtomRankParams {
            break_ties: break_ties.unwrap_or(defaults.break_ties),
            include_chirality: include_chirality.unwrap_or(defaults.include_chirality),
            include_isotopes: include_isotopes.unwrap_or(defaults.include_isotopes),
            include_atom_maps: include_atom_maps.unwrap_or(defaults.include_atom_maps),
            include_chiral_presence: include_chiral_presence
                .unwrap_or(defaults.include_chiral_presence),
        };
        self.inner
            .canonical_atom_ranks_in_fragment_with_params(
                &atoms.into_iter().map(ck::AtomId::new).collect::<Vec<_>>(),
                &bonds
                    .unwrap_or_default()
                    .into_iter()
                    .map(ck::BondId::new)
                    .collect::<Vec<_>>(),
                atom_symbols.as_deref(),
                &p,
            )
            .map_err(CanonicalRankError::from_rust)
    }

    /// Rank selected atoms using only the specified bonds. Unselected atoms have rank -1. Uses the supplied configuration.
    #[pyo3(signature = (atoms, bonds=None, atom_symbols=None, *, params))]
    fn canonical_atom_ranks_in_fragment_with_params(
        &self,
        py: Python<'_>,
        atoms: Vec<usize>,
        bonds: Option<Vec<usize>>,
        atom_symbols: Option<Vec<String>>,
        params: &AtomRankParams,
    ) -> PyResult<Vec<i64>> {
        self._canonical_atom_ranks_in_fragment_params(py, atoms, bonds, atom_symbols, Some(params))
    }

    /// Enumerate linear paths. Length counts bonds by default, or atoms in atom-index mode.
    #[pyo3(signature = (length, params=Argument::Missing, *, use_atom_indices=Argument::Missing, use_hydrogens=Argument::Missing, rooted_at_atom=Argument::Missing, only_shortest_paths=Argument::Missing))]
    fn paths_of_length(
        &self,
        py: Python<'_>,
        length: usize,
        params: Argument<Option<PyRef<'_, PathEnumerationParams>>>,
        use_atom_indices: Argument<Option<bool>>,
        use_hydrogens: Argument<Option<bool>>,
        rooted_at_atom: Argument<Option<usize>>,
        only_shortest_paths: Argument<Option<bool>>,
    ) -> PyResult<Vec<Vec<usize>>> {
        let has_keywords = use_atom_indices.is_supplied()
            || use_hydrogens.is_supplied()
            || rooted_at_atom.is_supplied()
            || only_shortest_paths.is_supplied();
        configuration_arguments!(params; use_atom_indices, use_hydrogens, rooted_at_atom, only_shortest_paths);
        if has_keywords {
            self._paths_of_length_keywords(
                length,
                use_atom_indices,
                use_hydrogens,
                rooted_at_atom,
                only_shortest_paths,
            )
        } else {
            self._paths_of_length_params(py, length, params.as_deref())
        }
    }

    #[doc = r#"{"role":"overload","target":"Molecule.paths_of_length"}"#]
    #[pyo3(signature=(length, params=None))]
    fn _paths_of_length_params(
        &self,
        py: Python<'_>,
        length: usize,
        params: Option<&PathEnumerationParams>,
    ) -> PyResult<Vec<Vec<usize>>> {
        let p = params
            .map(|p| p.to_rust(py))
            .transpose()?
            .unwrap_or_default();
        self.inner
            .paths_of_length_with_params(length, &p)
            .map_err(PathError::from_rust)
            .map(|v| v.into_iter().map(path_indices).collect())
    }

    #[doc = r#"{"role":"overload","target":"Molecule.paths_of_length"}"#]
    #[pyo3(signature = (length, *, use_atom_indices=None, use_hydrogens=None, rooted_at_atom=None, only_shortest_paths=None))]
    fn _paths_of_length_keywords(
        &self,
        length: usize,
        use_atom_indices: Option<bool>,
        use_hydrogens: Option<bool>,
        rooted_at_atom: Option<usize>,
        only_shortest_paths: Option<bool>,
    ) -> PyResult<Vec<Vec<usize>>> {
        let defaults = ck::PathEnumerationParams::default();
        let p = ck::PathEnumerationParams {
            use_atom_indices: use_atom_indices.unwrap_or(defaults.use_atom_indices),
            use_hydrogens: use_hydrogens.unwrap_or(defaults.use_hydrogens),
            rooted_at_atom: rooted_at_atom,
            only_shortest_paths: only_shortest_paths.unwrap_or(defaults.only_shortest_paths),
        };
        self.inner
            .paths_of_length_with_params(length, &p)
            .map_err(PathError::from_rust)
            .map(|v| v.into_iter().map(path_indices).collect())
    }

    /// Enumerate linear paths. Length counts bonds by default, or atoms in atom-index mode. Uses the supplied configuration.
    #[pyo3(signature = (length, params))]
    fn paths_of_length_with_params(
        &self,
        py: Python<'_>,
        length: usize,
        params: &PathEnumerationParams,
    ) -> PyResult<Vec<Vec<usize>>> {
        self._paths_of_length_params(py, length, Some(params))
    }

    /// Enumerate connected sets of the specified number of bonds, including branched sets.
    #[pyo3(signature = (length, params=Argument::Missing, *, use_hydrogens=Argument::Missing, rooted_at_atom=Argument::Missing))]
    fn subgraphs_of_length(
        &self,
        py: Python<'_>,
        length: usize,
        params: Argument<Option<PyRef<'_, SubgraphEnumerationParams>>>,
        use_hydrogens: Argument<Option<bool>>,
        rooted_at_atom: Argument<Option<usize>>,
    ) -> PyResult<Vec<Vec<usize>>> {
        let has_keywords = use_hydrogens.is_supplied() || rooted_at_atom.is_supplied();
        configuration_arguments!(params; use_hydrogens, rooted_at_atom);
        if has_keywords {
            self._subgraphs_of_length_keywords(length, use_hydrogens, rooted_at_atom)
        } else {
            self._subgraphs_of_length_params(py, length, params.as_deref())
        }
    }

    #[doc = r#"{"role":"overload","target":"Molecule.subgraphs_of_length"}"#]
    #[pyo3(signature=(length, params=None))]
    fn _subgraphs_of_length_params(
        &self,
        py: Python<'_>,
        length: usize,
        params: Option<&SubgraphEnumerationParams>,
    ) -> PyResult<Vec<Vec<usize>>> {
        let p = params
            .map(|p| p.to_rust(py))
            .transpose()?
            .unwrap_or_default();
        self.inner
            .subgraphs_of_length_with_params(length, &p)
            .map_err(PathError::from_rust)
            .map(|v| v.into_iter().map(bond_indices).collect())
    }

    #[doc = r#"{"role":"overload","target":"Molecule.subgraphs_of_length"}"#]
    #[pyo3(signature = (length, *, use_hydrogens=None, rooted_at_atom=None))]
    fn _subgraphs_of_length_keywords(
        &self,
        length: usize,
        use_hydrogens: Option<bool>,
        rooted_at_atom: Option<usize>,
    ) -> PyResult<Vec<Vec<usize>>> {
        let defaults = ck::SubgraphEnumerationParams::default();
        let p = ck::SubgraphEnumerationParams {
            use_hydrogens: use_hydrogens.unwrap_or(defaults.use_hydrogens),
            rooted_at_atom: rooted_at_atom,
        };
        self.inner
            .subgraphs_of_length_with_params(length, &p)
            .map_err(PathError::from_rust)
            .map(|v| v.into_iter().map(bond_indices).collect())
    }

    /// Enumerate connected sets of the specified number of bonds, including branched sets. Uses the supplied configuration.
    #[pyo3(signature = (length, params))]
    fn subgraphs_of_length_with_params(
        &self,
        py: Python<'_>,
        length: usize,
        params: &SubgraphEnumerationParams,
    ) -> PyResult<Vec<Vec<usize>>> {
        self._subgraphs_of_length_params(py, length, Some(params))
    }

    /// Enumerate connected bond sets for each length from min to max, inclusive.
    #[pyo3(signature = (min, max, params=Argument::Missing, *, use_hydrogens=Argument::Missing, rooted_at_atom=Argument::Missing))]
    fn subgraphs_in_range(
        &self,
        py: Python<'_>,
        min: usize,
        max: usize,
        params: Argument<Option<PyRef<'_, SubgraphEnumerationParams>>>,
        use_hydrogens: Argument<Option<bool>>,
        rooted_at_atom: Argument<Option<usize>>,
    ) -> PyResult<Vec<Vec<Vec<usize>>>> {
        let has_keywords = use_hydrogens.is_supplied() || rooted_at_atom.is_supplied();
        configuration_arguments!(params; use_hydrogens, rooted_at_atom);
        if has_keywords {
            self._subgraphs_in_range_keywords(min, max, use_hydrogens, rooted_at_atom)
        } else {
            self._subgraphs_in_range_params(py, min, max, params.as_deref())
        }
    }

    #[doc = r#"{"role":"overload","target":"Molecule.subgraphs_in_range"}"#]
    #[pyo3(signature=(min, max, params=None))]
    fn _subgraphs_in_range_params(
        &self,
        py: Python<'_>,
        min: usize,
        max: usize,
        params: Option<&SubgraphEnumerationParams>,
    ) -> PyResult<Vec<Vec<Vec<usize>>>> {
        let p = params
            .map(|p| p.to_rust(py))
            .transpose()?
            .unwrap_or_default();
        self.inner
            .subgraphs_in_range_with_params(min, max, &p)
            .map_err(PathError::from_rust)
            .map(|v| {
                v.into_iter()
                    .map(|group| group.into_iter().map(bond_indices).collect())
                    .collect()
            })
    }

    #[doc = r#"{"role":"overload","target":"Molecule.subgraphs_in_range"}"#]
    #[pyo3(signature = (min, max, *, use_hydrogens=None, rooted_at_atom=None))]
    fn _subgraphs_in_range_keywords(
        &self,
        min: usize,
        max: usize,
        use_hydrogens: Option<bool>,
        rooted_at_atom: Option<usize>,
    ) -> PyResult<Vec<Vec<Vec<usize>>>> {
        let defaults = ck::SubgraphEnumerationParams::default();
        let p = ck::SubgraphEnumerationParams {
            use_hydrogens: use_hydrogens.unwrap_or(defaults.use_hydrogens),
            rooted_at_atom: rooted_at_atom,
        };
        self.inner
            .subgraphs_in_range_with_params(min, max, &p)
            .map_err(PathError::from_rust)
            .map(|v| {
                v.into_iter()
                    .map(|group| group.into_iter().map(bond_indices).collect())
                    .collect()
            })
    }

    /// Enumerate connected bond sets for each length from min to max, inclusive. Uses the supplied configuration.
    #[pyo3(signature = (min, max, params))]
    fn subgraphs_in_range_with_params(
        &self,
        py: Python<'_>,
        min: usize,
        max: usize,
        params: &SubgraphEnumerationParams,
    ) -> PyResult<Vec<Vec<Vec<usize>>>> {
        self._subgraphs_in_range_params(py, min, max, Some(params))
    }

    /// Return the first representative of each distinct connected bond subgraph.
    #[pyo3(signature = (length, params=Argument::Missing, *, use_hydrogens=Argument::Missing, rooted_at_atom=Argument::Missing, use_bond_orders=Argument::Missing))]
    fn unique_subgraphs_of_length(
        &self,
        py: Python<'_>,
        length: usize,
        params: Argument<Option<PyRef<'_, UniqueSubgraphSearchParams>>>,
        use_hydrogens: Argument<Option<bool>>,
        rooted_at_atom: Argument<Option<usize>>,
        use_bond_orders: Argument<Option<bool>>,
    ) -> PyResult<Vec<Vec<usize>>> {
        let has_keywords = use_hydrogens.is_supplied()
            || rooted_at_atom.is_supplied()
            || use_bond_orders.is_supplied();
        configuration_arguments!(params; use_hydrogens, rooted_at_atom, use_bond_orders);
        if has_keywords {
            self._unique_subgraphs_of_length_keywords(
                length,
                use_hydrogens,
                rooted_at_atom,
                use_bond_orders,
            )
        } else {
            self._unique_subgraphs_of_length_params(py, length, params.as_deref())
        }
    }

    #[doc = r#"{"role":"overload","target":"Molecule.unique_subgraphs_of_length"}"#]
    #[pyo3(signature=(length, params=None))]
    fn _unique_subgraphs_of_length_params(
        &self,
        py: Python<'_>,
        length: usize,
        params: Option<&UniqueSubgraphSearchParams>,
    ) -> PyResult<Vec<Vec<usize>>> {
        let p = params
            .map(|p| p.to_rust(py))
            .transpose()?
            .unwrap_or_default();
        self.inner
            .unique_subgraphs_of_length_with_params(length, &p)
            .map_err(PathError::from_rust)
            .map(|v| v.into_iter().map(bond_indices).collect())
    }

    #[doc = r#"{"role":"overload","target":"Molecule.unique_subgraphs_of_length"}"#]
    #[pyo3(signature = (length, *, use_hydrogens=None, rooted_at_atom=None, use_bond_orders=None))]
    fn _unique_subgraphs_of_length_keywords(
        &self,
        length: usize,
        use_hydrogens: Option<bool>,
        rooted_at_atom: Option<usize>,
        use_bond_orders: Option<bool>,
    ) -> PyResult<Vec<Vec<usize>>> {
        let defaults = ck::UniqueSubgraphSearchParams::default();
        let p = ck::UniqueSubgraphSearchParams {
            use_hydrogens: use_hydrogens.unwrap_or(defaults.use_hydrogens),
            rooted_at_atom: rooted_at_atom,
            use_bond_orders: use_bond_orders.unwrap_or(defaults.use_bond_orders),
        };
        self.inner
            .unique_subgraphs_of_length_with_params(length, &p)
            .map_err(PathError::from_rust)
            .map(|v| v.into_iter().map(bond_indices).collect())
    }

    /// Return the first representative of each distinct connected bond subgraph. Uses the supplied configuration.
    #[pyo3(signature = (length, params))]
    fn unique_subgraphs_of_length_with_params(
        &self,
        py: Python<'_>,
        length: usize,
        params: &UniqueSubgraphSearchParams,
    ) -> PyResult<Vec<Vec<usize>>> {
        self._unique_subgraphs_of_length_params(py, length, Some(params))
    }

    /// Return source bonds and distances to atoms in a rooted environment.
    #[pyo3(signature = (radius, root, params=Argument::Missing, *, use_hydrogens=Argument::Missing, enforce_radius=Argument::Missing))]
    fn atom_environment(
        &self,
        py: Python<'_>,
        radius: usize,
        root: usize,
        params: Argument<Option<PyRef<'_, AtomEnvironmentParams>>>,
        use_hydrogens: Argument<Option<bool>>,
        enforce_radius: Argument<Option<bool>>,
    ) -> PyResult<AtomEnvironment> {
        let has_keywords = use_hydrogens.is_supplied() || enforce_radius.is_supplied();
        configuration_arguments!(params; use_hydrogens, enforce_radius);
        if has_keywords {
            self._atom_environment_keywords(radius, root, use_hydrogens, enforce_radius)
        } else {
            self._atom_environment_params(py, radius, root, params.as_deref())
        }
    }

    #[doc = r#"{"role":"overload","target":"Molecule.atom_environment"}"#]
    #[pyo3(signature=(radius, root, params=None))]
    fn _atom_environment_params(
        &self,
        py: Python<'_>,
        radius: usize,
        root: usize,
        params: Option<&AtomEnvironmentParams>,
    ) -> PyResult<AtomEnvironment> {
        let p = params
            .map(|p| p.to_rust(py))
            .transpose()?
            .unwrap_or_default();
        self.inner
            .atom_environment_with_params(radius, ck::AtomId::new(root), &p)
            .map_err(PathError::from_rust)
            .map(|v| AtomEnvironment { inner: v })
    }

    #[doc = r#"{"role":"overload","target":"Molecule.atom_environment"}"#]
    #[pyo3(signature = (radius, root, *, use_hydrogens=None, enforce_radius=None))]
    fn _atom_environment_keywords(
        &self,
        radius: usize,
        root: usize,
        use_hydrogens: Option<bool>,
        enforce_radius: Option<bool>,
    ) -> PyResult<AtomEnvironment> {
        let defaults = ck::AtomEnvironmentParams::default();
        let p = ck::AtomEnvironmentParams {
            use_hydrogens: use_hydrogens.unwrap_or(defaults.use_hydrogens),
            enforce_radius: enforce_radius.unwrap_or(defaults.enforce_radius),
        };
        self.inner
            .atom_environment_with_params(radius, ck::AtomId::new(root), &p)
            .map_err(PathError::from_rust)
            .map(|v| AtomEnvironment { inner: v })
    }

    /// Return source bonds and distances to atoms in a rooted environment. Uses the supplied configuration.
    #[pyo3(signature = (radius, root, params))]
    fn atom_environment_with_params(
        &self,
        py: Python<'_>,
        radius: usize,
        root: usize,
        params: &AtomEnvironmentParams,
    ) -> PyResult<AtomEnvironment> {
        self._atom_environment_params(py, radius, root, Some(params))
    }

    /// Return the bond index, or None for two valid atoms without a bond.
    #[pyo3(signature = (begin,end))]
    fn bond_between_atoms(&self, begin: usize, end: usize) -> PyResult<Option<usize>> {
        self.inner
            .bond_between_atoms(ck::AtomId::new(begin), ck::AtomId::new(end))
            .map_err(PathError::from_rust)
            .map(|v| v.map(|i| i.index()))
    }

    /// Return a shortest path as atom indices; disconnected atoms give an empty path.
    #[pyo3(signature = (begin,end))]
    fn shortest_path(&self, begin: usize, end: usize) -> PyResult<Vec<usize>> {
        self.inner
            .shortest_path(ck::AtomId::new(begin), ck::AtomId::new(end))
            .map_err(PathError::from_rust)
            .map(|v| v.into_iter().map(|i| i.index()).collect())
    }

    /// Return adjacency in atom order, optionally weighted by each bond's valence contribution to the row atom.
    #[pyo3(signature = (use_bond_order=false))]
    fn adjacency_matrix(&self, use_bond_order: bool) -> PyResult<DenseMatrix> {
        self.inner
            .adjacency_matrix(use_bond_order)
            .map_err(MatrixError::from_rust)
            .map(|v| DenseMatrix { inner: v })
    }
}
