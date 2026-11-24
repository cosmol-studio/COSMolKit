use crate::{
    arguments::{Argument, configuration_arguments},
    construction_values::{AtomSpec, BondSpec},
    native::Molecule,
    operation_errors::OperationError,
    params::SubgraphExtractionParams,
};
use cosmolkit as ck;
use pyo3::prelude::*;
crate::configuration::type_projection! {value(eq)
    /// Source atoms, source bonds, or both, selected for extraction or retention.
    pub struct SubgraphSelection {pub(crate) inner: ck::SubgraphSelection}
}
crate::configuration::type_projection! {value(eq)
    /// Source atom and bond index of each local row at the time of selection.
    pub struct SubgraphMapping {pub(crate) inner: ck::SubgraphMapping}
}
crate::configuration::type_projection! {value(eq)
    /// Extracted molecule and its independent local-to-source correspondence.
    pub struct MappedSubgraph {pub(crate) inner: ck::MappedSubgraph}
}
#[pymethods]
impl SubgraphSelection {
    /// Keep these atoms and all bonds between them.
    #[staticmethod]
    fn atoms(indices: Vec<usize>) -> Self {
        Self {
            inner: ck::SubgraphSelection::Atoms(indices.into_iter().map(ck::AtomId::new).collect()),
        }
    }
    /// Keep exactly these bonds and their endpoints.
    #[staticmethod]
    fn bonds(indices: Vec<usize>) -> Self {
        Self {
            inner: ck::SubgraphSelection::Bonds(indices.into_iter().map(ck::BondId::new).collect()),
        }
    }
    /// Keep exactly these atoms and bonds. Both endpoints of every bond must be selected.
    #[staticmethod]
    fn atoms_and_bonds(atoms: Vec<usize>, bonds: Vec<usize>) -> Self {
        Self {
            inner: ck::SubgraphSelection::AtomsAndBonds(
                atoms.into_iter().map(ck::AtomId::new).collect(),
                bonds.into_iter().map(ck::BondId::new).collect(),
            ),
        }
    }
}
#[pymethods]
impl SubgraphMapping {
    /// Source atom index of each local atom at call entry.
    #[getter]
    fn source_atom_indices(&self) -> Vec<usize> {
        self.inner
            .source_atom_indices()
            .iter()
            .map(|i| i.index())
            .collect()
    }
    /// Source bond index of each local bond at call entry.
    #[getter]
    fn source_bond_indices(&self) -> Vec<usize> {
        self.inner
            .source_bond_indices()
            .iter()
            .map(|i| i.index())
            .collect()
    }
}
#[pymethods]
impl MappedSubgraph {
    /// Each read yields an independently editable molecule.
    #[getter]
    fn molecule(&self) -> Molecule {
        Molecule {
            inner: self.inner.molecule().clone(),
        }
    }
    /// Correspondence to the source at call entry, independent of later edits.
    #[getter]
    fn mapping(&self) -> SubgraphMapping {
        SubgraphMapping {
            inner: self.inner.mapping().clone(),
        }
    }
}
#[pymethods]
impl Molecule {
    /// Apply a complete new-to-old permutation to atom and coordinate rows. Bonds retain their row order. Returns a new molecule; this molecule is unchanged.
    #[pyo3(signature=(order))]
    fn with_atoms_renumbered(&self, order: Vec<usize>) -> PyResult<Molecule> {
        self.inner
            .with_atoms_renumbered(&order.into_iter().map(ck::AtomId::new).collect::<Vec<_>>())
            .map(|inner| Molecule { inner })
            .map_err(OperationError::from_rust)
    }

    /// Apply a complete new-to-old permutation to atom and coordinate rows. Bonds retain their row order. Changes this molecule.
    #[pyo3(signature=(order))]
    fn renumber_atoms_(&mut self, order: Vec<usize>) -> PyResult<()> {
        self.inner
            .renumber_atoms_(&order.into_iter().map(ck::AtomId::new).collect::<Vec<_>>())
            .map_err(OperationError::from_rust)
    }

    /// Remove an atom and its incident bonds; surviving rows stay in source order. Returns a new molecule; this molecule is unchanged.
    #[pyo3(signature=(atom))]
    fn without_atom(&self, atom: usize) -> PyResult<Molecule> {
        self.inner
            .without_atom(ck::AtomId::new(atom))
            .map(|inner| Molecule { inner })
            .map_err(OperationError::from_rust)
    }

    /// Remove an atom and its incident bonds; surviving rows stay in source order. Changes this molecule.
    #[pyo3(signature=(atom))]
    fn remove_atom_(&mut self, atom: usize) -> PyResult<()> {
        self.inner
            .remove_atom_(ck::AtomId::new(atom))
            .map_err(OperationError::from_rust)
    }

    /// Replace one atom slot without moving its coordinates. Set preserve_props to use the original property dictionary. Returns a new molecule; this molecule is unchanged.
    #[pyo3(signature=(atom,replacement,preserve_props=false))]
    fn with_atom_replaced(
        &self,
        atom: usize,
        replacement: &AtomSpec,
        preserve_props: bool,
    ) -> PyResult<Molecule> {
        self.inner
            .with_atom_replaced(ck::AtomId::new(atom), &replacement.inner, preserve_props)
            .map(|inner| Molecule { inner })
            .map_err(OperationError::from_rust)
    }

    /// Replace one atom slot without moving its coordinates. Set preserve_props to use the original property dictionary. Changes this molecule.
    #[pyo3(signature=(atom,replacement,preserve_props=false))]
    fn replace_atom_(
        &mut self,
        atom: usize,
        replacement: &AtomSpec,
        preserve_props: bool,
    ) -> PyResult<()> {
        self.inner
            .replace_atom_(ck::AtomId::new(atom), &replacement.inner, preserve_props)
            .map_err(OperationError::from_rust)
    }

    /// Replace a bond, keeping its endpoints; supplied stereo neighbors belong to the target graph. Returns a new molecule; this molecule is unchanged.
    #[pyo3(signature=(bond,replacement,preserve_props=false,keep_sgroups=true))]
    fn with_bond_replaced(
        &self,
        bond: usize,
        replacement: &BondSpec,
        preserve_props: bool,
        keep_sgroups: bool,
    ) -> PyResult<Molecule> {
        self.inner
            .with_bond_replaced(
                ck::BondId::new(bond),
                &replacement.inner,
                preserve_props,
                keep_sgroups,
            )
            .map(|inner| Molecule { inner })
            .map_err(OperationError::from_rust)
    }

    /// Replace a bond, keeping its endpoints; supplied stereo neighbors belong to the target graph. Changes this molecule.
    #[pyo3(signature=(bond,replacement,preserve_props=false,keep_sgroups=true))]
    fn replace_bond_(
        &mut self,
        bond: usize,
        replacement: &BondSpec,
        preserve_props: bool,
        keep_sgroups: bool,
    ) -> PyResult<()> {
        self.inner
            .replace_bond_(
                ck::BondId::new(bond),
                &replacement.inner,
                preserve_props,
                keep_sgroups,
            )
            .map_err(OperationError::from_rust)
    }

    /// Append the other graph without an inter-molecular bond. Match coordinate sets by dimension and ID. Returns a new molecule; this molecule is unchanged.
    #[pyo3(signature=(right))]
    fn with_molecule_appended(&self, right: &Molecule) -> PyResult<Molecule> {
        self.inner
            .with_molecule_appended(&right.inner)
            .map(|inner| Molecule { inner })
            .map_err(OperationError::from_rust)
    }

    /// Append the other graph without an inter-molecular bond. Match coordinate sets by dimension and ID. Changes this molecule.
    fn insert_molecule_(slf: &Bound<'_, Self>, right: &Bound<'_, Self>) -> PyResult<()> {
        let snapshot = right.try_borrow()?.inner.clone();
        slf.try_borrow_mut()?
            .inner
            .insert_molecule_(&snapshot)
            .map_err(OperationError::from_rust)
    }

    /// Return a new molecule and local-to-source mappings for the selected atoms and bonds. By default all conformers follow the selected rows and the subset is not sanitized.
    #[pyo3(signature=(selection,params=Argument::Missing,*,sanitize=Argument::Missing,copy_coordinates=Argument::Missing,clear_computed_props=Argument::Missing))]
    fn extract_subgraph(
        &self,
        py: Python<'_>,
        selection: &SubgraphSelection,
        params: Argument<Option<PyRef<'_, SubgraphExtractionParams>>>,
        sanitize: Argument<Option<bool>>,
        copy_coordinates: Argument<Option<bool>>,
        clear_computed_props: Argument<Option<bool>>,
    ) -> PyResult<MappedSubgraph> {
        let keywords = sanitize.is_supplied()
            || copy_coordinates.is_supplied()
            || clear_computed_props.is_supplied();
        configuration_arguments!(params; sanitize,copy_coordinates,clear_computed_props);
        if keywords {
            self._extract_subgraph_keywords(
                selection,
                sanitize,
                copy_coordinates,
                clear_computed_props,
            )
        } else {
            self._extract_subgraph_params(py, selection, params.as_deref())
        }
    }
    #[doc = r#"{"role":"overload","target":"Molecule.extract_subgraph"}"#]
    #[pyo3(signature=(selection,params=None))]
    fn _extract_subgraph_params(
        &self,
        py: Python<'_>,
        selection: &SubgraphSelection,
        params: Option<&SubgraphExtractionParams>,
    ) -> PyResult<MappedSubgraph> {
        let p = params
            .map(|p| p.to_rust(py))
            .transpose()?
            .unwrap_or_default();
        self.inner
            .extract_subgraph_with_params(&selection.inner, &p)
            .map(|inner| MappedSubgraph { inner })
            .map_err(OperationError::from_rust)
    }
    #[doc = r#"{"role":"overload","target":"Molecule.extract_subgraph"}"#]
    #[pyo3(signature=(selection,*,sanitize=None,copy_coordinates=None,clear_computed_props=None))]
    fn _extract_subgraph_keywords(
        &self,
        selection: &SubgraphSelection,
        sanitize: Option<bool>,
        copy_coordinates: Option<bool>,
        clear_computed_props: Option<bool>,
    ) -> PyResult<MappedSubgraph> {
        let defaults = ck::SubgraphExtractionParams::default();
        let p = ck::SubgraphExtractionParams {
            sanitize: sanitize.unwrap_or(defaults.sanitize),
            copy_coordinates: copy_coordinates.unwrap_or(defaults.copy_coordinates),
            clear_computed_props: clear_computed_props.unwrap_or(defaults.clear_computed_props),
        };
        self.inner
            .extract_subgraph_with_params(&selection.inner, &p)
            .map(|inner| MappedSubgraph { inner })
            .map_err(OperationError::from_rust)
    }
    /// Return a new molecule and local-to-source mappings for the selected atoms and bonds.
    fn extract_subgraph_with_params(
        &self,
        py: Python<'_>,
        selection: &SubgraphSelection,
        params: &SubgraphExtractionParams,
    ) -> PyResult<MappedSubgraph> {
        self._extract_subgraph_params(py, selection, Some(params))
    }

    /// Retain the selected graph in this molecule; return its local-to-source mapping. By default all conformers follow the selected rows and the subset is not sanitized.
    #[pyo3(signature=(selection,params=Argument::Missing,*,sanitize=Argument::Missing,copy_coordinates=Argument::Missing,clear_computed_props=Argument::Missing))]
    fn retain_subgraph_(
        &mut self,
        py: Python<'_>,
        selection: &SubgraphSelection,
        params: Argument<Option<PyRef<'_, SubgraphExtractionParams>>>,
        sanitize: Argument<Option<bool>>,
        copy_coordinates: Argument<Option<bool>>,
        clear_computed_props: Argument<Option<bool>>,
    ) -> PyResult<SubgraphMapping> {
        let keywords = sanitize.is_supplied()
            || copy_coordinates.is_supplied()
            || clear_computed_props.is_supplied();
        configuration_arguments!(params; sanitize,copy_coordinates,clear_computed_props);
        if keywords {
            self._retain_subgraph_keywords(
                selection,
                sanitize,
                copy_coordinates,
                clear_computed_props,
            )
        } else {
            self._retain_subgraph_params(py, selection, params.as_deref())
        }
    }
    #[doc = r#"{"role":"overload","target":"Molecule.retain_subgraph_"}"#]
    #[pyo3(signature=(selection,params=None))]
    fn _retain_subgraph_params(
        &mut self,
        py: Python<'_>,
        selection: &SubgraphSelection,
        params: Option<&SubgraphExtractionParams>,
    ) -> PyResult<SubgraphMapping> {
        let p = params
            .map(|p| p.to_rust(py))
            .transpose()?
            .unwrap_or_default();
        self.inner
            .retain_subgraph_with_params_(&selection.inner, &p)
            .map(|inner| SubgraphMapping { inner })
            .map_err(OperationError::from_rust)
    }
    #[doc = r#"{"role":"overload","target":"Molecule.retain_subgraph_"}"#]
    #[pyo3(signature=(selection,*,sanitize=None,copy_coordinates=None,clear_computed_props=None))]
    fn _retain_subgraph_keywords(
        &mut self,
        selection: &SubgraphSelection,
        sanitize: Option<bool>,
        copy_coordinates: Option<bool>,
        clear_computed_props: Option<bool>,
    ) -> PyResult<SubgraphMapping> {
        let defaults = ck::SubgraphExtractionParams::default();
        let p = ck::SubgraphExtractionParams {
            sanitize: sanitize.unwrap_or(defaults.sanitize),
            copy_coordinates: copy_coordinates.unwrap_or(defaults.copy_coordinates),
            clear_computed_props: clear_computed_props.unwrap_or(defaults.clear_computed_props),
        };
        self.inner
            .retain_subgraph_with_params_(&selection.inner, &p)
            .map(|inner| SubgraphMapping { inner })
            .map_err(OperationError::from_rust)
    }
    /// Retain the selected graph in this molecule; return its local-to-source mapping.
    fn retain_subgraph_with_params_(
        &mut self,
        py: Python<'_>,
        selection: &SubgraphSelection,
        params: &SubgraphExtractionParams,
    ) -> PyResult<SubgraphMapping> {
        self._retain_subgraph_params(py, selection, Some(params))
    }
}
