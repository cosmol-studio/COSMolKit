//! Canonical checked detached construction. Chemistry remains a Molecule operation.
use ::cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

/// Atom construction specification: element, charge, hydrogen count, isotope and map information. with_* methods return updated specifications.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct AtomSpec {
    pub(crate) inner: ck::AtomSpec,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl AtomSpec {
    /// Construct a AtomSpec value from the supplied inputs.
    #[new]
    fn py_new(element: &crate::canonical_element_metadata::Element) -> Self {
        Self::new(element)
    }
    /// Construct a AtomSpec value from the supplied inputs.
    #[staticmethod]
    fn new(element: &crate::canonical_element_metadata::Element) -> Self {
        Self {
            inner: ck::AtomSpec::new(element.inner),
        }
    }
    /// Return a new AtomSpec with the requested formal charge.
    fn with_formal_charge(&self, value: i8) -> Self {
        Self {
            inner: self.inner.clone().with_formal_charge(value),
        }
    }
    /// Return a new AtomSpec with the requested atom-stored explicit hydrogen count.
    fn with_explicit_hydrogens(&self, value: u8) -> Self {
        Self {
            inner: self.inner.clone().with_explicit_hydrogens(value),
        }
    }
    /// Return a new AtomSpec with the requested reaction atom-map number.
    fn with_atom_map(&self, value: u32) -> Self {
        Self {
            inner: self.inner.clone().with_atom_map(value),
        }
    }
    /// Return a new AtomSpec with the requested isotope mass number.
    fn with_isotope(&self, value: u16) -> Self {
        Self {
            inner: self.inner.clone().with_isotope(value),
        }
    }
    /// Return a new AtomSpec controlling whether implicit hydrogens are allowed.
    fn with_no_implicit(&self, value: bool) -> Self {
        Self {
            inner: self.inner.clone().with_no_implicit(value),
        }
    }
    fn __repr__(&self) -> String {
        format!("AtomSpec(element='{}')", self.inner.element().symbol())
    }
}
/// Bond construction specification containing begin/end atom indices and bond order.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct BondSpec {
    pub(crate) inner: ck::BondSpec,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl BondSpec {
    /// Construct a BondSpec value from the supplied inputs.
    #[new]
    fn py_new(
        begin: usize,
        end: usize,
        #[gen_stub(override_type(type_repr = "BondOrder | builtins.str | builtins.int"))]
        order: &Bound<'_, PyAny>,
    ) -> PyResult<Self> {
        Self::new(begin, end, order)
    }
    /// Construct a BondSpec value from the supplied inputs.
    #[staticmethod]
    fn new(
        begin: usize,
        end: usize,
        #[gen_stub(override_type(type_repr = "BondOrder | builtins.str | builtins.int"))]
        order: &Bound<'_, PyAny>,
    ) -> PyResult<Self> {
        let order = crate::canonical_atom_bond::enum_code(order, "BondOrder")?;
        let order = ck::BondOrder::from_rdkit_code(order)
            .ok_or_else(|| PyValueError::new_err(format!("invalid BondOrder code: {order}")))?;
        Ok(Self {
            inner: ck::BondSpec::new(ck::AtomId::new(begin), ck::AtomId::new(end), order),
        })
    }
    fn __repr__(&self) -> String {
        format!(
            "BondSpec(begin={}, end={}, order='{}')",
            self.inner.begin(),
            self.inner.end(),
            self.inner.order().rdkit_name()
        )
    }
}
/// Editable molecular construction state.
///
/// Add atoms/bonds and set coordinates or properties, then call build() to validate
/// and obtain a Molecule. Editing a builder created from a molecule does not change
/// the source molecule.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct MoleculeBuilder {
    pub(crate) inner: ck::MoleculeBuilder,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl MoleculeBuilder {
    /// Construct a value from explicit detached parts; required structural consistency is checked at the public boundary.
    #[staticmethod]
    fn from_parts(
        topology: &crate::canonical_detached_blocks::TopologyBlock,
        coordinates: &crate::canonical_detached_blocks::CoordinateBlock,
        properties: &crate::canonical_property_values::MoleculeProperties,
    ) -> Self {
        Self {
            inner: ck::MoleculeBuilder::from_parts(
                topology.inner.clone(),
                coordinates.inner.clone(),
                properties.inner.clone(),
            ),
        }
    }
    /// Return atom rows in stored graph/hierarchy order.
    fn atoms(&self, py: Python<'_>) -> PyResult<Vec<crate::canonical_atom_bond::Atom>> {
        self.inner
            .atoms()
            .iter()
            .map(|atom| {
                self.inner
                    .degree(atom.id())
                    .map(|degree| {
                        crate::canonical_atom_bond::Atom::from_detached(atom.clone(), degree)
                    })
                    .map_err(|error| crate::drawing_binding::operation_pyerr(py, error))
            })
            .collect()
    }
    /// Return the builder coordinate block, keeping 2D and 3D storage separate.
    fn coordinates(&self) -> crate::canonical_detached_blocks::CoordinateBlock {
        crate::canonical_detached_blocks::CoordinateBlock {
            inner: self.inner.coordinates().clone(),
        }
    }
    /// Return bond rows in graph order.
    fn bonds(&self) -> Vec<crate::canonical_atom_bond::Bond> {
        self.inner
            .bonds()
            .iter()
            .cloned()
            .map(|inner| crate::canonical_atom_bond::Bond { inner })
            .collect()
    }
    /// Substance-group annotations referencing graph atoms and bonds.
    fn substance_groups(&self) -> Vec<crate::canonical_group_values::SubstanceGroup> {
        self.inner
            .substance_groups()
            .iter()
            .cloned()
            .map(|inner| crate::canonical_group_values::SubstanceGroup { inner })
            .collect()
    }
    /// Enhanced stereochemistry groups referencing graph atoms.
    fn stereo_groups(&self) -> Vec<crate::canonical_group_values::StereoGroup> {
        self.inner
            .stereo_groups()
            .iter()
            .cloned()
            .map(|inner| crate::canonical_group_values::StereoGroup { inner })
            .collect()
    }
    /// Append a substance-group annotation using local atom/bond identifiers.
    fn add_substance_group(
        &mut self,
        py: Python<'_>,
        group: &crate::canonical_group_values::SubstanceGroup,
    ) -> PyResult<crate::canonical_group_values::SubstanceGroupId> {
        self.inner
            .add_substance_group(group.inner.clone())
            .map(|inner| crate::canonical_group_values::SubstanceGroupId { inner })
            .map_err(|error| crate::drawing_binding::operation_pyerr(py, error))
    }
    /// Append an enhanced-stereochemistry group referencing local atoms.
    fn add_stereo_group(
        &mut self,
        py: Python<'_>,
        group: &crate::canonical_group_values::StereoGroup,
    ) -> PyResult<usize> {
        self.inner
            .add_stereo_group(group.inner.clone())
            .map_err(|error| crate::drawing_binding::operation_pyerr(py, error))
    }
    /// Construct a MoleculeBuilder value from the supplied inputs.
    #[staticmethod]
    fn new() -> Self {
        Self {
            inner: ck::MoleculeBuilder::new(),
        }
    }
    /// Append an atom using AtomSpec and return its zero-based atom index.
    fn add_atom(&mut self, spec: &AtomSpec) -> usize {
        self.inner.add_atom(spec.inner.clone()).index()
    }
    /// Append a bond using BondSpec; validate its endpoint atom indices.
    fn add_bond(&mut self, py: Python<'_>, spec: &BondSpec) -> PyResult<usize> {
        self.inner
            .add_bond(spec.inner.clone())
            .map(ck::BondId::index)
            .map_err(|e| crate::drawing_binding::operation_pyerr(py, e))
    }
    /// Update the formal charge of the specified builder atom.
    fn set_atom_formal_charge(&mut self, py: Python<'_>, atom: usize, charge: i8) -> PyResult<()> {
        self.inner
            .set_atom_formal_charge(ck::AtomId::new(atom), charge)
            .map_err(|e| crate::drawing_binding::operation_pyerr(py, e))
    }

    /// Remove the bond joining the two specified atom indices, if one is present.
    fn remove_bond_between_atoms(
        &mut self,
        py: Python<'_>,
        begin_atom: usize,
        end_atom: usize,
    ) -> PyResult<bool> {
        self.inner
            .remove_bond_between_atoms(ck::AtomId::new(begin_atom), ck::AtomId::new(end_atom))
            .map_err(|e| crate::drawing_binding::operation_pyerr(py, e))
    }

    /// Number of graph neighbors, including explicit hydrogen vertices.
    fn degree(&self, py: Python<'_>, atom_id: usize) -> PyResult<usize> {
        self.inner
            .degree(ck::AtomId::new(atom_id))
            .map_err(|e| crate::drawing_binding::operation_pyerr(py, e))
    }

    /// Return the bond indices incident to the specified builder atom.
    fn neighbor_bonds(&self, py: Python<'_>, atom_id: usize) -> PyResult<Vec<usize>> {
        self.inner
            .neighbor_bonds(ck::AtomId::new(atom_id))
            .map(|ids| ids.into_iter().map(ck::BondId::index).collect())
            .map_err(|e| crate::drawing_binding::operation_pyerr(py, e))
    }

    /// Return the bond joining the two specified builder atoms, when present.
    fn bond_between_atoms(
        &self,
        py: Python<'_>,
        begin_atom: usize,
        end_atom: usize,
    ) -> PyResult<Option<usize>> {
        self.inner
            .bond_between_atoms(ck::AtomId::new(begin_atom), ck::AtomId::new(end_atom))
            .map(|id| id.map(ck::BondId::index))
            .map_err(|e| crate::drawing_binding::operation_pyerr(py, e))
    }

    /// Return a detached snapshot of the stored properties.
    fn properties(&self) -> crate::canonical_property_values::MoleculeProperties {
        crate::canonical_property_values::MoleculeProperties {
            inner: self.inner.properties().clone(),
        }
    }

    /// Return a builder with the requested molecule name.
    fn with_name(&self, name: String) -> Self {
        Self {
            inner: self.inner.clone().with_name(name),
        }
    }

    /// Return a builder with the supplied molecule property key/value.
    fn with_property(&self, py: Python<'_>, key: String, value: String) -> PyResult<Self> {
        self.inner
            .clone()
            .with_property(key, value)
            .map(|inner| Self { inner })
            .map_err(|e| crate::drawing_binding::operation_pyerr(py, e))
    }

    /// Return a builder with an appended ordered SDF data field.
    fn with_sdf_data_field(&self, key: String, value: String) -> Self {
        Self {
            inner: self.inner.clone().with_sdf_data_field(key, value),
        }
    }

    /// Return a builder with the supplied detached molecule properties.
    fn with_properties(
        &self,
        properties: &crate::canonical_property_values::MoleculeProperties,
    ) -> Self {
        Self {
            inner: self.inner.clone().with_properties(properties.inner.clone()),
        }
    }
    /// Update the bond order of the specified builder bond.
    fn set_bond_order(
        &mut self,
        py: Python<'_>,
        bond: usize,
        #[gen_stub(override_type(type_repr = "BondOrder | builtins.str | builtins.int"))]
        order: &Bound<'_, PyAny>,
    ) -> PyResult<()> {
        let order = crate::canonical_atom_bond::enum_code(order, "BondOrder")?;
        let order = ck::BondOrder::from_rdkit_code(order)
            .ok_or_else(|| PyValueError::new_err(format!("invalid BondOrder code: {order}")))?;
        self.inner
            .set_bond_order(ck::BondId::new(bond), order)
            .map_err(|e| crate::drawing_binding::operation_pyerr(py, e))
    }
    /// Set atom-ordered 2D positions on the builder; the row count must equal the atom count.
    fn set_2d_coordinates(
        &mut self,
        py: Python<'_>,
        coordinates: &Bound<'_, PyAny>,
    ) -> PyResult<()> {
        let rows = crate::canonical_coordinate_input::matrix(
            py,
            coordinates,
            self.inner.atoms().len(),
            true,
        )?;
        // This already registered builder accepts exactly XY; reject extra
        // columns at numeric ingress rather than creating a binding z policy.
        let rows = rows
            .into_iter()
            .enumerate()
            .map(|(row, values)| {
                if values.len() != 2 {
                    return Err(crate::canonical_coordinate_input::error_pyerr(
                        py,
                        &ck::CoordinateInputError::Shape {
                            dimension: "2D",
                            row,
                            columns: values.len(),
                            expected: "2",
                        },
                    ));
                }
                Ok([values[0], values[1]])
            })
            .collect::<PyResult<Vec<_>>>()?;
        self.inner
            .set_2d_coordinates(rows)
            .map_err(|e| crate::drawing_binding::operation_pyerr(py, e))
    }
    /// Adds one checked 2D conformer and returns its dimension-local id.
    fn add_2d_conformer(&mut self, py: Python<'_>, coordinates: Vec<[f64; 2]>) -> PyResult<usize> {
        self.inner
            .add_2d_conformer(coordinates)
            .map_err(|e| crate::drawing_binding::operation_pyerr(py, e))
    }
    /// Return a new result that will append a conformer using caller-supplied atom-ordered 3D coordinates. The source molecule is unchanged.
    fn add_3d_conformer(&mut self, py: Python<'_>, coordinates: Vec<[f64; 3]>) -> PyResult<usize> {
        self.inner
            .add_3d_conformer(coordinates)
            .map_err(|e| crate::drawing_binding::operation_pyerr(py, e))
    }
    /// Construct and validate an owned Molecule from the builder state; invalid graph or coordinate parts raise OperationError.
    fn build(&self, py: Python<'_>) -> PyResult<crate::drawing_binding::Molecule> {
        // Python retains its editing value. Pass an owned detached snapshot
        // to the exact Rust consuming signature; no sanitation/default flag.
        self.inner
            .clone()
            .build()
            .map(|inner| crate::drawing_binding::Molecule { inner })
            .map_err(|e| crate::drawing_binding::operation_pyerr(py, e))
    }
    fn __repr__(&self) -> String {
        format!(
            "MoleculeBuilder(num_atoms={}, num_bonds={})",
            self.inner.atoms().len(),
            self.inner.bonds().len()
        )
    }
}
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<AtomSpec>()?;
    module.add_class::<BondSpec>()?;
    module.add_class::<MoleculeBuilder>()?;
    Ok(())
}
