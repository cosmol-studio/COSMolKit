//! Detached construction state is owned entirely by the canonical Rust builder.
use crate::atom_bond_values::{Atom, Bond};
use crate::construction_values::AtomSpec;
use crate::detached_values::CoordinateBlock;
use crate::group_values::SubstanceGroup;
use crate::metadata_values::MoleculeProperties;
use crate::stereo_groups::StereoGroup;
use crate::topology_values::TopologyBlock;
use cosmolkit as ck;
use pyo3::prelude::*;

#[pyclass(module = "cosmolkit")]
/// Mutable construction state for assembling atoms, bonds, coordinates and metadata.
/// build() produces a validated Molecule without consuming this Python builder.
pub struct MoleculeBuilder {
    pub(crate) inner: ck::MoleculeBuilder,
}

#[pymethods]
impl MoleculeBuilder {
    #[new]
    fn construct() -> Self {
        Self::new()
    }
    #[staticmethod]
    /// Create an empty molecular builder.
    fn new() -> Self {
        Self {
            inner: ck::MoleculeBuilder::new(),
        }
    }
    #[staticmethod]
    /// Initialize a builder from detached topology, coordinate and property values.
    fn from_parts(
        topology: &TopologyBlock,
        coordinates: &CoordinateBlock,
        properties: &MoleculeProperties,
    ) -> Self {
        Self {
            inner: ck::MoleculeBuilder::from_parts(
                topology.inner().clone(),
                coordinates.inner.clone(),
                properties.inner.clone(),
            ),
        }
    }
    /// Append an atom specification and return its zero-based atom index.
    fn add_atom(&mut self, atom: &AtomSpec) -> usize {
        self.inner.add_atom(atom.inner.clone()).index()
    }
    /// Append a bond specification and return its zero-based bond index; invalid endpoints raise OperationError.
    fn add_bond(&mut self, bond: &crate::construction_values::BondSpec) -> PyResult<usize> {
        self.inner
            .add_bond(bond.inner.clone())
            .map(ck::BondId::index)
            .map_err(crate::operation_errors::OperationError::from_rust)
    }
    /// Update the specified atom's formal charge in this builder.
    fn set_atom_formal_charge(&mut self, atom_id: usize, formal_charge: i8) -> PyResult<()> {
        self.inner
            .set_atom_formal_charge(ck::AtomId::new(atom_id), formal_charge)
            .map_err(crate::operation_errors::OperationError::from_rust)
    }
    /// Update the specified bond's order in this builder.
    fn set_bond_order(
        &mut self,
        bond_id: usize,
        order: crate::chemical_values::BondOrderInput,
    ) -> PyResult<()> {
        self.inner
            .set_bond_order(ck::BondId::new(bond_id), order.into_rust()?)
            .map_err(crate::operation_errors::OperationError::from_rust)
    }
    /// Remove a bond between the two atom indices; return whether a bond was removed.
    fn remove_bond_between_atoms(&mut self, begin_atom: usize, end_atom: usize) -> PyResult<bool> {
        self.inner
            .remove_bond_between_atoms(ck::AtomId::new(begin_atom), ck::AtomId::new(end_atom))
            .map_err(crate::operation_errors::OperationError::from_rust)
    }
    /// Return the number of neighbors of the specified atom in the current builder topology.
    fn degree(&self, atom_id: usize) -> PyResult<usize> {
        self.inner
            .degree(ck::AtomId::new(atom_id))
            .map_err(crate::operation_errors::OperationError::from_rust)
    }
    /// Return bond indices incident to the specified atom.
    fn neighbor_bonds(&self, atom_id: usize) -> PyResult<Vec<usize>> {
        self.inner
            .neighbor_bonds(ck::AtomId::new(atom_id))
            .map(|values| values.into_iter().map(ck::BondId::index).collect())
            .map_err(crate::operation_errors::OperationError::from_rust)
    }
    /// Return the connecting bond index, or None when no bond connects the specified atoms.
    fn bond_between_atoms(&self, begin_atom: usize, end_atom: usize) -> PyResult<Option<usize>> {
        self.inner
            .bond_between_atoms(ck::AtomId::new(begin_atom), ck::AtomId::new(end_atom))
            .map(|value| value.map(ck::BondId::index))
            .map_err(crate::operation_errors::OperationError::from_rust)
    }
    /// Set the builder's two-dimensional coordinates from one (x, y) pair per atom.
    fn set_2d_coordinates(&mut self, coordinates: Vec<[f64; 2]>) -> PyResult<()> {
        self.inner
            .set_2d_coordinates(coordinates)
            .map_err(crate::operation_errors::OperationError::from_rust)
    }
    /// Append a two-dimensional conformer and return its identifier.
    fn add_2d_conformer(&mut self, coordinates: Vec<[f64; 2]>) -> PyResult<usize> {
        self.inner
            .add_2d_conformer(coordinates)
            .map_err(crate::operation_errors::OperationError::from_rust)
    }
    /// Append a three-dimensional conformer and return its identifier.
    fn add_3d_conformer(&mut self, coordinates: Vec<[f64; 3]>) -> PyResult<usize> {
        self.inner
            .add_3d_conformer(coordinates)
            .map_err(crate::operation_errors::OperationError::from_rust)
    }
    /// Append a substance group and return its index; native structural validation errors propagate.
    fn add_substance_group(&mut self, group: &SubstanceGroup) -> PyResult<usize> {
        self.inner
            .add_substance_group(group.inner.clone())
            .map(|id| id.index())
            .map_err(crate::operation_errors::OperationError::from_rust)
    }
    /// Append an enhanced stereochemistry group and return its index.
    fn add_stereo_group(&mut self, group: &StereoGroup) -> PyResult<usize> {
        self.inner
            .add_stereo_group(group.inner.clone())
            .map_err(crate::operation_errors::OperationError::from_rust)
    }
    #[getter]
    /// Independent list of detached atom values in builder order.
    fn atoms(&self) -> Vec<Atom> {
        self.inner
            .atoms()
            .iter()
            .cloned()
            .map(|inner| Atom { inner })
            .collect()
    }
    #[getter]
    /// Independent list of detached bond values in builder order.
    fn bonds(&self) -> Vec<Bond> {
        self.inner
            .bonds()
            .iter()
            .cloned()
            .map(|inner| Bond { inner })
            .collect()
    }
    #[getter]
    /// Independent list of substance-group values.
    fn substance_groups(&self) -> Vec<SubstanceGroup> {
        self.inner
            .substance_groups()
            .iter()
            .cloned()
            .map(|inner| SubstanceGroup { inner })
            .collect()
    }
    #[getter]
    /// Independent list of enhanced stereochemistry-group values.
    fn stereo_groups(&self) -> Vec<StereoGroup> {
        self.inner
            .stereo_groups()
            .iter()
            .cloned()
            .map(|inner| StereoGroup { inner })
            .collect()
    }
    #[getter]
    /// Detached coordinate-block value from the builder.
    fn coordinates(&self) -> CoordinateBlock {
        CoordinateBlock {
            inner: self.inner.coordinates().clone(),
        }
    }
    #[getter]
    /// Detached molecule-property value from the builder.
    fn properties(&self) -> MoleculeProperties {
        MoleculeProperties {
            inner: self.inner.properties().clone(),
        }
    }

    /// Validate and construct a Molecule while keeping this Python builder reusable.
    fn build(&self) -> PyResult<crate::native::Molecule> {
        self.inner
            .clone()
            .build()
            .map(|inner| crate::native::Molecule { inner })
            .map_err(crate::operation_errors::OperationError::from_rust)
    }

    /// Return an independent builder value with a new molecule name.
    fn with_name(&self, name: String) -> Self {
        Self {
            inner: self.inner.clone().with_name(name),
        }
    }

    /// Return an independent builder value with the specified string property.
    fn with_property(&self, key: String, value: String) -> PyResult<Self> {
        self.inner
            .clone()
            .with_property(key, value)
            .map(|inner| Self { inner })
            .map_err(crate::operation_errors::OperationError::from_rust)
    }

    /// Return an independent builder value with an SDF data field.
    fn with_sdf_data_field(&self, key: String, value: String) -> Self {
        Self {
            inner: self.inner.clone().with_sdf_data_field(key, value),
        }
    }

    /// Return an independent builder value with the supplied molecule properties.
    fn with_properties(&self, properties: &MoleculeProperties) -> Self {
        Self {
            inner: self.inner.clone().with_properties(properties.inner.clone()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::native::Molecule;
    #[test]
    fn consuming_builder_calls_preserve_receiver_and_forward_native_cow_values() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let source = ck::Molecule::from_smiles("CCO").unwrap();
            let native = source.to_builder();
            let actual = Py::new(
                py,
                MoleculeBuilder {
                    inner: native.clone(),
                },
            )?;
            let built = actual
                .bind(py)
                .call_method0("build")?
                .extract::<Py<Molecule>>()?;
            assert_eq!(built.borrow(py).inner, native.clone().build().unwrap());
            assert!(std::ptr::eq(
                built.borrow(py).inner.topology(),
                source.topology()
            ));
            for (method, args, expected) in [
                (
                    "with_name",
                    vec!["renamed"],
                    native.clone().with_name("renamed".into()),
                ),
                (
                    "with_property",
                    vec!["key", "value"],
                    native
                        .clone()
                        .with_property("key".into(), "value".into())
                        .unwrap(),
                ),
                (
                    "with_sdf_data_field",
                    vec!["tag", "text"],
                    native
                        .clone()
                        .with_sdf_data_field("tag".into(), "text".into()),
                ),
            ] {
                let args = pyo3::types::PyTuple::new(py, args)?;
                let result = actual
                    .bind(py)
                    .call_method1(method, args)?
                    .extract::<Py<MoleculeBuilder>>()?;
                assert_eq!(result.borrow(py).inner, expected);
                assert_eq!(actual.borrow(py).inner, native);
                let result = result
                    .bind(py)
                    .call_method0("build")?
                    .extract::<Py<Molecule>>()?;
                assert!(std::ptr::eq(
                    result.borrow(py).inner.topology(),
                    source.topology()
                ));
            }
            let props = Py::new(
                py,
                MoleculeProperties {
                    inner: ck::MoleculeProperties::default().with_name("replacement"),
                },
            )?;
            let result = actual
                .bind(py)
                .call_method1("with_properties", (props.clone_ref(py),))?
                .extract::<Py<MoleculeBuilder>>()?;
            assert_eq!(
                result.borrow(py).inner,
                native
                    .clone()
                    .with_properties(props.borrow(py).inner.clone())
            );
            let failure = actual
                .bind(py)
                .call_method1("with_property", ("", "x"))
                .unwrap_err();
            assert!(failure.is_instance(py, <crate::operation_errors::OperationError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
            assert_eq!(actual.borrow(py).inner, native);
            actual
                .bind(py)
                .call_method1("set_atom_formal_charge", (0, 1))?;
            assert_eq!(source.topology().atoms[0].formal_charge(), 0);
            assert_eq!(
                built.borrow(py).inner.topology().atoms[0].formal_charge(),
                0
            );
            Ok(())
        })
        .unwrap();
    }
    #[test]
    fn non_consuming_builder_edit_query_and_from_parts_calls_match_native_values_and_errors() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let source = ck::Molecule::from_smiles("CCO").unwrap();
            let before = source.to_binary().unwrap();
            let mut expected = source.to_builder();
            let actual = Py::new(
                py,
                MoleculeBuilder {
                    inner: source.to_builder(),
                },
            )?;
            let bond = ck::BondSpec::new(
                ck::AtomId::new(0),
                ck::AtomId::new(2),
                ck::BondOrder::Single,
            );
            let expected_id = expected.add_bond(bond.clone()).unwrap().index();
            let projected_bond = Py::new(py, crate::construction_values::BondSpec { inner: bond })?;
            assert_eq!(
                actual
                    .bind(py)
                    .call_method1("add_bond", (projected_bond,))?
                    .extract::<usize>()?,
                expected_id
            );
            assert_eq!(actual.borrow(py).inner, expected);
            expected
                .set_atom_formal_charge(ck::AtomId::new(0), -1)
                .unwrap();
            actual
                .bind(py)
                .call_method1("set_atom_formal_charge", (0, -1))?;
            assert_eq!(actual.borrow(py).inner, expected);
            expected
                .set_bond_order(ck::BondId::new(0), ck::BondOrder::Double)
                .unwrap();
            actual
                .bind(py)
                .call_method1("set_bond_order", (0, "double"))?;
            assert_eq!(actual.borrow(py).inner, expected);
            for atom in 0..source.num_atoms() {
                assert_eq!(
                    actual
                        .bind(py)
                        .call_method1("degree", (atom,))?
                        .extract::<usize>()?,
                    expected.degree(ck::AtomId::new(atom)).unwrap()
                );
                assert_eq!(
                    actual
                        .bind(py)
                        .call_method1("neighbor_bonds", (atom,))?
                        .extract::<Vec<usize>>()?,
                    expected
                        .neighbor_bonds(ck::AtomId::new(atom))
                        .unwrap()
                        .into_iter()
                        .map(ck::BondId::index)
                        .collect::<Vec<_>>()
                );
            }
            for pair in [(0, 2), (1, 2)] {
                assert_eq!(
                    actual
                        .bind(py)
                        .call_method1("bond_between_atoms", pair)?
                        .extract::<Option<usize>>()?,
                    expected
                        .bond_between_atoms(ck::AtomId::new(pair.0), ck::AtomId::new(pair.1))
                        .unwrap()
                        .map(ck::BondId::index)
                );
            }
            let removed = expected
                .remove_bond_between_atoms(ck::AtomId::new(0), ck::AtomId::new(2))
                .unwrap();
            assert_eq!(
                actual
                    .bind(py)
                    .call_method1("remove_bond_between_atoms", (0, 2))?
                    .extract::<bool>()?,
                removed
            );
            assert_eq!(actual.borrow(py).inner, expected);
            let xy = vec![[0.0, 1.0], [2.0, 3.0], [4.0, 5.0]];
            expected.set_2d_coordinates(xy.clone()).unwrap();
            actual
                .bind(py)
                .call_method1("set_2d_coordinates", (xy.clone(),))?;
            let next = expected.add_2d_conformer(xy.clone()).unwrap();
            assert_eq!(
                actual
                    .bind(py)
                    .call_method1("add_2d_conformer", (xy,))?
                    .extract::<usize>()?,
                next
            );
            let xyz = vec![[0.0, 1.0, 2.0], [3.0, 4.0, 5.0], [6.0, 7.0, 8.0]];
            let next = expected.add_3d_conformer(xyz.clone()).unwrap();
            assert_eq!(
                actual
                    .bind(py)
                    .call_method1("add_3d_conformer", (xyz,))?
                    .extract::<usize>()?,
                next
            );
            assert_eq!(actual.borrow(py).inner, expected);
            let group = ck::SubstanceGroup::new(
                ck::SubstanceGroupId::new(99),
                ck::SubstanceGroupKind::Data,
            )
            .with_atoms(vec![ck::AtomId::new(0)]);
            let projected = Py::new(
                py,
                SubstanceGroup {
                    inner: group.clone(),
                },
            )?;
            let group_id = expected.add_substance_group(group).unwrap().index();
            assert_eq!(
                actual
                    .bind(py)
                    .call_method1("add_substance_group", (projected.clone_ref(py),))?
                    .extract::<usize>()?,
                group_id
            );
            assert_eq!(projected.borrow(py).inner.id().index(), 99);
            let stereo =
                ck::StereoGroup::new(ck::StereoGroupKind::Or, vec![ck::AtomId::new(0)], vec![])
                    .unwrap()
                    .with_id(7);
            let projected = Py::new(
                py,
                StereoGroup {
                    inner: stereo.clone(),
                },
            )?;
            let group_id = expected.add_stereo_group(stereo).unwrap();
            assert_eq!(
                actual
                    .bind(py)
                    .call_method1("add_stereo_group", (projected.clone_ref(py),))?
                    .extract::<usize>()?,
                group_id
            );
            assert_eq!(projected.borrow(py).inner.id(), Some(7));
            assert_eq!(actual.borrow(py).inner, expected);
            let original = actual.borrow(py).inner.clone();
            let failure = actual
                .bind(py)
                .call_method1("set_2d_coordinates", (vec![[0.0, 0.0]],))
                .unwrap_err();
            assert!(failure.is_instance(py, <crate::operation_errors::OperationError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
            assert_eq!(actual.borrow(py).inner, original);
            let failure = actual
                .bind(py)
                .call_method1("degree", (usize::MAX,))
                .unwrap_err();
            let native_failure = expected.degree(ck::AtomId::new(usize::MAX)).unwrap_err();
            assert_eq!(
                failure.value(py).str()?.to_str()?,
                native_failure.to_string()
            );
            assert!(failure.cause(py).is_none());
            assert_eq!(actual.borrow(py).inner, original);
            let native = ck::Molecule::from_parts(
                ck::TopologyBlock::try_from_parts(
                    expected.atoms().to_vec(),
                    expected.bonds().to_vec(),
                    expected.substance_groups().to_vec(),
                    expected.stereo_groups().to_vec(),
                )
                .unwrap(),
                expected.coordinates().clone(),
                expected.properties().clone(),
            )
            .unwrap();
            let topology = Py::new(py, TopologyBlock::from_molecule(&native))?;
            let coordinates = Py::new(
                py,
                CoordinateBlock {
                    inner: expected.coordinates().clone(),
                },
            )?;
            let properties = Py::new(
                py,
                MoleculeProperties {
                    inner: expected.properties().clone(),
                },
            )?;
            let built = py
                .get_type::<Molecule>()
                .call_method1("from_parts", (topology, coordinates, properties))?
                .extract::<Py<Molecule>>()?;

            assert_eq!(built.borrow(py).inner, native);
            assert_eq!(actual.borrow(py).inner, expected);
            assert_eq!(source.to_binary().unwrap(), before);
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn real_builder_keeps_detached_inputs_rows_and_source_molecule() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let source = ck::Molecule::from_smiles("CCO").unwrap();
            let mut native = source.to_builder();
            native
                .add_3d_conformer(vec![[0., 0., 0.], [1., 0.2, 0.], [2., 0.4, 0.1]])
                .unwrap();
            let native = native
                .with_name("native source".to_owned())
                .with_sdf_data_field("tag".into(), "value".into());
            let coordinates = CoordinateBlock {
                inner: native.coordinates().clone(),
            };
            let properties = MoleculeProperties {
                inner: native.properties().clone(),
            };
            let projected = MoleculeBuilder::from_parts(
                &TopologyBlock::from_molecule(&source),
                &coordinates,
                &properties,
            );
            assert_eq!(projected.inner.atoms(), native.atoms());
            assert_eq!(projected.inner.bonds(), native.bonds());
            assert_eq!(projected.inner.coordinates(), native.coordinates());
            assert_eq!(projected.inner.properties(), native.properties());
            let projected = Py::new(py, projected)?;
            assert_eq!(projected.bind(py).getattr("atoms")?.len()?, 3);
            assert_eq!(projected.bind(py).getattr("bonds")?.len()?, 2);
            assert_eq!(
                projected
                    .bind(py)
                    .getattr("properties")?
                    .getattr("name")?
                    .extract::<String>()?,
                "native source"
            );
            let atom = Py::new(
                py,
                AtomSpec {
                    inner: ck::AtomSpec::new(ck::Element::N),
                },
            )?;
            assert_eq!(
                projected
                    .bind(py)
                    .call_method1("add_atom", (atom,))?
                    .extract::<usize>()?,
                3
            );
            assert_eq!(projected.borrow(py).inner.atoms().len(), 4);
            assert_eq!(source.num_atoms(), 3);
            assert_eq!(native.atoms().len(), 3);
            let source = Py::new(py, Molecule { inner: source })?;
            let detached = source.bind(py).call_method0("to_builder")?;
            assert_eq!(detached.getattr("atoms")?.len()?, 3);
            assert_eq!(
                py.get_type::<MoleculeBuilder>()
                    .call0()?
                    .getattr("atoms")?
                    .len()?,
                0
            );
            Ok(())
        })
        .unwrap();
    }
}
