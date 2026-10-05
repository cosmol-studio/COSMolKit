//! Canonical checked detached construction. Chemistry remains a Molecule operation.
use ::cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct AtomSpec {
    pub(crate) inner: ck::AtomSpec,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl AtomSpec {
    #[new]
    fn new(element: &crate::canonical_element_metadata::Element) -> Self {
        Self {
            inner: ck::AtomSpec::new(element.inner),
        }
    }
    fn with_formal_charge(&self, value: i8) -> Self {
        Self {
            inner: self.inner.clone().with_formal_charge(value),
        }
    }
    fn with_explicit_hydrogens(&self, value: u8) -> Self {
        Self {
            inner: self.inner.clone().with_explicit_hydrogens(value),
        }
    }
    fn with_atom_map(&self, value: u32) -> Self {
        Self {
            inner: self.inner.clone().with_atom_map(value),
        }
    }
    fn with_isotope(&self, value: u16) -> Self {
        Self {
            inner: self.inner.clone().with_isotope(value),
        }
    }
    fn with_no_implicit(&self, value: bool) -> Self {
        Self {
            inner: self.inner.clone().with_no_implicit(value),
        }
    }
    fn __repr__(&self) -> String {
        format!("AtomSpec(element='{}')", self.inner.element().symbol())
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct BondSpec {
    pub(crate) inner: ck::BondSpec,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BondSpec {
    #[new]
    fn new(begin: usize, end: usize, order: i64) -> PyResult<Self> {
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
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct MoleculeBuilder {
    pub(crate) inner: ck::MoleculeBuilder,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl MoleculeBuilder {
    #[staticmethod]
    fn new() -> Self {
        Self {
            inner: ck::MoleculeBuilder::new(),
        }
    }
    fn add_atom(&mut self, spec: &AtomSpec) -> usize {
        self.inner.add_atom(spec.inner.clone()).index()
    }
    fn add_bond(&mut self, py: Python<'_>, spec: &BondSpec) -> PyResult<usize> {
        self.inner
            .add_bond(spec.inner.clone())
            .map(ck::BondId::index)
            .map_err(|e| crate::drawing_binding::operation_pyerr(py, e))
    }
    fn set_atom_formal_charge(&mut self, py: Python<'_>, atom: usize, charge: i8) -> PyResult<()> {
        self.inner
            .set_atom_formal_charge(ck::AtomId::new(atom), charge)
            .map_err(|e| crate::drawing_binding::operation_pyerr(py, e))
    }
    fn set_bond_order(&mut self, py: Python<'_>, bond: usize, order: i64) -> PyResult<()> {
        let order = ck::BondOrder::from_rdkit_code(order)
            .ok_or_else(|| PyValueError::new_err(format!("invalid BondOrder code: {order}")))?;
        self.inner
            .set_bond_order(ck::BondId::new(bond), order)
            .map_err(|e| crate::drawing_binding::operation_pyerr(py, e))
    }
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
    fn add_2d_conformer(&mut self, py: Python<'_>, coordinates: Vec<[f64; 2]>) -> PyResult<usize> {
        self.inner
            .add_2d_conformer(coordinates)
            .map_err(|e| crate::drawing_binding::operation_pyerr(py, e))
    }
    fn add_3d_conformer(&mut self, py: Python<'_>, coordinates: Vec<[f64; 3]>) -> PyResult<usize> {
        self.inner
            .add_3d_conformer(coordinates)
            .map_err(|e| crate::drawing_binding::operation_pyerr(py, e))
    }
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
