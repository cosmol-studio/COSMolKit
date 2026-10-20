use crate::chemical_values::BondOrderInput;
use crate::elements::Element;
use cosmolkit as ck;
use pyo3::prelude::*;

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen)]
    pub struct AtomSpec {
        pub(crate) inner: ck::AtomSpec,
    }
}

#[pymethods]
impl AtomSpec {
    #[new]
    fn py_new(element: &Element) -> Self {
        Self::new(element)
    }
    #[staticmethod]
    fn new(element: &Element) -> Self {
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
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen)]
    pub struct BondSpec {
        pub(crate) inner: ck::BondSpec,
    }
}

#[pymethods]
impl BondSpec {
    #[new]
    fn py_new(begin: usize, end: usize, order: BondOrderInput) -> PyResult<Self> {
        Self::new(begin, end, order)
    }
    #[staticmethod]
    fn new(begin: usize, end: usize, order: BondOrderInput) -> PyResult<Self> {
        Ok(Self {
            inner: ck::BondSpec::new(
                ck::AtomId::new(begin),
                ck::AtomId::new(end),
                order.into_rust()?,
            ),
        })
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    pub struct ReactionCoordinateSelection {
        pub(crate) inner: ck::ReactionCoordinateSelection,
    }
}

#[pymethods]
impl ReactionCoordinateSelection {
    #[staticmethod]
    fn auto() -> Self {
        Self {
            inner: ck::ReactionCoordinateSelection::auto(),
        }
    }
    #[staticmethod]
    fn two_d(id: usize) -> Self {
        Self {
            inner: ck::ReactionCoordinateSelection::two_d(id),
        }
    }
    #[staticmethod]
    fn three_d(id: usize) -> Self {
        Self {
            inner: ck::ReactionCoordinateSelection::three_d(id),
        }
    }
    #[getter]
    fn id(&self) -> Option<usize> {
        self.inner.id()
    }
    #[getter]
    fn is_auto(&self) -> bool {
        self.inner.is_auto()
    }
    #[getter]
    fn is_2d(&self) -> bool {
        self.inner.is_2d()
    }
    #[getter]
    fn is_3d(&self) -> bool {
        self.inner.is_3d()
    }
}
