use crate::chemical_values::BondOrderInput;
use crate::elements::Element;
use cosmolkit as ck;
use pyo3::prelude::*;

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen)]
    /// Detached inputs for constructing an atom; with_* methods return independent specifications.
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
    /// Create an atom specification for an element with native default atom attributes.
    fn new(element: &Element) -> Self {
        Self {
            inner: ck::AtomSpec::new(element.inner),
        }
    }
    /// Return an independent atom specification with the requested formal charge.
    fn with_formal_charge(&self, value: i8) -> Self {
        Self {
            inner: self.inner.clone().with_formal_charge(value),
        }
    }
    /// Return an independent atom specification with the explicit hydrogen count replaced.
    fn with_explicit_hydrogens(&self, value: u8) -> Self {
        Self {
            inner: self.inner.clone().with_explicit_hydrogens(value),
        }
    }
    /// Return an independent atom specification with the atom-map number replaced.
    fn with_atom_map(&self, value: u32) -> Self {
        Self {
            inner: self.inner.clone().with_atom_map(value),
        }
    }
    /// Return an independent atom specification with the isotope mass number replaced.
    fn with_isotope(&self, value: u16) -> Self {
        Self {
            inner: self.inner.clone().with_isotope(value),
        }
    }
    /// Return an independent atom specification with implicit-hydrogen inference enabled or disabled.
    fn with_no_implicit(&self, value: bool) -> Self {
        Self {
            inner: self.inner.clone().with_no_implicit(value),
        }
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen)]
    /// Detached bond-construction inputs: zero-based endpoint atom identifiers and bond order.
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
    /// Create a bond specification with zero-based endpoint indices and the requested bond order.
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
    /// Explicit or automatic selection of the coordinates used by a reaction.
    pub struct ReactionCoordinateSelection {
        pub(crate) inner: ck::ReactionCoordinateSelection,
    }
}

#[pymethods]
impl ReactionCoordinateSelection {
    #[staticmethod]
    /// Use the canonical automatic reaction coordinate selection.
    fn auto() -> Self {
        Self {
            inner: ck::ReactionCoordinateSelection::auto(),
        }
    }
    #[staticmethod]
    /// Select the stored 2D conformer with this identifier.
    fn two_d(id: usize) -> Self {
        Self {
            inner: ck::ReactionCoordinateSelection::two_d(id),
        }
    }
    #[staticmethod]
    /// Select the stored 3D conformer with this identifier.
    fn three_d(id: usize) -> Self {
        Self {
            inner: ck::ReactionCoordinateSelection::three_d(id),
        }
    }
    #[getter]
    /// Selected conformer identifier, or None for automatic selection.
    fn id(&self) -> Option<usize> {
        self.inner.id()
    }
    #[getter]
    /// Whether automatic coordinate selection is requested.
    fn is_auto(&self) -> bool {
        self.inner.is_auto()
    }
    #[getter]
    /// Whether an explicit 2D conformer is selected.
    fn is_2d(&self) -> bool {
        self.inner.is_2d()
    }
    #[getter]
    /// Whether an explicit 3D conformer is selected.
    fn is_3d(&self) -> bool {
        self.inner.is_3d()
    }
}
