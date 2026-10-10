//! Immutable element values and result-field projections of the public facade.

use ::cosmolkit as ck;
use pyo3::prelude::*;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

/// Stable chemical-element identity.
///
/// The private atomic number is always in the inclusive range ``0..=118``.
/// Zero is the source-compatible dummy atom (``*``); ``1..=118`` are H through
/// Og.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, eq, hash, from_py_object)]
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) struct Element {
    pub(crate) inner: ck::Element,
}

impl Element {
    pub(crate) fn from_inner(inner: ck::Element) -> Self {
        Self { inner }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl Element {
    /// Construct an element from a checked atomic number.
    #[staticmethod]
    fn from_atomic_number(atomic_number: u8) -> Option<Self> {
        ck::Element::from_atomic_number(atomic_number).map(|inner| Self { inner })
    }

    /// Construct an element from a source-recognized symbol.
    #[staticmethod]
    fn from_symbol(symbol: &str) -> Option<Self> {
        ck::Element::from_symbol(symbol).map(|inner| Self { inner })
    }

    /// Atomic number (proton count); zero denotes a dummy atom.
    fn atomic_number(&self) -> u8 {
        self.inner.atomic_number()
    }

    /// Chemical element symbol represented by this value.
    fn symbol(&self) -> &'static str {
        self.inner.symbol()
    }

    fn __str__(&self) -> String {
        self.inner.to_string()
    }

    fn __repr__(&self) -> String {
        format!(
            "Element(symbol='{}', atomic_number={})",
            self.inner.symbol(),
            self.inner.atomic_number()
        )
    }
}

/// Periodic-table metadata for one element, including its symbol, atomic number, atomic weights and allowed valences.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct ElementInfo {
    inner: ck::ElementInfo,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl ElementInfo {
    // Read-callable transport of the existing eight public result fields.
    /// Chemical element as an Element value.
    fn element(&self) -> Element {
        Element {
            inner: self.inner.element,
        }
    }

    /// Chemical element symbol represented by this value.
    fn symbol(&self) -> &'static str {
        self.inner.symbol
    }

    /// Atomic number (proton count); zero denotes a dummy atom.
    fn atomic_number(&self) -> u8 {
        self.inner.atomic_number
    }

    /// Periodic-table period of the element.
    fn period(&self) -> u8 {
        self.inner.period
    }

    /// Number of outer-shell electrons from the element table.
    fn outer_electrons(&self) -> i32 {
        self.inner.outer_electrons
    }

    /// Allowed valence values from the element table.
    fn valences(&self) -> Vec<i32> {
        self.inner.valences.to_vec()
    }

    /// RDKit's source-defined ``Rb0`` bond radius in angstroms.
    fn rb0(&self) -> f64 {
        self.inner.rb0
    }

    /// Average atomic weight from the element table.
    fn atomic_weight(&self) -> f64 {
        self.inner.atomic_weight
    }

    fn __repr__(&self) -> String {
        format!(
            "ElementInfo(symbol='{}', atomic_number={}, period={})",
            self.inner.symbol, self.inner.atomic_number, self.inner.period
        )
    }
}

/// Returns RDKit periodic-table metadata for an element, including the dummy (``*``).
///
/// All elements with atomic numbers ``0..=118`` are supported. The record contains
/// the canonical symbol, period, outer electron count, complete ordered valence
/// list, ``Rb0`` bond radius in angstroms, and atomic weight. Source-defined zeros
/// and the unrestricted-valence sentinel ``-1`` are preserved.
///
/// The returned symbol and valence slice borrow immutable shared table data.
/// This query has no options or errors and does not change any molecule.
#[cfg_attr(feature = "stubgen", pyo3_stub_gen::derive::gen_stub_pyfunction)]
#[pyfunction]
fn element_info(element: &Element) -> ElementInfo {
    ElementInfo {
        inner: ck::element_info(element.inner),
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<Element>()?;
    module.add_class::<ElementInfo>()?;
    let class = module.py().get_type::<Element>();
    for inner in ck::Element::iter_with_dummy() {
        let name = if inner == ck::Element::DUMMY {
            "DUMMY".to_owned()
        } else {
            inner.symbol().to_ascii_uppercase()
        };
        class.setattr(name, Element { inner })?;
    }
    module.add_function(wrap_pyfunction!(element_info, module)?)?;
    Ok(())
}
