//! Immutable element values and result-field projections of the public facade.

use ::cosmolkit as ck;
use pyo3::prelude::*;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

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
    #[staticmethod]
    fn from_atomic_number(atomic_number: u8) -> Option<Self> {
        ck::Element::from_atomic_number(atomic_number).map(|inner| Self { inner })
    }

    #[staticmethod]
    fn from_symbol(symbol: &str) -> Option<Self> {
        ck::Element::from_symbol(symbol).map(|inner| Self { inner })
    }

    fn atomic_number(&self) -> u8 {
        self.inner.atomic_number()
    }

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

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct ElementInfo {
    inner: ck::ElementInfo,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl ElementInfo {
    // Read-callable transport of the existing eight public result fields.
    fn element(&self) -> Element {
        Element {
            inner: self.inner.element,
        }
    }

    fn symbol(&self) -> &'static str {
        self.inner.symbol
    }

    fn atomic_number(&self) -> u8 {
        self.inner.atomic_number
    }

    fn period(&self) -> u8 {
        self.inner.period
    }

    fn outer_electrons(&self) -> i32 {
        self.inner.outer_electrons
    }

    fn valences(&self) -> Vec<i32> {
        self.inner.valences.to_vec()
    }

    fn rb0(&self) -> f64 {
        self.inner.rb0
    }

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
