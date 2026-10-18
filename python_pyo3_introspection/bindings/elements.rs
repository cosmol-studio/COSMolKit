//! Element values and periodic metadata from the public Rust facade.
use cosmolkit as ck;
use pyo3::prelude::*;

#[pyclass(module = "fixed_projection._native", frozen, from_py_object)]
#[derive(Clone)]
pub struct Element {
    pub(crate) inner: ck::Element,
}

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
}

#[pyclass(module = "fixed_projection._native", frozen)]
pub struct ElementInfo {
    inner: ck::ElementInfo,
}

#[pymethods]
impl ElementInfo {
    #[getter]
    fn element(&self) -> Element {
        Element {
            inner: self.inner.element,
        }
    }
    #[getter]
    fn symbol(&self) -> &'static str {
        self.inner.symbol
    }
    #[getter]
    fn atomic_number(&self) -> u8 {
        self.inner.atomic_number
    }
    #[getter]
    fn period(&self) -> u8 {
        self.inner.period
    }
    #[getter]
    fn outer_electrons(&self) -> i32 {
        self.inner.outer_electrons
    }
    #[getter]
    fn valences(&self) -> Vec<i32> {
        self.inner.valences.to_vec()
    }
    #[getter]
    fn rb0(&self) -> f64 {
        self.inner.rb0
    }
    #[getter]
    fn atomic_weight(&self) -> f64 {
        self.inner.atomic_weight
    }
}

#[pyfunction]
pub fn element_info(element: &Element) -> ElementInfo {
    ElementInfo {
        inner: ck::element_info(element.inner),
    }
}
