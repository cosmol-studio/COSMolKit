//! Native ordered-value collections. Business-specific exports stay in bindings.
use pyo3::prelude::*;

pub(crate) trait Projection: Sized {
    type Wrapper: pyo3::PyClass + From<Self> + Into<PyClassInitializer<Self::Wrapper>>;
}

/// Return conversion resolves from the native result type, not a method table.
pub(crate) struct Output<T: Projection>(pub T);
impl<'py, T: Projection> IntoPyObject<'py> for Output<T> {
    type Target = T::Wrapper;
    type Output = Bound<'py, T::Wrapper>;
    type Error = PyErr;
    const OUTPUT_TYPE: pyo3::inspect::PyStaticExpr =
        <Py<T::Wrapper> as IntoPyObject<'py>>::OUTPUT_TYPE;
    fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
        Bound::new(py, T::Wrapper::from(self.0))
    }
}

macro_rules! project_sequence {
    ($wrapper:ident, $native:ty, $item:ty, { $($methods:tt)* }) => {
        #[doc = concat!("Ordered, immutable collection of optional ", stringify!($item),
            " values.\n\nNone retains a missing result at its original position. ",
            "len() counts all records. Integer indexing returns a record; slices, ",
            "index sequences and boolean masks return a new collection. ",
            "Selections preserve order and duplicate indices.")]
        #[pyo3::pyclass(module = "cosmolkit", frozen)]
        pub struct $wrapper {
            pub(crate) inner: $native,
            // Immutable object identities are cached; Rust remains the sole value storage.
            objects: std::sync::Mutex<Vec<Option<pyo3::Py<$item>>>>,
        }
        impl From<$native> for $wrapper {
            fn from(inner: $native) -> Self {
                let objects = std::iter::repeat_with(|| None).take(inner.len()).collect();
                Self { inner, objects: std::sync::Mutex::new(objects) }
            }
        }
        impl $crate::sequence::Projection for $native { type Wrapper = $wrapper; }

        impl $crate::indexing::Collection for $wrapper {
            type Item = Option<pyo3::Py<$item>>;
            type Selected = pyo3::Py<Self>;
            fn length(&self) -> usize { self.inner.len() }
            fn item(&self, py: pyo3::Python<'_>, at: usize) -> pyo3::PyResult<Self::Item> {
                let mut objects = self.objects.lock().map_err(|error| pyo3::exceptions::PyRuntimeError::new_err(error.to_string()))?;
                if let Some(inner) = &self.inner.values()[at] {
                    if objects[at].is_none() {
                        objects[at] = Some(pyo3::Py::new(py, <$item as $crate::parameters::Configuration>::from_native(py, inner.clone())?)?);
                    }
                }
                Ok(objects[at].as_ref().map(|v| v.clone_ref(py)))
            }
            fn select(&self, py: pyo3::Python<'_>, selected: &[usize]) -> pyo3::PyResult<Self::Selected> {
                let inner = self.inner.select(selected).ok_or_else(|| pyo3::exceptions::PyIndexError::new_err("invalid selection"))?;
                let objects = self.objects.lock().map_err(|error| pyo3::exceptions::PyRuntimeError::new_err(error.to_string()))?;
                let objects = selected.iter().map(|&i| objects[i].as_ref().map(|v| v.clone_ref(py))).collect();
                pyo3::Py::new(py, Self { inner, objects: std::sync::Mutex::new(objects) })
            }
        }
        #[pyo3::pymethods]
        impl $wrapper {
            $($methods)*
            #[new]
            fn new(py: pyo3::Python<'_>, values: Vec<Option<pyo3::Py<$item>>>) -> pyo3::PyResult<Self> {
                let rows = values.iter().map(|item| item.as_ref().map(|item| {
                    <$item as $crate::parameters::Configuration>::to_native(&*item.try_borrow(py)?, py)
                }).transpose()).collect::<pyo3::PyResult<Vec<_>>>()?;
                let inner = <$native>::new(rows);
                Ok(Self { inner, objects: std::sync::Mutex::new(values) })
            }
            /// Number of records, including missing (None) records.
            #[getter]
            fn len(&self) -> usize { self.inner.len() }
            fn __len__(&self) -> usize { self.len() }
            fn __repr__(&self) -> String { format!("{}(rows={})", stringify!($wrapper), self.inner.len()) }
            #[doc = r#"{"role":"indexing","protocol":"integer_selection","selection":"Self"}"#]
            /// Read a record or select records without changing this collection.
            /// Negative indices and Python slices are supported. Integer selections
            /// preserve duplicates; boolean masks must match the collection length.
            fn __getitem__(&self, py: pyo3::Python<'_>, index: $crate::indexing::Index<'_>)
                -> pyo3::PyResult<$crate::indexing::Indexed<Option<pyo3::Py<$item>>, pyo3::Py<Self>>> {
                $crate::indexing::get(self, py, index)
            }
        }
    };
}
pub(crate) use project_sequence;
