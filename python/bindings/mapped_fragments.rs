use crate::arguments::{Argument, configuration_arguments};
use crate::native::Molecule;
use crate::operation_errors::OperationError;
use crate::params::FragmentExtractionParams;
use cosmolkit as ck;
use pyo3::prelude::*;
use pyo3::types::PyTuple;

crate::configuration::type_projection! {
    value(eq)
    /// A connected component and its local-to-source atom correspondence.
    pub struct MappedFragment { pub(crate) inner: ck::MappedFragment }
}

struct SourceAtomIndices<'a>(&'a [ck::AtomId]);

impl<'py> IntoPyObject<'py> for SourceAtomIndices<'_> {
    type Target = PyTuple;
    type Output = Bound<'py, PyTuple>;
    type Error = PyErr;
    const OUTPUT_TYPE: pyo3::inspect::PyStaticExpr =
        <Vec<usize> as crate::parameters::Container>::TYPE;
    fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
        PyTuple::new(py, self.0.iter().map(|id| id.index()))
    }
}

#[pymethods]
impl MappedFragment {
    /// Independently editable molecule in fragment-local atom order.
    #[getter]
    fn molecule(&self) -> Molecule {
        Molecule {
            inner: self.inner.molecule().clone(),
        }
    }

    /// Entry i identifies local atom i in the input at extraction time.
    #[getter]
    fn source_atom_indices(&self) -> SourceAtomIndices<'_> {
        SourceAtomIndices(self.inner.source_atom_indices())
    }
}

impl Molecule {
    fn extract_mapped_fragments(
        &self,
        py: Python<'_>,
        params: &ck::FragmentExtractionParams,
    ) -> PyResult<Vec<MappedFragment>> {
        py.detach(|| self.inner.fragments_with_atom_mapping_with_params(params))
            .map(|parts| {
                parts
                    .into_iter()
                    .map(|inner| MappedFragment { inner })
                    .collect()
            })
            .map_err(OperationError::from_rust)
    }
}

#[pymethods]
impl Molecule {
    /// Extract connected components with their atom indices in this molecule.
    ///
    /// Parameters
    /// ----------
    /// params : FragmentExtractionParams, optional
    ///     Extraction settings. Supply this object or the keyword option.
    /// sanitize : bool, optional
    ///     Sanitize each component before returning it. Defaults to True.
    ///
    /// Returns
    /// -------
    /// list of MappedFragment
    ///     One record per component in input order, including repeated components.
    ///     Each record's indices are in the same order as its molecule's atoms.
    ///
    /// Raises
    /// ------
    /// OperationError
    ///     If extraction or the requested sanitization fails.
    /// TypeError
    ///     If a parameter object and a keyword option are supplied together.
    #[pyo3(signature = (params=Argument::Missing, *, sanitize=Argument::Missing))]
    fn fragments_with_atom_mapping(
        &self,
        py: Python<'_>,
        params: Argument<Option<PyRef<'_, FragmentExtractionParams>>>,
        sanitize: Argument<bool>,
    ) -> PyResult<Vec<MappedFragment>> {
        let keyword_form = sanitize.is_supplied();
        configuration_arguments!(params; sanitize);
        if keyword_form {
            self._fragments_with_atom_mapping_keywords(py, sanitize)
        } else {
            self._fragments_with_atom_mapping_params(py, params.as_deref())
        }
    }

    #[doc = r#"{"role":"overload","target":"Molecule.fragments_with_atom_mapping"}"#]
    #[pyo3(signature = (params=None))]
    fn _fragments_with_atom_mapping_params(
        &self,
        py: Python<'_>,
        params: Option<&FragmentExtractionParams>,
    ) -> PyResult<Vec<MappedFragment>> {
        let params = params
            .map(|p| p.to_rust(py))
            .transpose()?
            .unwrap_or_default();
        self.extract_mapped_fragments(py, &params)
    }

    #[doc = r#"{"role":"overload","target":"Molecule.fragments_with_atom_mapping"}"#]
    #[pyo3(signature = (*, sanitize))]
    fn _fragments_with_atom_mapping_keywords(
        &self,
        py: Python<'_>,
        sanitize: bool,
    ) -> PyResult<Vec<MappedFragment>> {
        self.extract_mapped_fragments(py, &ck::FragmentExtractionParams { sanitize })
    }

    /// Extract mapped components using an explicit parameter object.
    fn fragments_with_atom_mapping_with_params(
        &self,
        py: Python<'_>,
        params: &FragmentExtractionParams,
    ) -> PyResult<Vec<MappedFragment>> {
        self._fragments_with_atom_mapping_params(py, Some(params))
    }
}
