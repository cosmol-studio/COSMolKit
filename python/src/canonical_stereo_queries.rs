//! Thin projections of the project-native read-only stereo values.
use ::cosmolkit as ck;
use pyo3::prelude::*;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

pyo3::create_exception!(cosmolkit, StereoReadError, pyo3::exceptions::PyValueError);

pub(crate) fn error_pyerr(py: Python<'_>, source: ck::StereoReadError) -> PyErr {
    let kind = match &source {
        ck::StereoReadError::InvalidTopology(_) => "InvalidTopology",
    };
    crate::canonical_values::annotate(
        py,
        StereoReadError::new_err(source.to_string()),
        "stereo",
        kind,
        &source,
    )
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct LigandRef {
    inner: ck::LigandRef,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl LigandRef {
    #[new]
    #[pyo3(signature = (atom=None))]
    fn new(atom: Option<usize>) -> Self {
        Self {
            inner: atom.map_or(ck::LigandRef::ImplicitHydrogen, |atom| {
                ck::LigandRef::Atom(ck::AtomId::new(atom))
            }),
        }
    }
    #[getter]
    fn atom(&self) -> Option<usize> {
        match self.inner {
            ck::LigandRef::Atom(atom) => Some(atom.index()),
            ck::LigandRef::ImplicitHydrogen => None,
        }
    }
    #[getter]
    fn is_implicit_hydrogen(&self) -> bool {
        matches!(self.inner, ck::LigandRef::ImplicitHydrogen)
    }
}

pub(crate) fn tetrahedral_row(py: Python<'_>, row: ck::TetrahedralStereo) -> PyResult<Py<PyAny>> {
    let ligands = row
        .ligands
        .into_iter()
        .map(|ligand| match ligand {
            ck::LigandRef::Atom(atom) => Some(atom.index()),
            ck::LigandRef::ImplicitHydrogen => None,
        })
        .collect::<Vec<_>>();
    Ok(py
        .import("cosmolkit")?
        .getattr("TetrahedralStereo")?
        .call1((row.center.index(), ligands))?
        .unbind())
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<LigandRef>()?;
    module.add("StereoReadError", module.py().get_type::<StereoReadError>())?;
    // Named records retain the original tuple/index/unpacking/equality protocol
    // while exposing the canonical result name and center/ligands fields.
    let kwargs = pyo3::types::PyDict::new(module.py());
    kwargs.set_item("module", "cosmolkit")?;
    let value = module
        .py()
        .import("collections")?
        .getattr("namedtuple")?
        .call(("TetrahedralStereo", ["center", "ligands"]), Some(&kwargs))?;
    module.add("TetrahedralStereo", value)?;
    Ok(())
}
