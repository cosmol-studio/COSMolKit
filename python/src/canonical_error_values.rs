//! Lossless projections of canonical error categories and source stages.
use ::cosmolkit as ck;
use pyo3::prelude::*;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pyclass_enum, gen_stub_pymethods};
/// Reading phase at which PDB input failed: stream access, record parsing or finalization.
///
/// Declared values: ``Stream``, ``Record``, ``Finalization``.
#[cosmolkit_macros::python_enum]
#[cfg_attr(feature = "stubgen", gen_stub_pyclass_enum)]
#[pyclass(module = "cosmolkit", eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum BioPdbReadStage {
    Stream,
    Record,
    Finalization,
}
impl From<ck::BioPdbReadStage> for BioPdbReadStage {
    fn from(value: ck::BioPdbReadStage) -> Self {
        match value {
            ck::BioPdbReadStage::Stream => Self::Stream,
            ck::BioPdbReadStage::Record => Self::Record,
            ck::BioPdbReadStage::Finalization => Self::Finalization,
        }
    }
}
/// Stage of structural mmCIF reading that produced a typed error.
///
/// Declared values: ``CifDocument``, ``CoordinateBlock``, ``CrystalCell``, ``Refinement``, ``Tls``, ``Experimental``, ``Reflections``, ``Software``, ``Ncs``, ``FractionalTransform``, ``Origx``, ``AnisotropicU``, ``AtomSites``, ``EntitySequence``, ``Helices``, ``Sheets``, ``Connections``, ``CisPeptides``, ``ModifiedResidues``, ``Assemblies``, ``SiftsUnp``, ``CcdRestoration``, ``Materialization``, ``StructureValidation``.
#[cosmolkit_macros::python_enum]
#[cfg_attr(feature = "stubgen", gen_stub_pyclass_enum)]
#[pyclass(module = "cosmolkit", eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum BioMmcifReadStage {
    CifDocument,
    CoordinateBlock,
    CrystalCell,
    Refinement,
    Tls,
    Experimental,
    Reflections,
    Software,
    Ncs,
    FractionalTransform,
    Origx,
    AnisotropicU,
    AtomSites,
    EntitySequence,
    Helices,
    Sheets,
    Connections,
    CisPeptides,
    ModifiedResidues,
    Assemblies,
    SiftsUnp,
    CcdRestoration,
    Materialization,
    StructureValidation,
}
impl From<ck::BioMmcifReadStage> for BioMmcifReadStage {
    fn from(value: ck::BioMmcifReadStage) -> Self {
        match value {
            ck::BioMmcifReadStage::CifDocument => Self::CifDocument,
            ck::BioMmcifReadStage::CoordinateBlock => Self::CoordinateBlock,
            ck::BioMmcifReadStage::CrystalCell => Self::CrystalCell,
            ck::BioMmcifReadStage::Refinement => Self::Refinement,
            ck::BioMmcifReadStage::Tls => Self::Tls,
            ck::BioMmcifReadStage::Experimental => Self::Experimental,
            ck::BioMmcifReadStage::Reflections => Self::Reflections,
            ck::BioMmcifReadStage::Software => Self::Software,
            ck::BioMmcifReadStage::Ncs => Self::Ncs,
            ck::BioMmcifReadStage::FractionalTransform => Self::FractionalTransform,
            ck::BioMmcifReadStage::Origx => Self::Origx,
            ck::BioMmcifReadStage::AnisotropicU => Self::AnisotropicU,
            ck::BioMmcifReadStage::AtomSites => Self::AtomSites,
            ck::BioMmcifReadStage::EntitySequence => Self::EntitySequence,
            ck::BioMmcifReadStage::Helices => Self::Helices,
            ck::BioMmcifReadStage::Sheets => Self::Sheets,
            ck::BioMmcifReadStage::Connections => Self::Connections,
            ck::BioMmcifReadStage::CisPeptides => Self::CisPeptides,
            ck::BioMmcifReadStage::ModifiedResidues => Self::ModifiedResidues,
            ck::BioMmcifReadStage::Assemblies => Self::Assemblies,
            ck::BioMmcifReadStage::SiftsUnp => Self::SiftsUnp,
            ck::BioMmcifReadStage::CcdRestoration => Self::CcdRestoration,
            ck::BioMmcifReadStage::Materialization => Self::Materialization,
            ck::BioMmcifReadStage::StructureValidation => Self::StructureValidation,
        }
    }
}
/// Category of UFF parameterization failure: preparation, parameter-table lookup or atom typing.
///
/// Declared values: ``Preparation``, ``ParameterTable``, ``Typing``.
#[cosmolkit_macros::python_enum]
#[cfg_attr(feature = "stubgen", gen_stub_pyclass_enum)]
#[pyclass(module = "cosmolkit", eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum UffParameterErrorKind {
    Preparation,
    ParameterTable,
    Typing,
}
impl From<ck::UffParameterErrorKind> for UffParameterErrorKind {
    fn from(value: ck::UffParameterErrorKind) -> Self {
        match value {
            ck::UffParameterErrorKind::Preparation => Self::Preparation,
            ck::UffParameterErrorKind::ParameterTable => Self::ParameterTable,
            ck::UffParameterErrorKind::Typing => Self::Typing,
        }
    }
}
/// Stable categories; unavailable geometry is not a guessed 2D fallback.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct UffOptimizationErrorKind {
    pub(crate) inner: ck::UffOptimizationErrorKind,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl UffOptimizationErrorKind {
    /// Named variant/category of the result.
    #[getter]
    fn variant(&self) -> &'static str {
        match self.inner {
            ck::UffOptimizationErrorKind::MissingConformer { .. } => "MissingConformer",
            ck::UffOptimizationErrorKind::Rings => "Rings",
            ck::UffOptimizationErrorKind::Optimization => "Optimization",
            ck::UffOptimizationErrorKind::ConformerOptimization => "ConformerOptimization",
            ck::UffOptimizationErrorKind::Evaluation => "Evaluation",
        }
    }
    /// Requested conformer identifier in a missing-conformer error, when available.
    #[getter]
    fn requested(&self) -> Option<usize> {
        match self.inner {
            ck::UffOptimizationErrorKind::MissingConformer { requested } => requested,
            _ => None,
        }
    }
    fn __eq__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        match other.extract::<PyRef<'_, Self>>() {
            Ok(other) => Ok((self.inner == other.inner)
                .into_pyobject(py)?
                .to_owned()
                .into_any()
                .unbind()),
            Err(_) => Ok(py.NotImplemented()),
        }
    }
    fn __repr__(&self) -> String {
        format!("UffOptimizationErrorKind({:?})", self.inner)
    }
}
