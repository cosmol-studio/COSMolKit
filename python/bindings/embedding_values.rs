use crate::chemical_values::{EnumValue, enum_value};
use crate::native::Molecule;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::collections::BTreeMap;
use std::error::Error;

enum_value!(CoordinateDimension, { TwoD, ThreeD });
enum_value!(EmbedTermination, { Completed, TimedOut, Interrupted });

enum_value!(EmbedFailureCause, {
    InitialCoords, FirstMinimization, CheckTetrahedralCenters, CheckChiralCenters,
    MinimizeFourthDimension, EtkMinimization, FinalChiralBounds, FinalCenterInVolume,
    LinearDoubleBond, BadDoubleBondStereo, CheckChiralCenters2, ExceededTimeout,
    Minimization, KtermViolation, Clash, EndOfEnum
});

#[derive(FromPyObject)]
pub enum EmbedFailureCauseInput {
    Enum(EmbedFailureCause),
    Name(String),
}
impl EmbedFailureCauseInput {
    fn into_rust(self) -> PyResult<ck::EmbedFailureCause> {
        match self {
            Self::Enum(value) => Ok(value.into()),
            Self::Name(name) => EmbedFailureCause::from_name(&name).map(Into::into),
        }
    }
}

crate::error_framework::error_projection!(
    EmbedFailureCountError, EmbedFailureCountErrorKind, ck::EmbedFailureCountError, "conformer", {
        InvalidCause,
        MissingCounter { index: usize, available: usize },
    }, cause = |_: &ck::EmbedFailureCountError| None::<PyErr>
);

pub struct FailureCounts<'a>(&'a [u32]);
impl<'py> IntoPyObject<'py> for FailureCounts<'_> {
    type Target = pyo3::types::PyTuple;
    type Output = Bound<'py, Self::Target>;
    type Error = PyErr;
    const OUTPUT_TYPE: pyo3::inspect::PyStaticExpr =
        <Vec<u32> as crate::parameters::Container>::TYPE;
    fn into_pyobject(self, py: Python<'py>) -> PyResult<Self::Output> {
        pyo3::types::PyTuple::new(py, self.0.iter().copied())
    }
}

crate::enums::binding_enum!(ConformerErrorKind, {
    Unsupported,
    CannotNormalizeZeroLengthVector,
    GenerationFailed,
    InvalidEmbedParametersJson,
});

crate::error_framework::error_projection! {record
pub struct ConformerError {
    inner: ck::ConformerError,
}
}
crate::error_framework::error_projection! {methods impl ConformerError {
    #[getter]
    /// Concrete native failure category; payload fields apply only to their corresponding variants.
    fn variant(&self) -> ConformerErrorKind {
        match &self.inner {
            ck::ConformerError::Unsupported => ConformerErrorKind::Unsupported,
            ck::ConformerError::CannotNormalizeZeroLengthVector => {
                ConformerErrorKind::CannotNormalizeZeroLengthVector
            }
            ck::ConformerError::GenerationFailed(_) => ConformerErrorKind::GenerationFailed,
            ck::ConformerError::InvalidEmbedParametersJson(_) => {
                ConformerErrorKind::InvalidEmbedParametersJson
            }
        }
    }
    #[getter]
    /// Native error domain identifying the conformer subsystem.
    fn domain(&self) -> &'static str {
        "conformer"
    }
    // These unnamed tuple payloads contain source messages, not nested errors.
    #[getter]
    /// Native failure message for generation or JSON parameter errors; None for variants without a message.
    fn message(&self) -> Option<&str> {
        match &self.inner {
            ck::ConformerError::GenerationFailed(message)
            | ck::ConformerError::InvalidEmbedParametersJson(message) => Some(message),
            ck::ConformerError::Unsupported
            | ck::ConformerError::CannotNormalizeZeroLengthVector => None,
        }
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
}
impl ConformerError {
    pub(crate) fn from_rust(inner: ck::ConformerError) -> PyErr {
        Python::attach(|py| match Py::new(py, Self { inner }) {
            Ok(error) => crate::error_framework::from_payload(error),
            Err(error) => error,
        })
    }
}

type CoordinateMap = crate::parameters::Dict<crate::parameters::Scalar<i32>, [f64; 3]>;
type CpciMap =
    crate::parameters::Dict<crate::parameters::Scalar<(u32, u32)>, crate::parameters::Scalar<f64>>;

/// Distance-geometry embedding options with native presets and optional coordinate constraints.
/// Set fields before execution; result objects carry execution diagnostics separately.
#[pyclass(module = "cosmolkit")]
pub struct EmbedParams {
    // The native owner retains input fields, private bounds and diagnostics.
    native: ck::EmbedParams,
}
impl EmbedParams {
    pub(crate) fn from_rust(_py: Python<'_>, native: ck::EmbedParams) -> PyResult<Self> {
        Ok(Self { native })
    }
    pub(crate) fn to_rust(&self, _py: Python<'_>) -> PyResult<ck::EmbedParams> {
        Ok(self.native.clone())
    }
}
macro_rules! embedding_parameters {
    ($($(#[$doc:meta])* $field:ident: $ty:ty),+ $(,)?) => {
        embedding_parameters!(fields {$($(#[$doc])* $field: $ty),+},
            presets { dg, kdg, etdg, etdg_v2, etkdg, etkdg_v2, etkdg_v3, sr_etkdg_v3 });
    };
    (fields {$($(#[$doc:meta])* $field:ident: $ty:ty),+ $(,)?}, presets {$($preset:ident),+ $(,)?}) => {
        paste::paste! {
            #[pymethods]
            impl EmbedParams {
                #[new]
                #[pyo3(signature=(*, $($field=ck::EmbedParams::default().$field,)+ coord_map=crate::parameters::FieldInput::Default, cpci=crate::parameters::FieldInput::Default))]
                fn construct(py: Python<'_>, $($field: $ty,)+ coord_map: crate::parameters::FieldInput<Option<CoordinateMap>>, cpci: crate::parameters::FieldInput<Option<CpciMap>>) -> PyResult<Self> {
                    let mut native = ck::EmbedParams::new();
                    $(native.$field = $field;)+
                    native.coord_map = coord_map.or_default(|| native.coord_map.take());
                    native.cpci = cpci.or_default(|| native.cpci.take());
                    Self::from_rust(py, native)
                }
                /// Create distance-geometry parameters using the native constructor defaults.
                #[staticmethod] fn new(py: Python<'_>) -> PyResult<Self> { Self::from_rust(py, ck::EmbedParams::new()) }
                $(#[doc = concat!("Create parameters using the native ", stringify!($preset), " embedding preset; subsequent field edits affect this configuration only.")]
                #[staticmethod] fn $preset(py: Python<'_>) -> PyResult<Self> { Self::from_rust(py, ck::EmbedParams::$preset()) })+
                $($(#[$doc])* #[getter] fn $field(&self) -> $ty { self.native.$field.clone() }
                #[setter] fn [<set_ $field>](&mut self, value: $ty) { self.native.$field = value; })+
                /// Readonly atom-index to three-dimensional coordinate constraints, or None. Replace the whole field to edit.
                #[getter] fn coord_map(&self) -> crate::parameters::ContainerOutput<'_, Option<CoordinateMap>> { crate::parameters::ContainerOutput(&self.native.coord_map) }
                #[setter] fn set_coord_map(&mut self, value: crate::parameters::ContainerInput<Option<CoordinateMap>>) -> PyResult<()> {
                    self.native.coord_map = value.0;
                    Ok(())
                }
                /// Readonly atom-pair to custom Coulomb-like interaction strengths, or None. Replace the whole field to edit.
                #[getter] fn cpci(&self) -> crate::parameters::ContainerOutput<'_, Option<CpciMap>> { crate::parameters::ContainerOutput(&self.native.cpci) }
                #[setter] fn set_cpci(&mut self, value: crate::parameters::ContainerInput<Option<CpciMap>>) -> PyResult<()> {
                    self.native.cpci = value.0;
                    Ok(())
                }
                /// Serialize the native JSON-supported embedding options; this is not a complete execution-state archive.
                fn to_json(&self, py: Python<'_>) -> PyResult<String> { Ok(self.to_rust(py)?.to_json()) }
                /// Return an independent configuration updated from JSON using the native parameter parser.
                fn with_json(&self, py: Python<'_>, json: &str) -> PyResult<Self> {
                    let native = self.to_rust(py)?.with_json(json).map_err(ConformerError::from_rust)?;
                    Self::from_rust(py, native)
                }
                /// Copy of recorded failure-stage counters; execution results also expose failure_counts directly.
                #[getter] fn failures(&self) -> Vec<u32> { self.native.failures().clone() }
            }
        }
    };
}
ck::project_configuration!(ck::EmbedParams, using embedding_parameters);

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// Single-conformer embedding result retaining the output molecule, optional generated ID and execution diagnostics.
    pub struct EmbedMoleculeResult {
        pub(crate) inner: ck::EmbedMoleculeResult,
    }
}
#[pymethods]
impl EmbedMoleculeResult {
    #[getter]
    /// Reason embedding stopped: completed, timed out or interrupted; inspect generated IDs separately.
    fn termination(&self) -> EmbedTermination {
        self.inner.termination().into()
    }
    #[getter]
    /// Whether failure-stage counters were enabled for this run.
    fn failures_tracked(&self) -> bool {
        self.inner.failures_tracked()
    }
    #[getter]
    /// Readonly counts in EmbedFailureCause order, or None when tracking was disabled.
    fn failure_counts(&self) -> Option<FailureCounts<'_>> {
        self.inner.failure_counts().map(FailureCounts)
    }
    /// Return the count for one failure cause, or None when tracking was disabled; invalid causes raise EmbedFailureCountError.
    fn failure_count(&self, cause: EmbedFailureCauseInput) -> PyResult<Option<u32>> {
        self.inner
            .failure_count(cause.into_rust()?)
            .map_err(EmbedFailureCountError::from_rust)
    }
    #[new]
    #[pyo3(signature = (molecule, conf_id, params, termination=EmbedTermination::Completed))]
    fn new(
        py: Python<'_>,
        molecule: &Molecule,
        conf_id: Option<u32>,
        params: &EmbedParams,
        termination: EmbedTermination,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: ck::EmbedMoleculeResult {
                molecule: molecule.inner.clone(),
                conf_id,
                termination: termination.into(),
                params: params.to_rust(py)?,
            },
        })
    }
    #[getter]
    /// Output molecule as a copy-on-write value; editing it does not alter this result.
    fn molecule(&self) -> Molecule {
        Molecule {
            inner: self.inner.molecule().clone(),
        }
    }
    #[getter]
    /// Generated conformer identifier, or None when no conformer was generated.
    fn conf_id(&self) -> Option<u32> {
        self.inner.conf_id()
    }
    #[getter]
    /// Independent configuration value carrying this run's diagnostics; modifying it does not change the result.
    fn params(&self, py: Python<'_>) -> PyResult<EmbedParams> {
        EmbedParams::from_rust(py, self.inner.params().clone())
    }
    #[getter]
    /// Whether a conformer identifier was produced by this single-conformer run.
    fn ok(&self) -> bool {
        self.inner.ok()
    }
}
crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// Multi-conformer embedding result retaining successful IDs, the requested count and execution diagnostics.
    pub struct EmbedMultipleConfsResult {
        pub(crate) inner: ck::EmbedMultipleConfsResult,
    }
}
#[pymethods]
impl EmbedMultipleConfsResult {
    #[getter]
    /// Reason embedding stopped: completed, timed out or interrupted; inspect generated IDs separately.
    fn termination(&self) -> EmbedTermination {
        self.inner.termination().into()
    }
    #[getter]
    /// Whether failure-stage counters were enabled for this run.
    fn failures_tracked(&self) -> bool {
        self.inner.failures_tracked()
    }
    #[getter]
    /// Readonly counts in EmbedFailureCause order, or None when tracking was disabled.
    fn failure_counts(&self) -> Option<FailureCounts<'_>> {
        self.inner.failure_counts().map(FailureCounts)
    }
    /// Return the count for one failure cause, or None when tracking was disabled; invalid causes raise EmbedFailureCountError.
    fn failure_count(&self, cause: EmbedFailureCauseInput) -> PyResult<Option<u32>> {
        self.inner
            .failure_count(cause.into_rust()?)
            .map_err(EmbedFailureCountError::from_rust)
    }
    #[new]
    #[pyo3(signature = (molecule, conf_ids, requested_num_confs, params, termination=EmbedTermination::Completed))]
    fn new(
        py: Python<'_>,
        molecule: &Molecule,
        conf_ids: Vec<u32>,
        requested_num_confs: u32,
        params: &EmbedParams,
        termination: EmbedTermination,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: ck::EmbedMultipleConfsResult {
                molecule: molecule.inner.clone(),
                conf_ids,
                termination: termination.into(),
                requested_num_confs,
                params: params.to_rust(py)?,
            },
        })
    }
    #[getter]
    /// Output molecule as a copy-on-write value; editing it does not alter this result.
    fn molecule(&self) -> Molecule {
        Molecule {
            inner: self.inner.molecule().clone(),
        }
    }
    #[getter]
    /// Identifiers of successfully generated conformers; failed attempts are not conformer IDs.
    fn conf_ids(&self) -> &[u32] {
        self.inner.conf_ids()
    }
    #[getter]
    /// Number of conformers requested, which may exceed the generated count after pruning or failure.
    fn requested_num_confs(&self) -> u32 {
        self.inner.requested_num_confs()
    }
    #[getter]
    /// Number of successfully generated conformer IDs, not the number requested or attempted.
    fn generated_count(&self) -> usize {
        self.inner.generated_count()
    }
    #[getter]
    /// Independent configuration value carrying this run's diagnostics; modifying it does not change the result.
    fn params(&self, py: Python<'_>) -> PyResult<EmbedParams> {
        EmbedParams::from_rust(py, self.inner.params().clone())
    }
}

#[cfg(test)]
mod diagnostic_tests {
    use super::*;

    #[test]
    fn actual_embedding_diagnostics_survive_result_parameter_projection() {
        let molecule = ck::Molecule::from_smiles("CCO").unwrap();
        let mut options = ck::EmbedParams::etkdg_v3();
        options.random_seed = 17;
        options.track_failures = true;
        let result = molecule
            .with_3d_conformer_result_with_params(&options)
            .unwrap();
        let expected = result.params().failures().clone();
        assert!(!expected.is_empty());
        assert!(options.failures().is_empty());
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let result = Py::new(py, EmbedMoleculeResult { inner: result })?;
            let options = result.bind(py).getattr("params")?;
            assert_eq!(
                options.getattr("failures")?.extract::<Vec<u32>>()?,
                expected
            );
            assert!(options.getattr("track_failures")?.extract::<bool>()?);
            assert!(
                options
                    .setattr("failures", vec![1_u32])
                    .unwrap_err()
                    .is_instance_of::<pyo3::exceptions::PyAttributeError>(py)
            );
            let copy = options.getattr("failures")?;
            copy.call_method1("append", (u32::MAX,))?;
            assert_eq!(
                options.getattr("failures")?.extract::<Vec<u32>>()?,
                expected
            );
            drop(result);
            assert_eq!(
                options.getattr("failures")?.extract::<Vec<u32>>()?,
                expected
            );
            Ok(())
        })
        .unwrap();
    }
}
