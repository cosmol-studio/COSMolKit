use crate::chemical_values::enum_value;
use crate::native::Molecule;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::collections::BTreeMap;

enum_value!(CoordinateDimension, { TwoD, ThreeD });

#[pyclass(module = "cosmolkit", eq, eq_int, from_py_object)]
#[derive(Clone, Copy, PartialEq)]
pub enum ConformerErrorKind {
    Unsupported,
    CannotNormalizeZeroLengthVector,
    GenerationFailed,
    InvalidEmbedParametersJson,
}

#[pyclass(module="cosmolkit", extends=pyo3::exceptions::PyValueError)]
pub struct ConformerError {
    inner: ck::ConformerError,
}
#[pymethods]
impl ConformerError {
    #[getter]
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
    fn domain(&self) -> &'static str {
        "conformer"
    }
    // These unnamed tuple payloads contain source messages, not nested errors.
    #[getter]
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
impl ConformerError {
    pub(crate) fn from_rust(inner: ck::ConformerError) -> PyErr {
        Python::attach(|py| match Py::new(py, Self { inner }) {
            Ok(error) => PyErr::from_value(error.into_bound(py).into_any()),
            Err(error) => error,
        })
    }
}

type CoordinateMap = crate::parameters::Dict<crate::parameters::Scalar<i32>, [f64; 3]>;
type CpciMap =
    crate::parameters::Dict<crate::parameters::Scalar<(u32, u32)>, crate::parameters::Scalar<f64>>;

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
    ($($field:ident: $ty:ty),+ $(,)?) => {
        embedding_parameters!(fields {$($field: $ty),+},
            presets { dg, kdg, etdg, etdg_v2, etkdg, etkdg_v2, etkdg_v3, sr_etkdg_v3 });
    };
    (fields {$($field:ident: $ty:ty),+ $(,)?}, presets {$($preset:ident),+ $(,)?}) => {
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
                #[staticmethod] fn new(py: Python<'_>) -> PyResult<Self> { Self::from_rust(py, ck::EmbedParams::new()) }
                $(#[staticmethod] fn $preset(py: Python<'_>) -> PyResult<Self> { Self::from_rust(py, ck::EmbedParams::$preset()) })+
                $(#[getter] fn $field(&self) -> $ty { self.native.$field.clone() }
                #[setter] fn [<set_ $field>](&mut self, value: $ty) { self.native.$field = value; })+
                #[getter] fn coord_map(&self) -> crate::parameters::ContainerOutput<'_, Option<CoordinateMap>> { crate::parameters::ContainerOutput(&self.native.coord_map) }
                #[setter] fn set_coord_map(&mut self, value: crate::parameters::ContainerInput<Option<CoordinateMap>>) -> PyResult<()> {
                    self.native.coord_map = value.0;
                    Ok(())
                }
                #[getter] fn cpci(&self) -> crate::parameters::ContainerOutput<'_, Option<CpciMap>> { crate::parameters::ContainerOutput(&self.native.cpci) }
                #[setter] fn set_cpci(&mut self, value: crate::parameters::ContainerInput<Option<CpciMap>>) -> PyResult<()> {
                    self.native.cpci = value.0;
                    Ok(())
                }
                fn to_json(&self, py: Python<'_>) -> PyResult<String> { Ok(self.to_rust(py)?.to_json()) }
                fn with_json(&self, py: Python<'_>, json: &str) -> PyResult<Self> {
                    let native = self.to_rust(py)?.with_json(json).map_err(ConformerError::from_rust)?;
                    Self::from_rust(py, native)
                }
                #[getter] fn failures(&self) -> Vec<u32> { self.native.failures().clone() }
            }
        }
    };
}
ck::project_configuration!(ck::EmbedParams, using embedding_parameters);

#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub struct EmbedMoleculeResult {
    pub(crate) inner: ck::EmbedMoleculeResult,
}
#[pymethods]
impl EmbedMoleculeResult {
    #[new]
    fn new(
        py: Python<'_>,
        molecule: &Molecule,
        conf_id: i32,
        params: &EmbedParams,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: ck::EmbedMoleculeResult {
                molecule: molecule.inner.clone(),
                conf_id,
                params: params.to_rust(py)?,
            },
        })
    }
    fn molecule(&self) -> Molecule {
        Molecule {
            inner: self.inner.molecule().clone(),
        }
    }
    fn conf_id(&self) -> i32 {
        self.inner.conf_id()
    }
    fn params(&self, py: Python<'_>) -> PyResult<EmbedParams> {
        EmbedParams::from_rust(py, self.inner.params().clone())
    }
    fn ok(&self) -> bool {
        self.inner.ok()
    }
}
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub struct EmbedMultipleConfsResult {
    pub(crate) inner: ck::EmbedMultipleConfsResult,
}
#[pymethods]
impl EmbedMultipleConfsResult {
    #[new]
    fn new(
        py: Python<'_>,
        molecule: &Molecule,
        conf_ids: Vec<i32>,
        requested_num_confs: u32,
        params: &EmbedParams,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: ck::EmbedMultipleConfsResult {
                molecule: molecule.inner.clone(),
                conf_ids,
                requested_num_confs,
                params: params.to_rust(py)?,
            },
        })
    }
    fn molecule(&self) -> Molecule {
        Molecule {
            inner: self.inner.molecule().clone(),
        }
    }
    fn conf_ids(&self) -> &[i32] {
        self.inner.conf_ids()
    }
    fn requested_num_confs(&self) -> u32 {
        self.inner.requested_num_confs()
    }
    fn generated_count(&self) -> usize {
        self.inner.generated_count()
    }
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
            let options = result.bind(py).call_method0("params")?;
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
