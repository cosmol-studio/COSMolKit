//! Configuration projects canonical native inputs; it owns no fingerprint logic.
use cosmolkit as ck;
use pyo3::prelude::*;

macro_rules! generator_parameters {
    ($name:ident) => {
        ck::project_configuration!(ck::$name, methods {
            fn info_string(&self, py: Python<'_>) -> PyResult<String> {
                Ok(self.to_rust(py)?.info_string())
            }
            fn to_json(&self, py: Python<'_>) -> PyResult<String> {
                Ok(self.to_rust(py)?.to_json())
            }
            fn with_json(&self, py: Python<'_>, json: &str) -> PyResult<Self> {
                let value = self.to_rust(py)?.with_json(json)
                    .map_err(crate::fingerprint_json::FingerprintJsonError::from_rust)?;
                Self::from_rust(py, value)
            }
        });
    };
}
generator_parameters!(MorganGeneratorParams);
generator_parameters!(AtomPairGeneratorParams);
generator_parameters!(TopologicalTorsionGeneratorParams);

ck::project_configuration!(ck::AtomPairFingerprintParams);
ck::project_configuration!(ck::TopologicalTorsionFingerprintParams);
