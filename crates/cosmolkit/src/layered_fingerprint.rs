//! Thin Layered facade over detached concrete/query inputs.
use crate::{
    Fingerprint, LayeredFingerprintError, LayeredFingerprintParams, LayeredFingerprintResult,
    Molecule, QueryGraph,
};
impl Molecule {
    pub fn layered_fingerprint(&self) -> Result<Fingerprint, LayeredFingerprintError> {
        self.layered_fingerprint_with_params(&LayeredFingerprintParams::default())
    }
    pub fn layered_fingerprint_with_params(
        &self,
        params: &LayeredFingerprintParams,
    ) -> Result<Fingerprint, LayeredFingerprintError> {
        cosmolkit_fingerprints::layered_fingerprint(
            self.topology(),
            self.derived_cache_runtime().valid_ring_info(),
            params,
        )
    }
    pub fn layered_fingerprint_with_output(
        &self,
    ) -> Result<LayeredFingerprintResult, LayeredFingerprintError> {
        self.layered_fingerprint_with_output_with_params(&LayeredFingerprintParams::default())
    }
    pub fn layered_fingerprint_with_output_with_params(
        &self,
        params: &LayeredFingerprintParams,
    ) -> Result<LayeredFingerprintResult, LayeredFingerprintError> {
        cosmolkit_fingerprints::layered_fingerprint_with_output(
            self.topology(),
            self.derived_cache_runtime().valid_ring_info(),
            params,
        )
    }
}

/// Interpret the canonical query graph without a concrete-molecule conversion.
pub fn layered_query_fingerprint_with_params(
    query: &QueryGraph,
    params: &LayeredFingerprintParams,
) -> Result<Fingerprint, LayeredFingerprintError> {
    cosmolkit_fingerprints::layered_query_fingerprint(query, None, params)
}

pub fn layered_query_fingerprint_with_output_with_params(
    query: &QueryGraph,
    params: &LayeredFingerprintParams,
) -> Result<LayeredFingerprintResult, LayeredFingerprintError> {
    cosmolkit_fingerprints::layered_query_fingerprint_with_output(query, None, params)
}
