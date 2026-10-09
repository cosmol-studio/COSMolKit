//! Thin Pattern API over detached concrete/query inputs.
use crate::{Fingerprint, Molecule, PatternFingerprintError, PatternFingerprintParams, QueryGraph};
impl Molecule {
    pub fn fingerprint_pattern(&self) -> Result<Fingerprint, PatternFingerprintError> {
        self.fingerprint_pattern_with_params(&PatternFingerprintParams::default())
    }
    pub fn fingerprint_pattern_with_params(
        &self,
        params: &PatternFingerprintParams,
    ) -> Result<Fingerprint, PatternFingerprintError> {
        cosmolkit_fingerprints::pattern_fingerprint(
            self.topology(),
            self.derived_cache_runtime().valid_ring_info(),
            params,
        )
    }
}
pub fn fingerprint_pattern_query(
    query: &QueryGraph,
) -> Result<Fingerprint, PatternFingerprintError> {
    fingerprint_pattern_query_with_params(query, &PatternFingerprintParams::default())
}
pub fn fingerprint_pattern_query_with_params(
    query: &QueryGraph,
    params: &PatternFingerprintParams,
) -> Result<Fingerprint, PatternFingerprintError> {
    cosmolkit_fingerprints::pattern_query_fingerprint(query, None, params)
}
