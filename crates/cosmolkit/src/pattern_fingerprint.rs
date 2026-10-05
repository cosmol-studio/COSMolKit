//! Thin Pattern API over detached concrete/query inputs.
use crate::{Fingerprint, Molecule, PatternFingerprintError, PatternFingerprintParams, QueryGraph};
impl Molecule {
    pub fn pattern_fingerprint(&self) -> Result<Fingerprint, PatternFingerprintError> {
        self.pattern_fingerprint_with_params(&PatternFingerprintParams::default())
    }
    pub fn pattern_fingerprint_with_params(
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
pub fn pattern_query_fingerprint(
    query: &QueryGraph,
) -> Result<Fingerprint, PatternFingerprintError> {
    pattern_query_fingerprint_with_params(query, &PatternFingerprintParams::default())
}
pub fn pattern_query_fingerprint_with_params(
    query: &QueryGraph,
    params: &PatternFingerprintParams,
) -> Result<Fingerprint, PatternFingerprintError> {
    cosmolkit_fingerprints::pattern_query_fingerprint(query, None, params)
}
