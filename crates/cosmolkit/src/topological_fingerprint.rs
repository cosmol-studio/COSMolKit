//! Thin source RDKFingerprint facade borrowing the sole detached graph values.
use crate::{
    Fingerprint, Molecule, QueryGraph, TopologicalFingerprintError,
    TopologicalFingerprintOutputRequest, TopologicalFingerprintParams,
    TopologicalFingerprintResult,
};
impl Molecule {
    pub fn topological_fingerprint(&self) -> Result<Fingerprint, TopologicalFingerprintError> {
        self.topological_fingerprint_with_params(&TopologicalFingerprintParams::default())
    }
    pub fn topological_fingerprint_with_params(
        &self,
        params: &TopologicalFingerprintParams,
    ) -> Result<Fingerprint, TopologicalFingerprintError> {
        cosmolkit_fingerprints::topological_fingerprint(self.topology(), params)
    }
    pub fn topological_fingerprint_with_output(
        &self,
    ) -> Result<TopologicalFingerprintResult, TopologicalFingerprintError> {
        self.topological_fingerprint_with_output_with_params(
            &TopologicalFingerprintParams::default(),
            TopologicalFingerprintOutputRequest::default(),
        )
    }
    pub fn topological_fingerprint_with_output_with_params(
        &self,
        params: &TopologicalFingerprintParams,
        request: TopologicalFingerprintOutputRequest,
    ) -> Result<TopologicalFingerprintResult, TopologicalFingerprintError> {
        cosmolkit_fingerprints::topological_fingerprint_with_output(
            self.topology(),
            params,
            request,
        )
    }
}
/// Interpret a canonical query directly; no query-bearing Molecule is constructed.
pub fn topological_query_fingerprint_with_params(
    query: &QueryGraph,
    params: &TopologicalFingerprintParams,
) -> Result<Fingerprint, TopologicalFingerprintError> {
    cosmolkit_fingerprints::topological_query_fingerprint(query, params)
}
pub fn topological_query_fingerprint_with_output_with_params(
    query: &QueryGraph,
    params: &TopologicalFingerprintParams,
    request: TopologicalFingerprintOutputRequest,
) -> Result<TopologicalFingerprintResult, TopologicalFingerprintError> {
    cosmolkit_fingerprints::topological_query_fingerprint_with_output(query, params, request)
}
