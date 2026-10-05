//! Thin read-only MACCS API; FP owns keys and SEARCH owns matching.
use crate::{Fingerprint, MaccsFingerprintError, MaccsFingerprintParams, Molecule};
impl Molecule {
    pub fn maccs_fingerprint(&self) -> Result<Fingerprint, MaccsFingerprintError> {
        self.maccs_fingerprint_with_params(&MaccsFingerprintParams::default())
    }
    pub fn maccs_fingerprint_with_params(
        &self,
        params: &MaccsFingerprintParams,
    ) -> Result<Fingerprint, MaccsFingerprintError> {
        let cache = self.derived_cache_runtime();
        cosmolkit_fingerprints::maccs_fingerprint(
            self.topology(),
            cache.valid_ring_info(),
            cache.valence_assignment(),
            params,
        )
    }
    /// Source raw 167-bit vector with unused bit zero, before public projection.
    pub fn maccs_fingerprint_raw(&self) -> Result<Fingerprint, MaccsFingerprintError> {
        let cache = self.derived_cache_runtime();
        cosmolkit_fingerprints::maccs_fingerprint_raw(
            self.topology(),
            cache.valid_ring_info(),
            cache.valence_assignment(),
        )
    }
}
