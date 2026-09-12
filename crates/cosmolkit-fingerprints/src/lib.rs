//! Detached molecular fingerprint boundaries.

use cosmolkit_model::TopologyBlock;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FingerprintError {
    Unsupported,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Fingerprint {
    bits: Vec<u32>,
}

impl Fingerprint {
    #[must_use]
    pub const fn empty() -> Self {
        Self { bits: Vec::new() }
    }

    #[must_use]
    pub fn on_bits(&self) -> &[u32] {
        &self.bits
    }
}

pub fn morgan(topology: &TopologyBlock) -> Result<Fingerprint, FingerprintError> {
    let _ = topology;
    Err(FingerprintError::Unsupported)
}

pub fn pattern(topology: &TopologyBlock) -> Result<Fingerprint, FingerprintError> {
    morgan(topology)
}
