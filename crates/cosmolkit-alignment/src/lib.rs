//! Detached molecular and conformer alignment boundaries.

use cosmolkit_model::{Conformer3D, TopologyBlock};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignmentError {
    Unsupported,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AlignmentResult {
    pub atom_pairs: Vec<(usize, usize)>,
    pub rmsd: f64,
}

pub fn align(
    reference: &TopologyBlock,
    probe: &TopologyBlock,
    reference_conformer: Option<&Conformer3D>,
    probe_conformer: Option<&Conformer3D>,
) -> Result<AlignmentResult, AlignmentError> {
    let _ = (reference, probe, reference_conformer, probe_conformer);
    Err(AlignmentError::Unsupported)
}
