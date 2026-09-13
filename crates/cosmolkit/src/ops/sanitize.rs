//! Thin sanitization and chemistry-problem projections over the detached core owner.

use cosmolkit_macros::mol_op_body;

use super::OperationError;
use crate::{
    ChemistryProblemReport, DerivedState, Molecule, PreservationProof, SanitizeError,
    SanitizeParams, TopologyEditKind,
};

#[mol_op_body(sanitize, parts)]
pub(crate) fn sanitize_impl(params: &SanitizeParams) -> Result<(), OperationError> {
    let topology = parts.checkout_topology()?;
    let assignment = match cosmolkit_core::sanitize_topology(&topology, params) {
        Ok(assignment) => assignment,
        Err(error) => {
            parts.install_topology(topology)?;
            return Err(OperationError::Sanitize(error));
        }
    };

    if assignment.topology.atoms.len() != topology.atoms.len() {
        let actual = assignment.topology.atoms.len();
        let expected = topology.atoms.len();
        parts.install_topology(topology)?;
        return Err(OperationError::InvalidAlgorithmResult {
            operation: "sanitize",
            field: "atom",
            actual,
            expected,
        });
    }
    if assignment.topology.bonds.len() != topology.bonds.len() {
        let actual = assignment.topology.bonds.len();
        let expected = topology.bonds.len();
        parts.install_topology(topology)?;
        return Err(OperationError::InvalidAlgorithmResult {
            operation: "sanitize",
            field: "bond",
            actual,
            expected,
        });
    }
    if let Some((candidate, source)) = assignment
        .topology
        .atoms
        .iter()
        .zip(&topology.atoms)
        .find(|(candidate, source)| candidate.id() != source.id())
    {
        let actual = candidate.id().index();
        let expected = source.id().index();
        parts.install_topology(topology)?;
        return Err(OperationError::InvalidAlgorithmResult {
            operation: "sanitize",
            field: "atom identity",
            actual,
            expected,
        });
    }
    if assignment
        .topology
        .bonds
        .iter()
        .zip(&topology.bonds)
        .any(|(candidate, source)| {
            candidate.id() != source.id()
                || candidate.begin() != source.begin()
                || candidate.end() != source.end()
        })
    {
        parts.install_topology(topology)?;
        return Err(OperationError::InvalidAlgorithmResult {
            operation: "sanitize",
            field: "bond identity or endpoints",
            actual: 1,
            expected: 0,
        });
    }
    if assignment.topology.adjacency != topology.adjacency {
        parts.install_topology(topology)?;
        return Err(OperationError::InvalidAlgorithmResult {
            operation: "sanitize",
            field: "adjacency",
            actual: 1,
            expected: 0,
        });
    }
    if assignment.topology.substance_groups != topology.substance_groups {
        parts.install_topology(topology)?;
        return Err(OperationError::InvalidAlgorithmResult {
            operation: "sanitize",
            field: "substance groups",
            actual: 1,
            expected: 0,
        });
    }
    if let Err(error) = assignment.topology.validate() {
        parts.install_topology(topology)?;
        return Err(OperationError::InvalidTopology(error));
    }

    parts.install_topology(assignment.topology)?;
    parts.record_topology_edit(TopologyEditKind::Local)?;
    parts.clear_cache(
        DerivedState::RINGS
            .union(DerivedState::RING_FAMILIES)
            .union(DerivedState::VALENCE)
            .union(DerivedState::AROMATICITY)
            .union(DerivedState::STEREO)
            .union(DerivedState::DRAWING)
            .union(DerivedState::FINGERPRINT),
    )?;
    parts.prove_preserved(
        DerivedState::COORDINATES,
        PreservationProof::SanitizeTopologyState,
    )?;
    parts.apply_cip_policy()
}

impl Molecule {
    /// Detect source-equivalent chemistry problems without mutating this molecule.
    pub fn detect_chemistry_problems(&self) -> Result<ChemistryProblemReport, SanitizeError> {
        self.detect_chemistry_problems_with_params(&SanitizeParams::default())
    }

    /// Detect source-equivalent chemistry problems using explicit stage selection.
    pub fn detect_chemistry_problems_with_params(
        &self,
        params: &SanitizeParams,
    ) -> Result<ChemistryProblemReport, SanitizeError> {
        cosmolkit_core::detect_chemistry_problems(self.topology(), params)
    }
}
