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
    // Only the owner's final PROPERTIES stage certifies valence for the
    // committed topology. Retain that result without calculating it again;
    // other stage selections must remove both stored values and validity.
    if let Some(valence) = assignment.final_valence {
        let mut cache = parts.checkout_derived_cache()?;
        cache.install_valence_assignment(valence);
        parts.install_derived_cache(cache)?;
        parts.mark_cache_updated(DerivedState::VALENCE)?;
    } else {
        parts.clear_cache(DerivedState::VALENCE)?;
    }
    parts.clear_cache(
        DerivedState::RINGS
            .union(DerivedState::RING_FAMILIES)
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

#[cfg(all(test, feature = "cap-valence"))]
mod tests {
    use super::*;
    use crate::{
        Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, CoordinateBlock, Element,
        MoleculeProperties, SanitizeOperations, TopologyBlock,
    };

    fn detached_source(aromatic: bool) -> Molecule {
        let atoms = if aromatic {
            (0..6)
                .map(|id| {
                    Atom::from_spec(
                        AtomId::new(id),
                        AtomSpec::new(Element::C).with_aromatic(true),
                    )
                })
                .collect()
        } else {
            [Element::C, Element::C, Element::O]
                .into_iter()
                .enumerate()
                .map(|(id, element)| Atom::from_spec(AtomId::new(id), AtomSpec::new(element)))
                .collect()
        };
        let bonds = if aromatic {
            (0..6)
                .map(|id| {
                    Bond::from_spec(
                        BondId::new(id),
                        BondSpec::new(
                            AtomId::new(id),
                            AtomId::new((id + 1) % 6),
                            BondOrder::Aromatic,
                        )
                        .with_aromatic(true),
                    )
                })
                .collect()
        } else {
            (0..2)
                .map(|id| {
                    Bond::from_spec(
                        BondId::new(id),
                        BondSpec::new(AtomId::new(id), AtomId::new(id + 1), BondOrder::Single),
                    )
                })
                .collect()
        };
        Molecule::from_parts(
            TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap(),
            CoordinateBlock::default(),
            MoleculeProperties::default()
                .with_prop("source", "retained")
                .unwrap(),
        )
        .unwrap()
    }

    #[test]
    fn sanitize_final_valence_runtime_stage_and_source_cache_product() {
        let mut calls = 0;
        for aromatic in [false, true] {
            for prepared in [false, true] {
                let cold = detached_source(aromatic);
                let source = if prepared {
                    cold.with_assigned_valence().unwrap()
                } else {
                    cold
                };
                let observer = source.clone();
                let original_cache = source.derived_cache_arc_runtime();
                let original_cache_value = original_cache.as_ref().clone();
                for operations in [
                    SanitizeOperations::ALL,
                    SanitizeOperations::PROPERTIES,
                    SanitizeOperations::NONE,
                    SanitizeOperations::KEKULIZE,
                ] {
                    let output = source
                        .sanitize_with_params(&SanitizeParams { operations })
                        .unwrap();
                    calls += 1;
                    let cache = output.derived_cache_runtime();
                    if operations.contains(SanitizeOperations::PROPERTIES) {
                        let expected = cosmolkit_core::assign_valence_with_options_for_topology(
                            output.topology(),
                            cosmolkit_core::ValenceModel::RdkitLike,
                            true,
                        )
                        .unwrap();
                        assert_eq!(cache.valence_assignment(), Some(&expected));
                        assert_eq!(cache.valid_states(), DerivedState::VALENCE);
                    } else {
                        assert_eq!(cache.valence_assignment(), None);
                        assert_eq!(cache.valid_states(), DerivedState::NONE);
                    }
                    assert_eq!(output.property("source"), Some("retained"));
                    assert_eq!(source, observer);
                    assert_eq!(source.derived_cache_runtime(), &original_cache_value);
                    assert!(std::sync::Arc::ptr_eq(
                        &source.derived_cache_arc_runtime(),
                        &original_cache,
                    ));
                    assert!(std::ptr::eq(
                        source.coordinate_block_runtime(),
                        output.coordinate_block_runtime(),
                    ));
                }
            }
        }
        assert_eq!(calls, 16);
    }

    #[test]
    fn sanitize_final_valence_failure_does_not_replace_source_topology_or_cache() {
        let source = Molecule::from_parts(
            TopologyBlock::try_from_parts(
                [Element::O, Element::C, Element::C, Element::C]
                    .into_iter()
                    .enumerate()
                    .map(|(id, element)| Atom::from_spec(AtomId::new(id), AtomSpec::new(element)))
                    .collect(),
                (1..4)
                    .map(|id| {
                        Bond::from_spec(
                            BondId::new(id - 1),
                            BondSpec::new(AtomId::new(0), AtomId::new(id), BondOrder::Single),
                        )
                    })
                    .collect(),
                Vec::new(),
                Vec::new(),
            )
            .unwrap(),
            CoordinateBlock::default(),
            MoleculeProperties::default(),
        )
        .unwrap()
        .with_assigned_valence_with_params(&cosmolkit_core::ValenceParams {
            strict: false,
            ..Default::default()
        })
        .unwrap();
        let observer = source.clone();
        let original_cache = source.derived_cache_arc_runtime();
        for operations in [SanitizeOperations::ALL, SanitizeOperations::PROPERTIES] {
            assert!(matches!(
                source.sanitize_with_params(&SanitizeParams { operations }),
                Err(OperationError::Sanitize(SanitizeError::Properties { .. })),
            ));
            assert_eq!(source, observer);
            assert!(std::sync::Arc::ptr_eq(
                &source.derived_cache_arc_runtime(),
                &original_cache,
            ));
            assert_eq!(source.derived_cache_runtime(), original_cache.as_ref());
        }
    }
}
