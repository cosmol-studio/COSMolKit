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
    // RDKit✔️✔️:     mol.updatePropertyCache(true);
    // RDKit✔️✔️:     mol.updatePropertyCache(false);
    // Both source branches leave actual atom scalar fields. Move the already
    // computed strict or nonstrict owner result; perform no second assignment.
    if let Some(valence) = assignment.final_valence.or(assignment.non_strict_valence) {
        let mut cache = parts.checkout_derived_cache()?;
        cache.install_valence_assignment(valence);
        parts.install_derived_cache(cache)?;
        parts.mark_cache_updated(DerivedState::VALENCE)?;
    } else {
        parts.clear_cache(DerivedState::VALENCE)?;
    }
    // Moved final ring state: the owner's final stage either supplies an
    // initialized result (stored and marked valid) or the state is
    // explicitly cleared. No local finder, no stale-row retention.
    match assignment.final_rings {
        Some(rings) => {
            let mut cache = parts.checkout_derived_cache()?;
            cache.install_ring_info(rings);
            parts.install_derived_cache(cache)?;
            parts.mark_cache_updated(DerivedState::RINGS)?;
        }
        None => parts.clear_cache(DerivedState::RINGS)?,
    }
    parts.clear_cache(
        DerivedState::RING_FAMILIES
            .union(DerivedState::AROMATICITY)
            .union(DerivedState::STEREO)
            .union(DerivedState::DRAWING)
            .union(DerivedState::FINGERPRINT),
    )?;
    parts.prove_preserved(
        DerivedState::COORDINATES,
        PreservationProof::SanitizeTopologyState,
    )?;
    parts.apply_cip_policy()?;
    // RDKit✔️✔️: mol.setProp(common_properties::numArom, narom, true);
    // The runtime's computed-property clearing must precede publication of
    // this already calculated source scalar. Use the declared property cap.
    if let Some(count) = assignment.aromatic_ring_count {
        let count = i32::try_from(count).map_err(|_| {
            OperationError::Sanitize(cosmolkit_core::SanitizeError::Aromaticity {
                stage: cosmolkit_core::SanitizeStage::SetAromaticity,
                source: cosmolkit_core::AromaticityError::IntegerOverflow {
                    field: "source numArom int",
                },
            })
        })?;
        let mut properties = parts.checkout_properties()?;
        let write = properties
            .set_computed_prop("numArom", count)
            .map_err(OperationError::InvalidProperty);
        parts.install_properties(properties)?;
        write?;
    }
    Ok(())
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

/// C2: 24 real sanitize operations (3 shapes x ALL/NONE/PROPERTIES/
/// SYMM_RINGS x raw/already-sanitized) proving moved final ring state,
/// explicit clearing with no stale rows, and a typed owner error that
/// preserves the source value including its installed rings.
#[cfg(all(
    test,
    feature = "cap-smiles",
    feature = "cap-sanitize",
    feature = "cap-rings"
))]
mod ring_live_tests {
    use super::*;
    use crate::{
        Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, CoordinateBlock, DerivedState,
        Element, MoleculeProperties, SanitizeOperations, TopologyBlock,
    };

    fn raw(input: &str) -> Molecule {
        Molecule::from_smiles_with_params(
            input,
            &cosmolkit_smiles::SmilesParseParams {
                sanitize: false,
                remove_hs: false,
                ..cosmolkit_smiles::SmilesParseParams::default()
            },
        )
        .unwrap()
    }

    #[test]
    fn ring_live_sanitize_state_product() {
        let mut calls = 0usize;
        for input in ["", "CC", "c1ccccc1"] {
            for (mask, operations) in [
                ("ALL", SanitizeOperations::ALL),
                ("NONE", SanitizeOperations::NONE),
                ("PROPERTIES", SanitizeOperations::PROPERTIES),
                ("SYMM_RINGS", SanitizeOperations::SYMM_RINGS),
            ] {
                for prepared in [false, true] {
                    let label = format!("{mask}/{input:?}/prepared={prepared}");
                    let source = if prepared {
                        raw(input)
                            .sanitize_with_params(&SanitizeParams {
                                operations: SanitizeOperations::ALL,
                            })
                            .unwrap()
                    } else {
                        raw(input)
                    };
                    let observer = source.clone();
                    let original_cache = source.derived_cache_arc_runtime();
                    let output = source
                        .sanitize_with_params(&SanitizeParams { operations })
                        .unwrap();
                    calls += 1;
                    let cache = output.derived_cache_runtime();

                    let symm = matches!(mask, "ALL" | "SYMM_RINGS");
                    assert_eq!(
                        cache.valid_states().contains(DerivedState::RINGS),
                        symm,
                        "{label}: RINGS validity"
                    );
                    assert_eq!(
                        cache.valid_states().contains(DerivedState::VALENCE),
                        true,
                        "{label}: VALENCE validity"
                    );
                    if symm {
                        let rings = cache.valid_ring_info().expect("{label}: installed");
                        assert!(rings.is_initialized(), "{label}: initialized");
                        assert_eq!(
                            rings.find_type(),
                            cosmolkit_core::RingFindType::SymmSssr,
                            "{label}: find type"
                        );
                        let (want_atoms, want_bonds) = match input {
                            "" => (0usize, 0usize),
                            "CC" => (2, 1),
                            _ => (6, 6),
                        };
                        assert_eq!(rings.atom_row_count(), want_atoms, "{label}: dims");
                        assert_eq!(rings.bond_row_count(), want_bonds, "{label}: dims");
                        if input == "c1ccccc1" {
                            assert_eq!(rings.atom_rings().len(), 1, "{label}: rows");
                            assert_eq!(rings.bond_rings().len(), 1, "{label}: bond rows");
                            let mut atoms_row: Vec<usize> = rings.atom_rings()[0]
                                .iter()
                                .map(|atom| atom.index())
                                .collect();
                            atoms_row.sort_unstable();
                            let mut bonds_row: Vec<usize> = rings.bond_rings()[0]
                                .iter()
                                .map(|bond| bond.index())
                                .collect();
                            bonds_row.sort_unstable();
                            assert_eq!(atoms_row, vec![0, 1, 2, 3, 4, 5], "{label}");
                            assert_eq!(bonds_row, vec![0, 1, 2, 3, 4, 5], "{label}");
                            for index in 0..6usize {
                                assert_eq!(
                                    rings.atom_members(AtomId::new(index)),
                                    &[0],
                                    "{label}: member {index}"
                                );
                                assert_eq!(
                                    rings.bond_members(BondId::new(index)),
                                    &[0],
                                    "{label}: bond member {index}"
                                );
                            }
                        } else {
                            assert!(rings.atom_rings().is_empty(), "{label}: rows");
                            assert!(rings.bond_rings().is_empty(), "{label}: bond rows");
                            for index in 0..want_atoms {
                                assert_eq!(
                                    rings.atom_members(AtomId::new(index)),
                                    &[] as &[usize],
                                    "{label}: member {index}"
                                );
                            }
                            for index in 0..want_bonds {
                                assert_eq!(
                                    rings.bond_members(BondId::new(index)),
                                    &[] as &[usize],
                                    "{label}: bond member {index}"
                                );
                            }
                        }
                    } else {
                        // Explicit clear: no stale old rows even when the
                        // already-sanitized source carried a valid Symm
                        // payload.
                        assert!(cache.valid_ring_info().is_none(), "{label}: absent");
                        assert!(cache.ring_info().is_none(), "{label}: storage cleared");
                    }
                    // Value semantics: the source and its cache are intact.
                    assert_eq!(source, observer, "{label}: source value");
                    assert!(std::sync::Arc::ptr_eq(
                        &source.derived_cache_arc_runtime(),
                        &original_cache
                    ));
                }
            }
        }
        assert_eq!(calls, 24, "exact census");
    }

    #[test]
    fn ring_live_sanitize_typed_error_preserves_source_ring_state() {
        // Real invalid valence (oxygen with three single bonds): the owner
        // fails in the PROPERTIES stage; the value form must preserve the
        // source molecule INCLUDING its installed Fast rings.
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
        .with_assigned_rings()
        .unwrap();
        let observer = source.clone();
        let original_cache = source.derived_cache_arc_runtime();
        let error = source
            .sanitize_with_params(&SanitizeParams {
                operations: SanitizeOperations::ALL,
            })
            .unwrap_err();
        assert!(
            matches!(
                &error,
                OperationError::Sanitize(SanitizeError::Properties { .. })
            ),
            "got {error:?}"
        );
        assert_eq!(source, observer);
        assert!(std::sync::Arc::ptr_eq(
            &source.derived_cache_arc_runtime(),
            &original_cache
        ));
        let cache = source.derived_cache_runtime();
        assert!(cache.valid_states().contains(DerivedState::RINGS));
        let rings = cache.valid_ring_info().unwrap();
        assert_eq!(rings.find_type(), cosmolkit_core::RingFindType::Fast);
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
                        if operations.contains(SanitizeOperations::SYMM_RINGS) {
                            // ALL runs the symmetrized ring stage: the final
                            // moved SymmSssr state (real rows on the aromatic
                            // cycle, initialized-empty on the acyclic form)
                            // is stored and marked valid.
                            assert_eq!(
                                cache.valid_states(),
                                DerivedState::VALENCE.union(DerivedState::RINGS)
                            );
                            let rings = cache.valid_ring_info().unwrap();
                            assert!(rings.is_initialized());
                            assert_eq!(rings.find_type(), cosmolkit_core::RingFindType::SymmSssr);
                            assert_eq!(rings.atom_rings().len(), usize::from(aromatic));
                        } else {
                            assert_eq!(cache.valid_states(), DerivedState::VALENCE);
                            assert!(cache.valid_ring_info().is_none());
                        }
                    } else {
                        let expected = cosmolkit_core::assign_valence_with_options_for_topology(
                            output.topology(),
                            cosmolkit_core::ValenceModel::RdkitLike,
                            false,
                        )
                        .unwrap();
                        assert_eq!(cache.valence_assignment(), Some(&expected));
                        // KEKULIZE-only: the ring-aware owner acquires SSSR
                        // exactly when aromatic kekulization work runs
                        // (K-RING frozen table: absent + mark=true acquires
                        // SSSR once); the acyclic form early-returns with
                        // the carrier absent. NONE never touches rings.
                        if operations == SanitizeOperations::KEKULIZE && aromatic {
                            assert_eq!(
                                cache.valid_states(),
                                DerivedState::VALENCE.union(DerivedState::RINGS),
                                "kekulize-only aromatic installs SSSR"
                            );
                            let rings = cache.valid_ring_info().unwrap();
                            assert!(rings.is_initialized());
                            assert_eq!(rings.find_type(), cosmolkit_core::RingFindType::Sssr);
                            assert_eq!(rings.atom_rings().len(), 1);
                        } else {
                            assert_eq!(cache.valid_states(), DerivedState::VALENCE);
                            assert!(cache.valid_ring_info().is_none());
                        }
                    }
                    assert_eq!(
                        output.property("source"),
                        Some(&crate::PropertyValue::String("retained".into()))
                    );
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
