//! Thin aromaticity-assignment projection over the detached core owner.

use cosmolkit_macros::mol_op_body;

use super::OperationError;
use crate::{AromaticityParams, DerivedState, PreservationProof, TopologyEditKind};

#[mol_op_body(with_assigned_aromaticity, parts)]
pub(crate) fn assign_aromaticity_impl(params: &AromaticityParams) -> Result<(), OperationError> {
    let topology = parts.checkout_topology()?;
    // Inspect the write-owned derived-cache working value for the source
    // ring guard below.
    let mut cache = parts.checkout_derived_cache()?;
    // Complete pinned source: Aromaticity.cpp setAromaticity ring guard.
    // RDKit✔️✔️:   VECT_INT_VECT srings;
    // RDKit✔️✔️:   if (mol.getRingInfo()->isInitialized()) {
    // RDKit✔️✔️:     srings = mol.getRingInfo()->atomRings();
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     MolOps::symmetrizeSSSR(mol, srings);
    // RDKit✔️✔️:   }
    // Behavior review: the pinned guard is INITIALIZED-only — existing
    // initialized rows (Fast/Sssr/Symm, including initialized-empty) are
    // BORROWED unchanged, never silently upgraded; ONLY absent state calls
    // the canonical symmetrized_sssr once, and those acquired rows are then
    // MOVED into the cache. Cost review: exactly zero finds on the borrowed
    // path and exactly one symmetrizeSSSR on the absent path; no row clone.
    let mut acquired: Option<cosmolkit_core::RingInfo> = None;
    let rings = match cache.valid_ring_info() {
        Some(existing) => existing,
        None => {
            let rings = match cosmolkit_core::symmetrized_sssr(
                &topology,
                &cosmolkit_core::RingSearchParams::default(),
            ) {
                Ok(rings) => rings,
                Err(error) => {
                    parts.install_topology(topology)?;
                    parts.install_derived_cache(cache)?;
                    return Err(OperationError::Rings(error));
                }
            };
            #[cfg(test)]
            crate::ops::cow_tests::ring_aromaticity_probe::record_acquisition();
            acquired.insert(rings)
        }
    };
    let assignment = match cosmolkit_core::assign_aromaticity(&topology, rings, params) {
        Ok(assignment) => assignment,
        Err(error) => {
            parts.install_topology(topology)?;
            parts.install_derived_cache(cache)?;
            return Err(OperationError::Aromaticity(error));
        }
    };

    if assignment.topology.atoms.len() != topology.atoms.len() {
        let actual = assignment.topology.atoms.len();
        let expected = topology.atoms.len();
        parts.install_topology(topology)?;
        parts.install_derived_cache(cache)?;
        return Err(OperationError::InvalidAlgorithmResult {
            operation: "with_assigned_aromaticity",
            field: "atom",
            actual,
            expected,
        });
    }
    if assignment.topology.bonds.len() != topology.bonds.len() {
        let actual = assignment.topology.bonds.len();
        let expected = topology.bonds.len();
        parts.install_topology(topology)?;
        parts.install_derived_cache(cache)?;
        return Err(OperationError::InvalidAlgorithmResult {
            operation: "with_assigned_aromaticity",
            field: "bond",
            actual,
            expected,
        });
    }
    if assignment.topology.adjacency != topology.adjacency {
        parts.install_topology(topology)?;
        parts.install_derived_cache(cache)?;
        return Err(OperationError::InvalidAlgorithmResult {
            operation: "with_assigned_aromaticity",
            field: "adjacency",
            actual: 1,
            expected: 0,
        });
    }
    if let Err(error) = assignment.topology.validate() {
        parts.install_topology(topology)?;
        parts.install_derived_cache(cache)?;
        return Err(OperationError::InvalidTopology(error));
    }

    parts.install_topology(assignment.topology)?;
    parts.record_topology_edit(TopologyEditKind::Local)?;
    // Borrowed rows stay exactly as supplied (validity explicitly
    // reaffirmed); acquired rows are MOVED into the cache and marked valid.
    if let Some(rings) = acquired {
        cache.install_ring_info(rings);
    }
    parts.install_derived_cache(cache)?;
    parts.mark_cache_updated(DerivedState::RINGS)?;
    parts.clear_cache(
        DerivedState::VALENCE
            .union(DerivedState::STEREO)
            .union(DerivedState::DRAWING)
            .union(DerivedState::FINGERPRINT),
    )?;
    parts.mark_cache_updated(DerivedState::AROMATICITY)?;
    parts.prove_preserved(
        DerivedState::RING_FAMILIES.union(DerivedState::COORDINATES),
        PreservationProof::AromaticityAssignment,
    )?;
    parts.apply_cip_policy()
}

/// C5: 12 real aromaticity operations (cyclohexane/benzene x absent/real
/// Fast/real Symm x value/in-place). The pinned initialized-only guard
/// borrows supplied rows with ZERO root acquisitions; only absent state
/// acquires Symm once.
#[cfg(all(
    test,
    feature = "cap-smiles",
    feature = "cap-aromaticity",
    feature = "cap-rings"
))]
mod ring_live_tests {
    use crate::{DerivedState, Molecule};

    fn raw(input: &str) -> Molecule {
        Molecule::from_smiles_with_params(
            input,
            &cosmolkit_smiles::SmilesParseParams {
                sanitize: false,
                remove_hydrogens: false,
                ..cosmolkit_smiles::SmilesParseParams::default()
            },
        )
        .unwrap()
    }

    fn supplied(input: &str, rings: Option<cosmolkit_core::RingInfo>) -> Molecule {
        let base = raw(input);
        crate::Molecule::from_smiles_parts_with_derived_state(
            base.topology().clone(),
            base.coordinate_block_runtime().clone(),
            base.properties().clone(),
            None,
            rings,
        )
        .unwrap()
    }

    #[test]
    fn ring_live_aromaticity_source_guard_move_and_peer_product() {
        let mut calls = 0usize;
        for input in ["C1CCCCC1", "c1ccccc1"] {
            let benzene = input.starts_with('c');
            let topology = raw(input).topology().clone();
            for (supply, state) in [
                ("absent", None),
                (
                    "fast",
                    Some(cosmolkit_core::fast_find_rings(&topology).unwrap()),
                ),
                (
                    "symm",
                    Some(
                        cosmolkit_core::symmetrized_sssr(
                            &topology,
                            &cosmolkit_core::RingSearchParams::default(),
                        )
                        .unwrap(),
                    ),
                ),
            ] {
                for in_place in [false, true] {
                    let label = format!("{input}/{supply}/ip={in_place}");
                    let baseline = state.clone();
                    let source = supplied(input, state.clone());
                    let observer = source.clone();
                    let original_cache = source.derived_cache_arc_runtime();
                    let acquisitions_before =
                        crate::ops::cow_tests::ring_aromaticity_probe::acquisitions();
                    let output = if in_place {
                        let mut target = source.clone();
                        target
                            .assign_aromaticity_with_params_(&Default::default())
                            .unwrap();
                        target
                    } else {
                        source
                            .with_assigned_aromaticity_with_params(&Default::default())
                            .unwrap()
                    };
                    calls += 1;
                    let acquisitions =
                        crate::ops::cow_tests::ring_aromaticity_probe::acquisitions()
                            - acquisitions_before;

                    // Root acquisition discipline: absent => exactly one
                    // Symm; supplied initialized rows => zero.
                    match supply {
                        "absent" => assert_eq!(acquisitions, 1, "{label}: one Symm"),
                        _ => assert_eq!(acquisitions, 0, "{label}: no acquisition"),
                    }

                    // Final cache state.
                    let cache = output.derived_cache_runtime();
                    assert!(
                        cache.valid_states().contains(DerivedState::RINGS),
                        "{label}"
                    );
                    assert!(
                        cache.valid_states().contains(DerivedState::AROMATICITY),
                        "{label}: aromaticity valid"
                    );
                    let rings = cache.valid_ring_info().unwrap();
                    assert!(rings.is_initialized(), "{label}");
                    match supply {
                        "absent" => {
                            assert_eq!(
                                rings.find_type(),
                                cosmolkit_core::RingFindType::SymmSssr,
                                "{label}: acquired Symm"
                            );
                        }
                        "fast" => {
                            assert_eq!(
                                rings.find_type(),
                                cosmolkit_core::RingFindType::Fast,
                                "{label}: Fast stays Fast"
                            );
                        }
                        _ => {
                            assert_eq!(
                                rings.find_type(),
                                cosmolkit_core::RingFindType::SymmSssr,
                                "{label}: Symm stays Symm"
                            );
                        }
                    }
                    // Ordered rows: both six-cycles have exactly one ring
                    // over all six atoms/bonds.
                    assert_eq!(rings.atom_rings().len(), 1, "{label}: rows");
                    let mut atoms_row: Vec<usize> =
                        rings.atom_rings()[0].iter().map(|a| a.index()).collect();
                    atoms_row.sort_unstable();
                    assert_eq!(atoms_row, vec![0, 1, 2, 3, 4, 5], "{label}: cycle");
                    // Supplied states keep their exact ordered value.
                    if supply != "absent" {
                        assert_eq!(Some(rings), baseline.as_ref(), "{label}: exact reuse");
                    }

                    // Topology outcome: benzene perceives aromatic, the
                    // cyclohexane stays non-aromatic; validity is enforced
                    // by the runtime install path.
                    assert_eq!(
                        output
                            .atoms()
                            .iter()
                            .filter(|atom| atom.is_aromatic())
                            .count(),
                        usize::from(benzene) * 6,
                        "{label}: aromatic atoms"
                    );
                    // Value semantics and peer safety.
                    assert_eq!(source, observer, "{label}: source value");
                    assert!(std::sync::Arc::ptr_eq(
                        &source.derived_cache_arc_runtime(),
                        &original_cache
                    ));
                }
            }
        }
        assert_eq!(calls, 12, "exact census");
    }
}
