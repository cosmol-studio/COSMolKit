//! Thin Kekule-bond assignment projection over the detached core owner.

use cosmolkit_macros::mol_op_body;

use super::OperationError;
use crate::{DerivedState, KekulizeParams, PreservationProof, TopologyEditKind};

#[mol_op_body(with_kekulized_bonds, parts)]
pub(crate) fn kekulize_bonds_impl(params: &KekulizeParams) -> Result<(), OperationError> {
    let topology = parts.checkout_topology()?;
    // Inspect the write-owned derived-cache working value and BORROW its
    // valid ordinary rings into the existing ring-aware core Kekulize owner
    // with the existing query-state input (None here, exactly as the plain
    // entry). No local finder and no row clone: the borrow is read-only and
    // the owner moves its own result back.
    let mut cache = parts.checkout_derived_cache()?;
    let assignment = match cosmolkit_core::kekulize_with_query_state_and_ring_info(
        &topology,
        params,
        None,
        cache.valid_ring_info(),
    ) {
        Ok(assignment) => assignment,
        Err(error) => {
            parts.install_topology(topology)?;
            parts.install_derived_cache(cache)?;
            return Err(OperationError::Kekulize(error));
        }
    };

    if assignment.topology.atoms.len() != topology.atoms.len() {
        let actual = assignment.topology.atoms.len();
        let expected = topology.atoms.len();
        parts.install_topology(topology)?;
        parts.install_derived_cache(cache)?;
        return Err(OperationError::InvalidAlgorithmResult {
            operation: "with_kekulized_bonds",
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
            operation: "with_kekulized_bonds",
            field: "bond",
            actual,
            expected,
        });
    }
    if assignment.topology.adjacency != topology.adjacency {
        parts.install_topology(topology)?;
        parts.install_derived_cache(cache)?;
        return Err(OperationError::InvalidAlgorithmResult {
            operation: "with_kekulized_bonds",
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
    // ring_update None => the old ordinary state is unchanged, explicitly
    // handled (validity reaffirmed when previously valid, explicit clear
    // when absent); Some(initialized) => move replacement; Some(reset) =>
    // clear. No local row copy.
    let old_rings_valid = cache.valid_states().contains(DerivedState::RINGS);
    match assignment.ring_update {
        Some(update) if update.is_initialized() => {
            // Actual-site observation of the MOVED acquisition buffers
            // immediately before installation.
            #[cfg(test)]
            if !update.atom_rings().is_empty() {
                crate::ops::cow_tests::ring_kekulize_probe::record(&update);
            }
            cache.install_ring_info(update);
            parts.install_derived_cache(cache)?;
            parts.mark_cache_updated(DerivedState::RINGS)?;
        }
        Some(_) => {
            parts.install_derived_cache(cache)?;
            parts.clear_cache(DerivedState::RINGS)?;
        }
        None => {
            // Actual-site observation of the REUSED (checked-out) cache
            // buffers immediately before re-installation.
            #[cfg(test)]
            if let Some(rings) = cache.valid_ring_info() {
                if !rings.atom_rings().is_empty() {
                    crate::ops::cow_tests::ring_kekulize_probe::record(rings);
                }
            }
            parts.install_derived_cache(cache)?;
            if old_rings_valid {
                parts.mark_cache_updated(DerivedState::RINGS)?;
            } else {
                parts.clear_cache(DerivedState::RINGS)?;
            }
        }
    }
    parts.clear_cache(
        DerivedState::VALENCE
            .union(DerivedState::AROMATICITY)
            .union(DerivedState::STEREO)
            .union(DerivedState::DRAWING)
            .union(DerivedState::FINGERPRINT),
    )?;
    parts.prove_preserved(
        DerivedState::RING_FAMILIES.union(DerivedState::COORDINATES),
        PreservationProof::KekulizeBondAssignment,
    )?;
    parts.apply_cip_policy()
}

/// C4: 24 real whole-graph Kekulize calls on the hand-built aromatic
/// six-cycle with six supplied ring states x canonical x mark_atoms_bonds,
/// checked against the Step-30 receipt table frozen from the pinned K-RING
/// dispositions BEFORE execution.
#[cfg(all(
    test,
    feature = "cap-kekulize",
    feature = "cap-rings",
    feature = "cap-smiles"
))]
mod ring_live_tests {
    use crate::{
        Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, CoordinateBlock, DerivedState,
        Element, KekulizeParams, Molecule, MoleculeProperties, OperationError, TopologyBlock,
    };

    fn aromatic_cycle() -> TopologyBlock {
        TopologyBlock::try_from_parts(
            (0..6)
                .map(|id| {
                    Atom::from_spec(
                        AtomId::new(id),
                        AtomSpec::new(Element::C).with_aromatic(true),
                    )
                })
                .collect(),
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
                .collect(),
            Vec::new(),
            Vec::new(),
        )
        .unwrap()
    }

    fn supplied(topology: &TopologyBlock, rings: Option<cosmolkit_core::RingInfo>) -> Molecule {
        Molecule::from_parsed_parts_with_derived_state(
            topology.clone(),
            CoordinateBlock::default(),
            MoleculeProperties::default(),
            None,
            rings,
        )
        .unwrap()
    }

    fn states(topology: &TopologyBlock) -> Vec<(&'static str, Option<cosmolkit_core::RingInfo>)> {
        vec![
            ("none", None),
            (
                "other-empty",
                Some(cosmolkit_core::RingInfo::new(
                    cosmolkit_core::RingFindType::OtherOrUnknown,
                    6,
                    6,
                )),
            ),
            (
                "fast-empty",
                Some(cosmolkit_core::RingInfo::new(
                    cosmolkit_core::RingFindType::Fast,
                    6,
                    6,
                )),
            ),
            (
                "fast",
                Some(cosmolkit_core::fast_find_rings(topology).unwrap()),
            ),
            (
                "sssr",
                Some(
                    cosmolkit_core::find_sssr(
                        topology,
                        &cosmolkit_core::RingSearchParams::default(),
                    )
                    .unwrap(),
                ),
            ),
            (
                "symm",
                Some(
                    cosmolkit_core::symmetrized_sssr(
                        topology,
                        &cosmolkit_core::RingSearchParams::default(),
                    )
                    .unwrap(),
                ),
            ),
        ]
    }

    #[test]
    fn ring_live_kekulize_transport_error_reuse_reset_product() {
        let mut calls = 0usize;
        let topology = aromatic_cycle();
        for (supply, state) in states(&topology) {
            for canonical in [false, true] {
                for mark in [false, true] {
                    let label = format!("{supply}/c={canonical}/m={mark}");
                    let baseline = state.clone();
                    let source = supplied(&topology, state.clone());
                    let observer = source.clone();
                    let original_cache = source.derived_cache_arc_runtime();
                    let params = KekulizeParams {
                        mark_atoms_bonds: mark,
                        canonical,
                        max_backtracks: KekulizeParams::default().max_backtracks,
                    };
                    let probe_before = crate::ops::cow_tests::ring_kekulize_probe::len();
                    let result = source.with_kekulized_bonds_with_params(&params);
                    calls += 1;

                    // Expected outcome per the corrected Step-30 table.
                    let expect_error =
                        mark && (supply == "fast-empty" || (supply == "other-empty" && !canonical));
                    let expect_replaced_sssr = !expect_error
                        && (supply == "none" || (supply == "other-empty" && canonical));
                    if expect_error {
                        let error = result
                            .err()
                            .unwrap_or_else(|| panic!("{label}: expected AromaticAtomOutsideRing"));
                        assert!(
                            matches!(
                                &error,
                                OperationError::Kekulize(
                                    cosmolkit_core::KekulizeError::AromaticAtomOutsideRing { atom }
                                ) if *atom == AtomId::new(0)
                            ),
                            "{label}: got {error:?}"
                        );
                        // Value failure preserves the source molecule and
                        // its supplied ring state.
                        assert_eq!(source, observer, "{label}: source value");
                        assert!(std::sync::Arc::ptr_eq(
                            &source.derived_cache_arc_runtime(),
                            &original_cache
                        ));
                        assert_eq!(
                            source.derived_cache_runtime().ring_info(),
                            baseline.as_ref()
                        );
                        continue;
                    }
                    let output = result.unwrap_or_else(|error| panic!("{label}: {error:?}"));
                    let cache = output.derived_cache_runtime();
                    assert!(
                        cache.valid_states().contains(DerivedState::RINGS),
                        "{label}"
                    );
                    let rings = cache.valid_ring_info().unwrap();
                    if expect_replaced_sssr {
                        assert!(rings.is_initialized(), "{label}");
                        assert_eq!(
                            rings.find_type(),
                            cosmolkit_core::RingFindType::Sssr,
                            "{label}: replaced with SSSR"
                        );
                        assert_eq!(rings.atom_rings().len(), 1, "{label}: rows");
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
                        }
                    } else {
                        // ring_update None: the supplied state is unchanged.
                        assert_eq!(Some(rings), baseline.as_ref(), "{label}: reused");
                    }
                    // Actual-site buffer proof: on nonempty reuse/move the
                    // buffers observed inside the transaction are EXACTLY
                    // the buffers committed in the live output cache (the
                    // caller's shared cache was COW-detached at checkout, so
                    // pre-call supply pointers are not expected to survive).
                    let observations =
                        crate::ops::cow_tests::ring_kekulize_probe::after(probe_before);
                    if !rings.atom_rings().is_empty() {
                        assert_eq!(observations.len(), 1, "{label}: one observation");
                        assert_eq!(
                            rings.atom_rings().as_ptr() as usize,
                            observations[0].0,
                            "{label}: atom buffers preserved through commit"
                        );
                        assert_eq!(
                            rings.bond_rings().as_ptr() as usize,
                            observations[0].1,
                            "{label}: bond buffers preserved through commit"
                        );
                    } else {
                        assert!(observations.is_empty(), "{label}: no observation");
                    }
                    // Topology outcome: trusted-empty borrowed rows yield
                    // zero ring candidates, so NO kekulization work runs and
                    // the whole cycle stays aromatic; real-row cells kekulize
                    // to alternating single/double bonds, and the source
                    // markAtomsBonds flag only controls clearing aromatic
                    // flags (mark=false retains ALL six, mark=true clears).
                    let borrowed_empty = !expect_error
                        && ((supply == "other-empty" && !canonical) || supply == "fast-empty");
                    if borrowed_empty {
                        assert!(
                            output
                                .bonds()
                                .iter()
                                .all(|bond| bond.order() == BondOrder::Aromatic),
                            "{label}: no work on trusted empty rows"
                        );
                        assert_eq!(
                            output
                                .atoms()
                                .iter()
                                .filter(|atom| atom.is_aromatic())
                                .count(),
                            6,
                            "{label}: atoms stay aromatic"
                        );
                    } else {
                        let singles = output
                            .bonds()
                            .iter()
                            .filter(|bond| bond.order() == BondOrder::Single)
                            .count();
                        let doubles = output
                            .bonds()
                            .iter()
                            .filter(|bond| bond.order() == BondOrder::Double)
                            .count();
                        assert_eq!((singles, doubles), (3, 3), "{label}: kekulized orders");
                        assert_eq!(
                            output
                                .atoms()
                                .iter()
                                .filter(|atom| atom.is_aromatic())
                                .count(),
                            if mark { 0 } else { 6 },
                            "{label}: atom flags follow mark"
                        );
                    }
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
}
