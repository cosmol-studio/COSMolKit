//! Hydrogen operation bodies; chemistry remains in its final algorithm owner.

use cosmolkit_macros::mol_op_body;

use super::OperationError;
use crate::{AddHsParams, DerivedState, PreservationProof, TopologyEditKind};

#[mol_op_body(with_hydrogens, parts)]
pub(crate) fn add_hydrogens_impl(params: &AddHsParams) -> Result<(), OperationError> {
    let topology = parts.checkout_topology()?;
    let coordinates = parts.checkout_coordinates()?;
    let properties = parts.checkout_properties()?;
    let mut cache = parts.checkout_derived_cache()?;
    let source_valence = cache.valence_assignment().cloned();
    let result = cosmolkit_core::add_hydrogens_with_source_valence(
        topology,
        coordinates,
        properties,
        params,
        source_valence,
    )
    .map_err(OperationError::Hydrogen)?;

    // ROOT CK-474bdce: publish actual selected-parent/new-H source scalars,
    // retaining untouched rows through the canonical validated cache path.
    cache.install_valence_assignment(result.final_valence);
    parts.install_derived_cache(cache)?;
    parts.mark_cache_updated(DerivedState::VALENCE)?;
    parts.install_topology(result.topology)?;
    parts.install_coordinates(result.coordinates)?;
    parts.install_properties(result.properties)?;
    parts.record_topology_edit(TopologyEditKind::Appending)?;
    parts.record_topology_mapping(result.mapping)?;
    parts.apply_runtime_remap()?;
    parts.clear_cache(
        DerivedState::AROMATICITY
            .union(DerivedState::STEREO)
            .union(DerivedState::DRAWING)
            .union(DerivedState::FINGERPRINT),
    )?;
    parts.prove_preserved(
        DerivedState::RINGS.union(DerivedState::RING_FAMILIES),
        PreservationProof::LeafAtomAppend,
    )?;
    parts.apply_cip_policy()
}

/// C3: 24 real RemoveHs operations (3 shapes x remove_nonimplicit x
/// sanitize x value/in-place). Only a NONEMPTY graph with BOTH flags true
/// installs the moved final SymmSssr rows (the source SIZE guard runs the
/// nested sanitize even when no hydrogen was removed); every other cell
/// explicitly clears ordinary rings.
#[cfg(all(
    test,
    feature = "cap-smiles",
    feature = "cap-hydrogens",
    feature = "cap-rings"
))]
mod ring_live_tests {
    use crate::{AtomId, BondId, DerivedState, Molecule, MoleculeProperties};

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

    fn params(remove_nonimplicit: bool, sanitize: bool) -> cosmolkit_core::RemoveHsParams {
        cosmolkit_core::RemoveHsParams {
            remove_nonimplicit,
            sanitize,
            update_explicit_count: false,
            ..cosmolkit_core::RemoveHsParams::default()
        }
    }

    #[test]
    fn ring_live_remove_hs_state_id_and_size_guard_product() {
        let mut calls = 0usize;
        for input in ["", "CC", "[H]C1CCCCC1"] {
            for remove_nonimplicit in [false, true] {
                for sanitize in [false, true] {
                    for in_place in [false, true] {
                        let label =
                            format!("{input:?}/ni={remove_nonimplicit}/s={sanitize}/ip={in_place}");
                        let source = raw(input);
                        let observer = source.clone();
                        let original_cache = source.derived_cache_arc_runtime();
                        let options = params(remove_nonimplicit, sanitize);
                        let output = if in_place {
                            let mut target = source.clone();
                            target.remove_hydrogens_with_params_(&options).unwrap();
                            target
                        } else {
                            source.without_hydrogens_with_params(&options).unwrap()
                        };
                        calls += 1;

                        // Literal final IDs and counts.
                        let hydrogen_cycle = input == "[H]C1CCCCC1";
                        let removed = hydrogen_cycle && remove_nonimplicit;
                        let want_atoms: usize = match input {
                            "" => 0,
                            "CC" => 2,
                            _ => usize::from(!removed) + 6,
                        };
                        let want_bonds: usize = match input {
                            "" => 0,
                            "CC" => 1,
                            _ => usize::from(!removed) + 6,
                        };
                        assert_eq!(output.num_atoms(), want_atoms, "{label}: atoms");
                        assert_eq!(output.num_bonds(), want_bonds, "{label}: bonds");
                        if hydrogen_cycle {
                            let carbon_start = usize::from(!removed);
                            for index in carbon_start..want_atoms {
                                assert_eq!(
                                    output.atoms()[index].element(),
                                    crate::Element::C,
                                    "{label}: element {index}"
                                );
                            }
                        }

                        // Only nonempty && remove_nonimplicit && sanitize
                        // installs the moved final SymmSssr state.
                        let symm = !input.is_empty() && remove_nonimplicit && sanitize;
                        let cache = output.derived_cache_runtime();
                        assert_eq!(
                            cache.valid_states().contains(DerivedState::RINGS),
                            symm,
                            "{label}: RINGS validity"
                        );
                        if symm {
                            let rings = cache.valid_ring_info().unwrap();
                            assert!(rings.is_initialized(), "{label}: initialized");
                            assert_eq!(
                                rings.find_type(),
                                cosmolkit_core::RingFindType::SymmSssr,
                                "{label}: find type"
                            );
                            assert_eq!(rings.atom_row_count(), want_atoms, "{label}: dims");
                            assert_eq!(rings.bond_row_count(), want_bonds, "{label}: dims");
                            if hydrogen_cycle {
                                // H0 removed and final IDs compacted to 0..5.
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
                                    assert_eq!(
                                        rings.bond_members(BondId::new(index)),
                                        &[0],
                                        "{label}: bond member {index}"
                                    );
                                }
                            } else {
                                // The SIZE guard regenerates even without any
                                // hydrogen removal: CC ends initialized-empty.
                                assert!(rings.atom_rings().is_empty(), "{label}: rows");
                                assert!(rings.bond_rings().is_empty(), "{label}: bond rows");
                                for index in 0..want_atoms {
                                    assert_eq!(
                                        rings.atom_members(AtomId::new(index)),
                                        &[] as &[usize],
                                        "{label}: member {index}"
                                    );
                                }
                            }
                        } else {
                            assert!(cache.valid_ring_info().is_none(), "{label}: absent");
                            assert!(cache.ring_info().is_none(), "{label}: storage cleared");
                        }
                        // Source scalar rows survive compaction, independently of rings.
                        assert_eq!(
                            cache.valid_states().contains(DerivedState::VALENCE),
                            true,
                            "{label}: VALENCE validity"
                        );
                        // Value semantics and value/in-place agreement.
                        assert_eq!(source, observer, "{label}: source value");
                        assert!(std::sync::Arc::ptr_eq(
                            &source.derived_cache_arc_runtime(),
                            &original_cache
                        ));
                        let value_output = source.without_hydrogens_with_params(&options);
                        let mut in_place_target = source.clone();
                        let in_place_result = in_place_target
                            .remove_hydrogens_with_params_(&options)
                            .map(|()| in_place_target);
                        match (&value_output, &in_place_result) {
                            (Ok(value), Ok(place)) => {
                                assert_eq!(value, place, "{label}: value == in-place");
                                assert_eq!(
                                    value.derived_cache_runtime(),
                                    place.derived_cache_runtime(),
                                    "{label}: cache equality"
                                );
                            }
                            _ => panic!("{label}: repeat calls diverged"),
                        }
                    }
                }
            }
        }
        assert_eq!(calls, 24, "exact census");
    }

    #[test]
    fn ring_live_remove_hs_isotope_retention_and_final_routing() {
        // D retention with the full guard: the isotope hydrogen survives the
        // non-implicit pass, and ONLY the final-pass routed Symm state is
        // installed (the preliminary isotope pass never publishes rows).
        for in_place in [false, true] {
            let label = format!("D/ip={in_place}");
            let source = raw("[2H]C1CCCCC1");
            let options = params(true, true);
            let output = if in_place {
                let mut target = source.clone();
                target.remove_hydrogens_with_params_(&options).unwrap();
                target
            } else {
                source.without_hydrogens_with_params(&options).unwrap()
            };
            assert_eq!(output.num_atoms(), 7, "{label}: D retained");
            assert_eq!(output.atoms()[0].element(), crate::Element::H, "{label}");
            assert_eq!(output.atoms()[0].isotope(), Some(2), "{label}: isotope");
            let cache = output.derived_cache_runtime();
            let rings = cache.valid_ring_info().expect("{label}: installed");
            assert_eq!(
                rings.find_type(),
                cosmolkit_core::RingFindType::SymmSssr,
                "{label}"
            );
            assert_eq!(rings.atom_rings().len(), 1, "{label}: rows");
            let mut atoms_row: Vec<usize> = rings.atom_rings()[0]
                .iter()
                .map(|atom| atom.index())
                .collect();
            atoms_row.sort_unstable();
            assert_eq!(atoms_row, vec![1, 2, 3, 4, 5, 6], "{label}: cycle atoms");
            assert_eq!(
                rings.atom_members(AtomId::new(0)),
                &[] as &[usize],
                "{label}: D membership empty"
            );
            for index in 1..7usize {
                assert_eq!(
                    rings.atom_members(AtomId::new(index)),
                    &[0],
                    "{label}: member {index}"
                );
            }
            let _ = MoleculeProperties::default();
        }
    }
}

#[mol_op_body(without_hydrogens, parts)]
pub(crate) fn remove_hydrogens_impl(
    params: &cosmolkit_core::RemoveHsParams,
) -> Result<(), OperationError> {
    let topology = parts.checkout_topology()?;
    let coordinates = parts.checkout_coordinates()?;
    let properties = parts.checkout_properties()?;
    let result =
        cosmolkit_core::remove_hydrogens_with_params(topology, coordinates, properties, params)
            .map_err(OperationError::Hydrogen)?;

    parts.install_topology(result.topology)?;
    parts.install_coordinates(result.coordinates)?;
    parts.install_properties(result.properties)?;
    parts.record_topology_edit(TopologyEditKind::Compacting)?;
    parts.record_topology_mapping(result.mapping)?;
    parts.apply_runtime_remap()?;

    // BEGIN RDKIT CPP FUNCTION MolOps::removeHs source scalar publication
    // RDKit✔️✔️:   for (auto atom : mol.atoms()) {
    // RDKit✔️✔️:     atom->updatePropertyCache(false);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   mol.clearComputedProps(true);
    // END RDKIT CPP FUNCTION MolOps::removeHs source scalar publication
    // ROOT CK-4c0846 withdraws CK-VALENCE-001's intentional deviation.
    // The owner's final, mapped scalar rows are source-observable materialized
    // state even when they differ from a fresh topology calculation. Install
    // them through the existing operation_defined effect and ordinary runtime
    // coherence checks. This moves the carrier; no extra chemistry or clone.
    if let Some(valence) = result.final_valence {
        let mut cache = parts.checkout_derived_cache()?;
        cache.install_valence_assignment(valence);
        parts.install_derived_cache(cache)?;
        parts.mark_cache_updated(DerivedState::VALENCE)?;
    } else {
        parts.clear_cache(DerivedState::VALENCE)?;
    }
    // Moved final ring state from the FINAL pass only — never a preliminary
    // isotope-pass value and never preserved old rows. The owner's source
    // SIZE guard (old_atom_count != 0 && remove_nonimplicit && sanitize)
    // decides whether the nested sanitize produced final rows; Some is
    // stored and marked valid, None is an explicit clear. No local finder.
    match result.final_rings {
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
    parts.apply_cip_policy()
}

#[cfg(all(
    test,
    feature = "cap-smiles",
    feature = "cap-hydrogens",
    feature = "cap-valence",
    feature = "cap-rings"
))]
mod add_source_scalar_live_tests {
    use crate::{
        AddHsParams, Atom, AtomId, AtomSpec, CoordinateBlock, DerivedState, Element, Molecule,
        MoleculeProperties, OperationError, TopologyBlock,
    };

    #[test]
    fn prepared_source_add_hs_publishes_source_scalars_and_preserves_peer_cow() {
        for in_place in [false, true] {
            let source = Molecule::from_smiles("CCO").unwrap();
            let peer = source.clone();
            let old_cache = source.derived_cache_arc_runtime();
            let old_assignment = old_cache.valence_assignment().unwrap().clone();
            let out = if in_place {
                let mut out = source.clone();
                out.add_hydrogens_().unwrap();
                out
            } else {
                source.with_hydrogens().unwrap()
            };
            assert_eq!(out.num_atoms(), 9);
            assert_eq!(out.num_bonds(), 8);
            let cache = out.derived_cache_runtime();
            assert!(cache.valid_states().contains(DerivedState::VALENCE));
            let v = cache.valence_assignment().unwrap();
            assert_eq!(v.explicit_valence, [4, 4, 2, 1, 1, 1, 1, 1, 1]);
            assert_eq!(v.implicit_hydrogens, [0; 9]);
            assert_eq!(peer.num_atoms(), 3);
            assert_eq!(source.num_atoms(), 3);
            assert_eq!(
                peer.derived_cache_runtime().valence_assignment(),
                Some(&old_assignment)
            );
            assert!(std::sync::Arc::ptr_eq(
                &old_cache,
                &source.derived_cache_arc_runtime()
            ));
            assert!(std::sync::Arc::ptr_eq(
                &old_cache,
                &peer.derived_cache_arc_runtime()
            ));
            assert_eq!(out.atom_metadata(false).unwrap().len(), 9);
        }
    }

    #[test]
    fn remove_then_add_reads_actual_stale_source_counts_before_refreshing_selected_parents() {
        for in_place in [false, true] {
            let expanded = Molecule::from_smiles("CCO")
                .unwrap()
                .with_hydrogens()
                .unwrap();
            let source = expanded
                .without_hydrogens_with_params(&cosmolkit_core::RemoveHsParams {
                    sanitize: false,
                    ..Default::default()
                })
                .unwrap();
            let peer = source.clone();
            let before = source.derived_cache_arc_runtime();
            assert_eq!(source.num_atoms(), 3);
            assert_eq!(
                before.valence_assignment().unwrap().explicit_valence,
                [4, 4, 2]
            );
            assert_eq!(
                before.valence_assignment().unwrap().implicit_hydrogens,
                [0, 0, 0]
            );
            // Source snapshots the stale zero counts before parent refresh.
            // A fresh-topology preparation would incorrectly append six Hs here.
            let first = if in_place {
                let mut out = source.clone();
                out.add_hydrogens_().unwrap();
                out
            } else {
                source.with_hydrogens().unwrap()
            };
            assert_eq!(first.num_atoms(), 3);
            assert_eq!(first.num_bonds(), 2);
            let valence = first.derived_cache_runtime().valence_assignment().unwrap();
            assert_eq!(valence.explicit_valence, [1, 2, 1]);
            assert_eq!(valence.implicit_hydrogens, [3, 2, 1]);
            // The NEXT source call uses the actually refreshed counts.
            let second = first.with_hydrogens().unwrap();
            assert_eq!(second.num_atoms(), 9);
            assert_eq!(second.num_bonds(), 8);
            assert_eq!(
                second
                    .derived_cache_runtime()
                    .valence_assignment()
                    .unwrap()
                    .implicit_hydrogens,
                [0; 9]
            );
            assert_eq!(source, peer);
            assert!(std::sync::Arc::ptr_eq(
                &before,
                &source.derived_cache_arc_runtime()
            ));
            assert!(std::sync::Arc::ptr_eq(
                &before,
                &peer.derived_cache_arc_runtime()
            ));
            assert_eq!(
                source
                    .derived_cache_runtime()
                    .valence_assignment()
                    .unwrap()
                    .explicit_valence,
                [4, 4, 2]
            );
        }
    }

    #[test]
    fn repeated_source_add_hs_retains_prepared_state_without_appending_again() {
        let source = Molecule::from_smiles("CCO")
            .unwrap()
            .with_hydrogens()
            .unwrap();
        let peer = source.clone();
        let out = source.with_hydrogens().unwrap();
        assert_eq!(out.num_atoms(), 9);
        assert_eq!(out.num_bonds(), 8);
        assert_eq!(
            out.derived_cache_runtime().valence_assignment(),
            source.derived_cache_runtime().valence_assignment()
        );
        assert_eq!(source, peer);
    }

    #[test]
    fn unprepared_implicit_source_is_a_structured_error_even_explicit_only() {
        let source = Molecule::from_parts(
            TopologyBlock::try_from_parts(
                vec![Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C))],
                vec![],
                vec![],
                vec![],
            )
            .unwrap(),
            CoordinateBlock::default(),
            MoleculeProperties::default(),
        )
        .unwrap();
        let peer = source.clone();
        for explicit in [false, true] {
            assert!(
                matches!(source.with_hydrogens_with_params(&AddHsParams { explicit_only:explicit,..Default::default() }),Err(OperationError::Hydrogen(cosmolkit_core::HydrogenError::Valence(cosmolkit_core::ValenceError::ImplicitValenceCacheNotInitialized { atom })) ) if atom==AtomId::new(0))
            );
            assert_eq!(source, peer);
            assert!(
                source
                    .derived_cache_runtime()
                    .valence_assignment()
                    .is_none()
            );
        }
    }

    #[test]
    fn leaf_append_preserves_exact_source_ring_prefix_and_new_h_members_are_empty() {
        let source = Molecule::from_smiles("C1CC1").unwrap();
        let peer = source.clone();
        let original = source
            .derived_cache_runtime()
            .valid_ring_info()
            .unwrap()
            .clone();
        let out = source.with_hydrogens().unwrap();
        assert_eq!(out.num_atoms(), 9);
        let rings = out.derived_cache_runtime().valid_ring_info().unwrap();
        assert_eq!(rings, &original);
        assert_eq!(rings.atom_row_count(), 3);
        assert_eq!(rings.bond_row_count(), 3);
        for i in 3..out.num_atoms() {
            assert_eq!(rings.atom_members(AtomId::new(i)), &[] as &[usize]);
            assert_eq!(rings.num_atom_rings(AtomId::new(i)), 0);
        }
        assert_eq!(source, peer);
    }
}
