//! Hydrogen operation bodies; chemistry remains in its final algorithm owner.

use cosmolkit_macros::mol_op_body;

use super::OperationError;
use crate::{AddHsParams, DerivedState, PreservationProof, TopologyEditKind};

#[mol_op_body(with_hydrogens, parts)]
pub(crate) fn add_hydrogens_impl(params: &AddHsParams) -> Result<(), OperationError> {
    let topology = parts.checkout_topology()?;
    let coordinates = parts.checkout_coordinates()?;
    let properties = parts.checkout_properties()?;
    let result =
        cosmolkit_core::add_hydrogens_with_params(topology, coordinates, properties, params)
            .map_err(OperationError::Hydrogen)?;

    parts.install_topology(result.topology)?;
    parts.install_coordinates(result.coordinates)?;
    parts.install_properties(result.properties)?;
    parts.record_topology_edit(TopologyEditKind::Appending)?;
    parts.record_topology_mapping(result.mapping)?;
    parts.apply_runtime_remap()?;
    parts.clear_cache(
        DerivedState::VALENCE
            .union(DerivedState::AROMATICITY)
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
                        // CK-VALENCE-001 stays independent of rings.
                        assert_eq!(
                            cache.valid_states().contains(DerivedState::VALENCE),
                            sanitize,
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

    // CK-VALENCE-001: sanitize=false allows an unsanitized molecule, not a stale
    // cache advertised as valid. Clear both the validity bit and stored value.
    // For sanitize=true, only a complete final-topology assignment may be
    // installed. None means unavailable, never a swallowed calculation error.
    // The existing operation_defined effect permits both update and clear;
    // it does not make historical RDKit cache values valid for current state.
    if !params.sanitize {
        parts.clear_cache(DerivedState::VALENCE)?;
    } else if let Some(valence) = result.final_valence {
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
