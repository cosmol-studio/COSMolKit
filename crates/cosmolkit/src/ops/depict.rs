//! Thin 2D-coordinate operation projection over the detached depict owner.

use cosmolkit_macros::mol_op_body;

use super::OperationError;
use crate::{Coordinate2DParams, DerivedState, PreservationProof};

fn next_2d_conformer_id(
    coordinates: &crate::CoordinateBlock,
) -> Result<usize, crate::Coordinate2DError> {
    // RDKit✔️✔️: if (assignId) {
    // RDKit✔️✔️:   int maxId = -1;
    // RDKit✔️✔️:   for (auto cptr : d_confs) {
    // RDKit✔️✔️:     maxId = std::max((int)(cptr->getId()), maxId);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   maxId++;
    // RDKit✔️✔️:   conf->setId((unsigned int)maxId);
    // RDKit✔️✔️: }
    // Behavior review: for the modeled ordinary nonnegative identifier range,
    // the maximum existing dimension-local 2D identifier is incremented once.
    // CK's approved coordinate contract deliberately keeps the independent 3D
    // namespace out of this scan. Rust's wider representational boundary is
    // checked instead of reproducing the source's narrowing/overflow behavior.
    // Complexity review: one borrowed pass over the 2D table is O(n), creates
    // no collection or clone, and matches the source scan asymptotically.
    let Some(max_id) = coordinates
        .conformers_2d
        .iter()
        .map(|conformer| conformer.id())
        .max()
    else {
        return Ok(0);
    };
    max_id
        .checked_add(1)
        .ok_or(crate::Coordinate2DError::ConformerIdOverflow { max_id })
}

#[mol_op_body(with_2d_coordinates, parts)]
pub(crate) fn with_2d_coordinates_impl(params: &Coordinate2DParams) -> Result<(), OperationError> {
    let topology = parts.topology()?;
    let properties = parts.properties()?;
    let atom_count = topology.atoms.len();
    let conformer = cosmolkit_depict::compute_2d_coordinates(topology, properties, params)
        .map_err(OperationError::Coordinate2D)?;

    // RDKit❗✔️: unsigned int copyCoordinate(RDKit::ROMol &mol,
    // RDKit❗✔️:                             const std::list<EmbeddedFrag> &efrags,
    // RDKit❗✔️:                             bool clearConfs) {
    // RDKit❗✔️:   auto *conf = new RDKit::Conformer(mol.getNumAtoms());
    // RDKit❗✔️:   conf->set3D(false);
    // RDKit❗✔️:   if (clearConfs) { mol.clearConformers(); }
    // RDKit❗✔️:   return mol.addConformer(conf, true);
    // RDKit❗✔️: }
    // Behavior review: CK stores 2D and 3D conformers independently, so the
    // approved coordinate contract applies source clearConfs only to the 2D
    // table and preserves every independent 3D row. Final identifier assignment
    // is reproduced by `next_2d_conformer_id` at the runtime-owned install edge.
    // Complexity review: installation materializes only the declared coordinate
    // block; identifier selection is the separate linear scan reviewed above.
    let mut coordinates = parts.checkout_coordinates()?;
    let id = if params.clear_existing_2d {
        coordinates.clear_2d_conformers();
        0
    } else {
        next_2d_conformer_id(&coordinates).map_err(OperationError::Coordinate2D)?
    };
    coordinates
        .record_source_conformer_append(crate::CoordinateDimension::TwoD)
        .map_err(OperationError::InvalidCoordinates)?;
    coordinates.conformers_2d.push(conformer.with_id(id));
    coordinates
        .validate_for_atom_count(atom_count)
        .map_err(OperationError::InvalidCoordinates)?;
    parts.install_coordinates(coordinates)?;
    parts.clear_cache(DerivedState::STEREO.union(DerivedState::DRAWING))?;
    parts.prove_preserved(
        DerivedState::RINGS
            .union(DerivedState::RING_FAMILIES)
            .union(DerivedState::VALENCE)
            .union(DerivedState::AROMATICITY)
            .union(DerivedState::FINGERPRINT),
        PreservationProof::CoordinateOnly,
    )?;
    parts.apply_cip_policy()
}

#[cfg(all(test, feature = "cap-smiles"))]
#[path = "../../tests/fixtures/d2_prepared_property.rs"]
mod prepared_property_fixture;

#[cfg(all(test, feature = "cap-smiles"))]
mod prepared_property_tests {
    use super::*;
    use crate::{
        Atom, AtomId, AtomSpec, CoordinateBlock, Element, Molecule, MoleculeProperties,
        TopologyBlock,
    };
    use std::sync::Arc;

    fn identities(molecule: &Molecule) -> [usize; 4] {
        [
            Arc::as_ptr(&molecule.topology_arc_runtime()) as usize,
            Arc::as_ptr(&molecule.coordinates_arc_runtime()) as usize,
            Arc::as_ptr(&molecule.properties_arc_runtime()) as usize,
            Arc::as_ptr(&molecule.derived_cache_arc_runtime()) as usize,
        ]
    }

    fn coordinate_bits(coordinates: &CoordinateBlock) -> (Vec<Vec<[u64; 2]>>, Vec<Vec<[u64; 3]>>) {
        (
            coordinates
                .conformers_2d
                .iter()
                .map(|conf| {
                    conf.coordinates()
                        .iter()
                        .map(|xy| xy.map(f64::to_bits))
                        .collect()
                })
                .collect(),
            coordinates
                .conformers_3d
                .iter()
                .map(|conf| {
                    conf.coordinates()
                        .iter()
                        .map(|xyz| xyz.map(f64::to_bits))
                        .collect()
                })
                .collect(),
        )
    }

    #[test]
    fn d2_prepared_property_internal_twelve_source_peer_storage_calls() {
        let mut calls = 0;
        let mut discrepancies = Vec::new();
        for (line, smiles) in prepared_property_fixture::CASES {
            for orientation in [false, true] {
                for repeat in 0..2 {
                    let source = Molecule::from_smiles(smiles).unwrap();
                    assert!(source.property("_StereochemDone").is_some());
                    let peer = source.clone();
                    let ids = identities(&source);
                    assert_eq!(identities(&peer), ids);
                    let topology = source.topology().clone();
                    let coordinates = source.coordinate_block_runtime().clone();
                    let properties = source.properties().clone();
                    let cache = source.derived_cache_runtime().clone();
                    let before_bits = coordinate_bits(&coordinates);
                    let valence = cache
                        .valence_assignment()
                        .expect("actual sanitized constructor cache")
                        .clone();
                    let valence_address =
                        source.derived_cache_runtime().valence_assignment().unwrap() as *const _;
                    let params = Coordinate2DParams {
                        canonical_orientation: orientation,
                        ..Default::default()
                    };
                    let before_params = params.clone();
                    let result = source.with_2d_coordinates_with_params(&params);
                    calls += 1;
                    // Both complete four-block source/peer snapshots precede Result inspection.
                    for molecule in [&source, &peer] {
                        assert_eq!(identities(molecule), ids);
                        assert_eq!(molecule.topology(), &topology);
                        assert_eq!(molecule.coordinate_block_runtime(), &coordinates);
                        assert_eq!(
                            coordinate_bits(molecule.coordinate_block_runtime()),
                            before_bits
                        );
                        assert_eq!(molecule.properties(), &properties);
                        assert_eq!(molecule.derived_cache_runtime(), &cache);
                        assert_eq!(
                            molecule.derived_cache_runtime().valence_assignment(),
                            Some(&valence)
                        );
                        assert_eq!(
                            molecule
                                .derived_cache_runtime()
                                .valence_assignment()
                                .unwrap() as *const _,
                            valence_address
                        );
                    }
                    assert_eq!(params, before_params);
                    match result {
                        Err(error) => discrepancies
                            .push(format!("line:{line}/{orientation}/{repeat}: {error:?}")),
                        Ok(output) => {
                            let output_ids = identities(&output);
                            assert_eq!(output_ids[0], ids[0]);
                            assert_eq!(output_ids[2], ids[2]);
                            assert_ne!(output_ids[1], ids[1]);
                            assert_ne!(output_ids[3], ids[3]);
                            assert_eq!(output.topology(), &topology);
                            assert_eq!(output.properties(), &properties);
                            let stored = output.coordinate_block_runtime();
                            assert_eq!(stored.conformers_2d.len(), 1);
                            assert_eq!(stored.conformers_2d[0].id(), 0);
                            assert_eq!(stored.conformers_3d, coordinates.conformers_3d);
                            assert_eq!(
                                stored.source_coordinate_dim,
                                Some(crate::CoordinateDimension::TwoD)
                            );
                            assert_eq!(
                                output.derived_cache_runtime().valence_assignment(),
                                Some(&valence)
                            );
                            let allowed_invalidations =
                                DerivedState::STEREO.union(DerivedState::DRAWING);
                            assert_eq!(
                                output.derived_cache_runtime().valid_states(),
                                cache.valid_states().difference(allowed_invalidations)
                            );
                            discrepancies.extend(prepared_property_fixture::check(
                                output.topology(),
                                output.properties(),
                                stored.conformers_2d[0].coordinates(),
                                line,
                                orientation,
                            ));
                        }
                    }
                }
            }
        }
        println!(
            "D2-PREPARED internal_calls={calls} source_peer_preservation=12/12 discrepancies={}",
            discrepancies.len()
        );
        assert_eq!(calls, 12);
        assert!(discrepancies.is_empty(), "{}", discrepancies.join("\n"));
    }

    #[test]
    fn d2_prepared_property_two_typed_errors_preserve_inputs_and_peer() {
        // Invalid detached topology cannot enter a valid live Molecule. Test
        // that existing domain boundary directly, without bypassing runtime
        // construction; the second control is the real public map error.
        let mut invalid = TopologyBlock::default();
        invalid
            .atoms
            .push(Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)));
        let peer_topology = invalid.clone();
        let properties = MoleculeProperties::default();
        let before_properties = properties.clone();
        let params = Coordinate2DParams::default();
        let before_params = params.clone();
        let result = cosmolkit_depict::compute_2d_coordinates(&invalid, &properties, &params);
        assert_eq!(invalid, peer_topology);
        assert_eq!(properties, before_properties);
        assert_eq!(params, before_params);
        assert!(
            matches!(result, Err(crate::Coordinate2DError::InvalidTopology(
            cosmolkit_model::TopologyValidationError::AtomIdMismatch { position: 0, id })
        ) if id == AtomId::new(1))
        );

        let source = Molecule::from_smiles("CC").unwrap();
        let peer = source.clone();
        let ids = identities(&source);
        let topology = source.topology().clone();
        let coordinates = source.coordinate_block_runtime().clone();
        let properties = source.properties().clone();
        let cache = source.derived_cache_runtime().clone();
        let bits = coordinate_bits(&coordinates);
        let mut params = Coordinate2DParams::default();
        params.coordinate_map.insert(2, [-0.0, 1.0]);
        let params_before = params.clone();
        let map_bits = params.coordinate_map[&2].map(f64::to_bits);
        let result = source.with_2d_coordinates_with_params(&params);
        for molecule in [&source, &peer] {
            assert_eq!(identities(molecule), ids);
            assert_eq!(molecule.topology(), &topology);
            assert_eq!(molecule.coordinate_block_runtime(), &coordinates);
            assert_eq!(coordinate_bits(molecule.coordinate_block_runtime()), bits);
            assert_eq!(molecule.properties(), &properties);
            assert_eq!(molecule.derived_cache_runtime(), &cache);
        }
        assert_eq!(params, params_before);
        assert_eq!(params.coordinate_map[&2].map(f64::to_bits), map_bits);
        assert!(matches!(
            result,
            Err(OperationError::Coordinate2D(
                crate::Coordinate2DError::Fragment(
                    crate::Coordinate2DLayoutError::AtomIndexOutOfRange {
                        atom: 2,
                        atom_count: 2
                    }
                )
            ))
        ));
        println!("D2-PREPARED typed_error_controls=2 preservation=2/2");
    }
}
