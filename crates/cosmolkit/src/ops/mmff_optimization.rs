//! Runtime commit adapter for source-backed detached MMFF optimization.
use crate::{DerivedState, Molecule, OperationError, PreservationProof, TopologyEditKind};
pub use cosmolkit_forcefields::{
    MmffConformerOptimizationParams, MmffOptimizationParams, MmffOptimizeMoleculeConfResult,
};
use cosmolkit_macros::{MoleculeResult, mol_op_body};
use std::{error::Error, fmt, sync::Arc};
#[derive(Clone, Debug, PartialEq, MoleculeResult)]
pub struct MmffOptimizeMoleculeResult<M = Molecule> {
    #[pending_molecule]
    pub molecule: M,
    pub needs_more: i32,
}
impl<M> MmffOptimizeMoleculeResult<M> {
    pub const fn needs_more(&self) -> bool {
        self.needs_more > 0
    }
    pub const fn status_code(&self) -> i32 {
        self.needs_more
    }
}
#[derive(Clone, Debug, PartialEq, MoleculeResult)]
pub struct MmffOptimizeMoleculeConfsResult<M = Molecule> {
    #[pending_molecule]
    pub molecule: M,
    pub conformer_results: Vec<MmffOptimizeMoleculeConfResult>,
}
/// Shared identity retains concrete owner causes when OperationError is cloned.
#[derive(Clone, Debug)]
pub struct MmffOptimizationError(Arc<cosmolkit_forcefields::MmffOptimizationError>);
impl PartialEq for MmffOptimizationError {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }
}
impl fmt::Display for MmffOptimizationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(&*self.0, f)
    }
}
impl Error for MmffOptimizationError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        Some(&*self.0)
    }
}
pub(crate) fn owner_error(cause: cosmolkit_forcefields::MmffOptimizationError) -> OperationError {
    OperationError::MmffOptimization(MmffOptimizationError(Arc::new(cause)))
}

#[mol_op_body(with_mmff_optimized, parts)]
pub(crate) fn with_mmff_optimized_impl(
    params: &MmffOptimizationParams,
) -> Result<
    MmffOptimizeMoleculeResult<crate::PendingMolecule<super::WithMmffOptimizedAccess>>,
    OperationError,
> {
    // Source wrappers delegate properties/builder/minimization to the owner.
    // This adapter acquires write capabilities only for the declared blocks.
    // Every successful checkout is restored before a fallible return.
    let mut topology = parts.checkout_topology()?;
    let mut coordinates = match parts.checkout_coordinates() {
        Ok(value) => value,
        Err(error) => {
            parts.install_topology(topology)?;
            return Err(error);
        }
    };
    let mut properties = match parts.checkout_properties() {
        Ok(value) => value,
        Err(error) => {
            parts.install_topology(topology)?;
            parts.install_coordinates(coordinates)?;
            return Err(error);
        }
    };
    let mut cache = match parts.checkout_derived_cache() {
        Ok(value) => value,
        Err(error) => {
            parts.install_topology(topology)?;
            parts.install_coordinates(coordinates)?;
            parts.install_properties(properties)?;
            return Err(error);
        }
    };
    let outcome = (|| {
        let outcome = cosmolkit_forcefields::optimize_mmff_single(
            &mut topology,
            &mut coordinates,
            &properties,
            params,
            cache.valid_ring_info(),
        )
        .map_err(owner_error)?;
        // RDKit❗✔️:   if (!mol.hasProp(common_properties::_MMFFSanitized)) {
        // RDKit❗✔️:     mol.setProp(common_properties::_MMFFSanitized, 1, true);
        // RDKit❗✔️:   }
        // AtomTyper.cpp constructor owns this exact presence guard. Rust
        // molecule properties retain the established textual representation.
        // This constant-key write allocates one property value, never chemistry.
        if properties.prop("_MMFFSanitized").is_none() {
            properties
                .set_computed_prop("_MMFFSanitized", 1_i32)
                .map_err(OperationError::InvalidProperty)?;
        }
        Ok::<_, OperationError>(outcome)
    })();
    parts.install_topology(topology)?;
    parts.install_coordinates(coordinates)?;
    parts.install_properties(properties)?;
    let outcome = match outcome {
        Ok(mut outcome) => {
            if let Some(rings) = outcome.acquired_rings.take() {
                cache.install_ring_info(rings);
            }
            parts.install_derived_cache(cache)?;
            outcome
        }
        Err(error) => {
            parts.install_derived_cache(cache)?;
            return Err(error);
        }
    };
    parts.record_topology_edit(TopologyEditKind::Local)?;
    parts.mark_cache_updated(DerivedState::RINGS.union(DerivedState::AROMATICITY))?;
    parts.clear_cache(
        DerivedState::VALENCE
            .union(DerivedState::STEREO)
            .union(DerivedState::DRAWING)
            .union(DerivedState::FINGERPRINT),
    )?;
    parts.prove_preserved(
        DerivedState::RING_FAMILIES.union(DerivedState::COORDINATES),
        PreservationProof::MmffPreparedOptimization,
    )?;
    parts.apply_cip_policy()?;
    Ok(MmffOptimizeMoleculeResult {
        molecule: parts.pending_molecule()?,
        needs_more: outcome.status,
    })
}

#[mol_op_body(with_mmff_optimized_confs, parts)]
pub(crate) fn with_mmff_optimized_confs_impl(
    params: &MmffConformerOptimizationParams,
) -> Result<
    MmffOptimizeMoleculeConfsResult<crate::PendingMolecule<super::WithMmffOptimizedConfsAccess>>,
    OperationError,
> {
    // Source wrappers delegate properties/builder/minimization to the owner.
    // This adapter acquires write capabilities only for the declared blocks.
    // Every successful checkout is restored before a fallible return.
    let mut topology = parts.checkout_topology()?;
    let mut coordinates = match parts.checkout_coordinates() {
        Ok(value) => value,
        Err(error) => {
            parts.install_topology(topology)?;
            return Err(error);
        }
    };
    let mut properties = match parts.checkout_properties() {
        Ok(value) => value,
        Err(error) => {
            parts.install_topology(topology)?;
            parts.install_coordinates(coordinates)?;
            return Err(error);
        }
    };
    let mut cache = match parts.checkout_derived_cache() {
        Ok(value) => value,
        Err(error) => {
            parts.install_topology(topology)?;
            parts.install_coordinates(coordinates)?;
            parts.install_properties(properties)?;
            return Err(error);
        }
    };
    let outcome = (|| {
        let outcome = cosmolkit_forcefields::optimize_mmff_conformers(
            &mut topology,
            &mut coordinates,
            &properties,
            params,
            cache.valid_ring_info(),
        )
        .map_err(owner_error)?;
        // RDKit❗✔️:   if (!mol.hasProp(common_properties::_MMFFSanitized)) {
        // RDKit❗✔️:     mol.setProp(common_properties::_MMFFSanitized, 1, true);
        // RDKit❗✔️:   }
        // AtomTyper.cpp constructor owns this exact presence guard. Rust
        // molecule properties retain the established textual representation.
        // This constant-key write allocates one property value, never chemistry.
        if properties.prop("_MMFFSanitized").is_none() {
            properties
                .set_computed_prop("_MMFFSanitized", 1_i32)
                .map_err(OperationError::InvalidProperty)?;
        }
        Ok::<_, OperationError>(outcome)
    })();
    parts.install_topology(topology)?;
    parts.install_coordinates(coordinates)?;
    parts.install_properties(properties)?;
    let outcome = match outcome {
        Ok(mut outcome) => {
            if let Some(rings) = outcome.acquired_rings.take() {
                cache.install_ring_info(rings);
            }
            parts.install_derived_cache(cache)?;
            outcome
        }
        Err(error) => {
            parts.install_derived_cache(cache)?;
            return Err(error);
        }
    };
    parts.record_topology_edit(TopologyEditKind::Local)?;
    parts.mark_cache_updated(DerivedState::RINGS.union(DerivedState::AROMATICITY))?;
    parts.clear_cache(
        DerivedState::VALENCE
            .union(DerivedState::STEREO)
            .union(DerivedState::DRAWING)
            .union(DerivedState::FINGERPRINT),
    )?;
    parts.prove_preserved(
        DerivedState::RING_FAMILIES.union(DerivedState::COORDINATES),
        PreservationProof::MmffPreparedOptimization,
    )?;
    parts.apply_cip_policy()?;
    Ok(MmffOptimizeMoleculeConfsResult {
        molecule: parts.pending_molecule()?,
        conformer_results: outcome.conformer_results,
    })
}

#[cfg(all(test, feature = "cap-rings", feature = "cap-valence"))]
mod original_public_optimization_tests {
    use super::*;
    use crate::{
        AtomSpec, BondOrder, BondSpec, Conformer3D, CoordinateBlock, Element, MoleculeBuilder,
    };
    fn bonded_pair_with_3d_conformer(
        first: AtomSpec,
        second: AtomSpec,
        order: BondOrder,
        coords: Vec<[f64; 3]>,
    ) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(first);
        let a1 = builder.add_atom(second);
        builder
            .add_bond(BondSpec::new(a0, a1, order))
            .expect("test bond");
        builder.add_3d_conformer(coords).expect("test 3d conformer");
        builder.build().expect("test bonded pair with conformer")
    }

    fn single_atom_with_3d_conformer(atom: AtomSpec, coord: [f64; 3]) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(atom);
        builder
            .add_3d_conformer(vec![coord])
            .expect("test 3d conformer");
        builder.build().expect("test single atom with conformer")
    }

    fn bonded_pair_with_named_3d_conformers(
        first: AtomSpec,
        second: AtomSpec,
        order: BondOrder,
        first_coords: Vec<[f64; 3]>,
        second_coords: Vec<[f64; 3]>,
    ) -> Molecule {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(first);
        let a1 = builder.add_atom(second);
        builder
            .add_bond(BondSpec::new(a0, a1, order))
            .expect("test bond");
        let base = builder.build().unwrap();
        Molecule::from_parts(
            base.topology().clone(),
            CoordinateBlock {
                conformers_3d: vec![
                    Conformer3D::new(0, first_coords, true),
                    Conformer3D::new(7, second_coords, true),
                ],
                ..CoordinateBlock::default()
            },
            base.properties().clone(),
        )
        .unwrap()
    }
    fn mmff_optimize_molecule(
        input: &Molecule,
        variant: &str,
        iters: i32,
        threshold: f64,
        id: isize,
        ignore: bool,
    ) -> Result<MmffOptimizeMoleculeResult, OperationError> {
        input.with_mmff_optimized_with_params(&MmffOptimizationParams {
            mmff_variant: variant.into(),
            max_iterations: iters,
            non_bonded_threshold: threshold,
            conformer_id: if id == -1 {
                None
            } else {
                Some(
                    usize::try_from(id)
                        .expect("signed negative selector is a language boundary case"),
                )
            },
            ignore_interfragment_interactions: ignore,
        })
    }
    fn mmff_optimize_molecule_confs(
        input: &Molecule,
        threads: i32,
        iters: i32,
        variant: &str,
        threshold: f64,
        ignore: bool,
    ) -> Result<MmffOptimizeMoleculeConfsResult, OperationError> {
        input.with_mmff_optimized_confs_with_params(&MmffConformerOptimizationParams {
            num_threads: threads,
            max_iterations: iters,
            mmff_variant: variant.into(),
            non_bonded_threshold: threshold,
            ignore_interfragment_interactions: ignore,
        })
    }
    #[test]
    fn mmff_public_api_mmff_optimize_molecule_returns_value_style_result_for_typed_molecule() {
        let molecule = bonded_pair_with_3d_conformer(
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            BondOrder::Single,
            vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
        );
        let original_coords = molecule.conformers_3d()[0].coordinates().to_vec();

        let result = mmff_optimize_molecule(&molecule, "MMFF94", 25, 100.0, -1, true)
            .expect("typed MMFF molecule should optimize");

        assert_eq!(result.needs_more, 0);
        assert_ne!(
            result.molecule.conformers_3d()[0].coordinates(),
            original_coords
        );
        assert_eq!(molecule.conformers_3d()[0].coordinates(), original_coords);
    }

    #[test]
    fn mmff_public_api_mmff_optimize_molecule_returns_minus_one_for_missing_atom_params() {
        let molecule = single_atom_with_3d_conformer(AtomSpec::new(Element::HE), [0.0, 0.0, 0.0]);
        let original_coords = molecule.conformers_3d()[0].coordinates().to_vec();

        let result = mmff_optimize_molecule(&molecule, "MMFF94", 25, 100.0, -1, true)
            .expect("missing MMFF atom typing should map to wrapper -1 result");

        assert_eq!(result.needs_more, -1);
        assert_eq!(
            result.molecule.conformers_3d()[0].coordinates(),
            original_coords
        );
        assert_eq!(molecule.conformers_3d()[0].coordinates(), original_coords);
    }

    #[test]
    fn mmff_public_api_mmff_optimize_molecule_uses_rdkit_variant_parser_for_invalid_variant() {
        let molecule = bonded_pair_with_3d_conformer(
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            BondOrder::Single,
            vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
        );

        let reference = mmff_optimize_molecule(&molecule, "MMFF94", 25, 100.0, -1, true)
            .expect("reference MMFF94 optimize should run");
        let result = mmff_optimize_molecule(&molecule, "MMFF94S", 25, 100.0, -1, true)
            .expect("invalid uppercase MMFF variant should fall back like RDKit parser");

        assert_eq!(result.needs_more, reference.needs_more);
        assert_eq!(
            result.molecule.conformers_3d()[0].coordinates(),
            reference.molecule.conformers_3d()[0].coordinates()
        );
    }

    #[test]
    fn mmff_public_api_mmff_optimize_molecule_reports_max_iteration_limit() {
        let molecule = bonded_pair_with_3d_conformer(
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            BondOrder::Single,
            vec![[0.0, 0.0, 0.0], [2.5, 0.0, 0.0]],
        );
        let original_coords = molecule.conformers_3d()[0].coordinates().to_vec();

        let result = mmff_optimize_molecule(&molecule, "MMFF94", 0, 100.0, -1, true)
            .expect("empty-typed MMFF optimize should run");

        assert_eq!(result.needs_more, 1);
        assert_eq!(
            result.molecule.conformers_3d()[0].coordinates(),
            original_coords
        );
    }

    #[test]
    fn mmff_public_api_mmff_optimize_molecule_updates_selected_named_conformer_only() {
        let molecule = bonded_pair_with_named_3d_conformers(
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            BondOrder::Single,
            vec![[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]],
            vec![[0.0, 0.0, 0.0], [2.5, 0.0, 0.0]],
        );
        let first_coords = molecule.conformers_3d()[0].coordinates().to_vec();
        let selected_coords = molecule.conformers_3d()[1].coordinates().to_vec();

        let result = mmff_optimize_molecule(&molecule, "MMFF94", 25, 100.0, 7, true)
            .expect("MMFF optimize should preserve unselected conformers");

        assert_eq!(result.needs_more, 0);
        assert_eq!(
            result.molecule.conformers_3d()[0].coordinates(),
            first_coords
        );
        assert_ne!(
            result.molecule.conformers_3d()[1].coordinates(),
            selected_coords
        );
        assert_eq!(molecule.conformers_3d()[1].coordinates(), selected_coords);
    }

    #[test]
    fn shared_forcefield_conformer_driver_mmff_returns_value_style_results_for_all_conformers() {
        let molecule = bonded_pair_with_named_3d_conformers(
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            BondOrder::Single,
            vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            vec![[0.0, 0.0, 0.0], [2.3, 0.0, 0.0]],
        );
        let original_first = molecule.conformers_3d()[0].coordinates().to_vec();
        let original_second = molecule.conformers_3d()[1].coordinates().to_vec();

        let result = mmff_optimize_molecule_confs(&molecule, 1, 25, "MMFF94", 100.0, true)
            .expect("missing MMFF atom typing should return wrapper-style -1 results");

        assert_eq!(result.conformer_results.len(), 2);
        assert!(
            result
                .conformer_results
                .iter()
                .all(|entry| entry.needs_more == 0 && entry.energy.is_finite())
        );
        assert_ne!(
            result.molecule.conformers_3d()[0].coordinates(),
            original_first
        );
        assert_ne!(
            result.molecule.conformers_3d()[1].coordinates(),
            original_second
        );
        assert_eq!(molecule.conformers_3d()[0].coordinates(), original_first);
        assert_eq!(molecule.conformers_3d()[1].coordinates(), original_second);
    }

    #[test]
    fn shared_forcefield_conformer_driver_mmff_returns_minus_one_for_missing_atom_params() {
        let molecule = single_atom_with_3d_conformer(AtomSpec::new(Element::HE), [0.0, 0.0, 0.0]);
        let original_coords = molecule.conformers_3d()[0].coordinates().to_vec();

        let result = mmff_optimize_molecule_confs(&molecule, 1, 25, "MMFF94", 100.0, true)
            .expect("missing MMFF atom typing should map to wrapper -1 result");

        assert_eq!(result.conformer_results.len(), 1);
        assert_eq!(result.conformer_results[0].needs_more, -1);
        assert_eq!(result.conformer_results[0].energy, -1.0);
        assert_eq!(
            result.molecule.conformers_3d()[0].coordinates(),
            original_coords
        );
        assert_eq!(molecule.conformers_3d()[0].coordinates(), original_coords);
    }

    #[test]
    fn mmff_public_api_mmff_optimize_molecule_confs_uses_rdkit_variant_parser_for_invalid_variant()
    {
        let molecule = bonded_pair_with_named_3d_conformers(
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            BondOrder::Single,
            vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            vec![[0.0, 0.0, 0.0], [2.2, 0.0, 0.0]],
        );

        let result = mmff_optimize_molecule_confs(&molecule, 1, 25, "MMFF94S", 100.0, true)
            .expect("invalid uppercase MMFF variant should fall back like RDKit parser");

        let reference = mmff_optimize_molecule_confs(&molecule, 1, 25, "MMFF94", 100.0, true)
            .expect("reference MMFF94 conformer optimize should run");

        assert_eq!(result.conformer_results.len(), 2);
        assert_eq!(result.conformer_results, reference.conformer_results);
        assert_eq!(
            result.molecule.conformers_3d()[0].coordinates(),
            reference.molecule.conformers_3d()[0].coordinates()
        );
        assert_eq!(
            result.molecule.conformers_3d()[1].coordinates(),
            reference.molecule.conformers_3d()[1].coordinates()
        );
    }

    #[test]
    fn shared_forcefield_conformer_driver_mmff_handles_non_positive_thread_request_like_non_threaded_rdkit_build()
     {
        let molecule = bonded_pair_with_named_3d_conformers(
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            BondOrder::Single,
            vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            vec![[0.0, 0.0, 0.0], [2.2, 0.0, 0.0]],
        );

        let zero_threads = mmff_optimize_molecule_confs(&molecule, 0, 5, "MMFF94", 100.0, true)
            .expect("zero-thread request should use non-threaded RDKit path");
        let negative_threads =
            mmff_optimize_molecule_confs(&molecule, -1, 5, "MMFF94", 100.0, true)
                .expect("negative-thread request should use non-threaded RDKit path");

        assert_eq!(zero_threads.conformer_results.len(), 2);
        assert_eq!(negative_threads.conformer_results.len(), 2);
        for (left, right) in zero_threads
            .conformer_results
            .iter()
            .zip(negative_threads.conformer_results.iter())
        {
            assert_eq!(left.needs_more, right.needs_more);
            assert!(left.energy.is_finite());
            assert!(right.energy.is_finite());
        }
    }

    #[test]
    fn shared_forcefield_conformer_driver_mmff_reports_max_iteration_limit() {
        let molecule = bonded_pair_with_3d_conformer(
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            BondOrder::Single,
            vec![[0.0, 0.0, 0.0], [2.5, 0.0, 0.0]],
        );
        let original_coords = molecule.conformers_3d()[0].coordinates().to_vec();

        let result = mmff_optimize_molecule_confs(&molecule, 1, 0, "MMFF94", 100.0, true)
            .expect("current modeled MMFF wrapper path should return -1 before minimization");

        assert_eq!(result.conformer_results.len(), 1);
        assert_eq!(result.conformer_results[0].needs_more, 1);
        assert!(result.conformer_results[0].energy.is_finite());
        assert_eq!(
            result.molecule.conformers_3d()[0].coordinates(),
            original_coords
        );
    }

    #[test]
    fn shared_forcefield_conformer_driver_mmff_preserves_all_named_conformers() {
        let molecule = bonded_pair_with_named_3d_conformers(
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            BondOrder::Single,
            vec![[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]],
            vec![[0.0, 0.0, 0.0], [2.5, 0.0, 0.0]],
        );
        let first_coords = molecule.conformers_3d()[0].coordinates().to_vec();
        let second_coords = molecule.conformers_3d()[1].coordinates().to_vec();

        let result = mmff_optimize_molecule_confs(&molecule, 1, 25, "MMFF94", 100.0, true)
            .expect("MMFF conformer optimization should preserve current modeled coordinates");

        assert_eq!(result.conformer_results.len(), 2);
        assert_ne!(
            result.molecule.conformers_3d()[0].coordinates(),
            first_coords
        );
        assert_ne!(
            result.molecule.conformers_3d()[1].coordinates(),
            second_coords
        );
        assert_eq!(molecule.conformers_3d()[0].coordinates(), first_coords);
        assert_eq!(molecule.conformers_3d()[1].coordinates(), second_coords);
    }

    #[test]
    fn mmff_live_missing_conformer_preserves_all_source_and_peer_blocks() {
        let source = bonded_pair_with_3d_conformer(
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            BondOrder::Single,
            vec![[0., 0., 0.], [2., 0., 0.]],
        );
        let peer = source.clone();
        let original = source.clone();
        let error = source
            .with_mmff_optimized_with_params(&MmffOptimizationParams {
                conformer_id: Some(999),
                ..Default::default()
            })
            .unwrap_err();
        assert!(matches!(&error, OperationError::MmffOptimization(_)));
        assert!(
            Error::source(&error)
                .unwrap()
                .source()
                .unwrap()
                .source()
                .is_some()
        );
        assert_eq!(error, error.clone());
        assert_eq!(source, original);
        assert_eq!(peer, original);
        for candidate in [&source, &peer] {
            assert!(Arc::ptr_eq(
                &original.topology_arc_runtime(),
                &candidate.topology_arc_runtime()
            ));
            assert!(Arc::ptr_eq(
                &original.coordinates_arc_runtime(),
                &candidate.coordinates_arc_runtime()
            ));
            assert!(Arc::ptr_eq(
                &original.properties_arc_runtime(),
                &candidate.properties_arc_runtime()
            ));
            assert!(Arc::ptr_eq(
                &original.derived_cache_arc_runtime(),
                &candidate.derived_cache_arc_runtime()
            ));
        }
    }
    #[test]
    fn mmff_live_source_property_guard_ring_rows_and_unchanged_block_sharing() {
        let base = bonded_pair_with_3d_conformer(
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            BondOrder::Single,
            vec![[0., 0., 0.], [2., 0., 0.]],
        )
        .with_assigned_rings()
        .unwrap()
        .with_assigned_valence()
        .unwrap();
        for existing in [false, true] {
            let mut props = base.properties().clone();
            props.set_prop("unchanged", "metadata").unwrap();
            if existing {
                props.set_prop("_MMFFSanitized", "already").unwrap();
            }
            let source = Molecule::from_runtime_parts(
                base.topology_arc_runtime(),
                base.coordinates_arc_runtime(),
                Arc::new(props),
                base.derived_cache_arc_runtime(),
            )
            .unwrap();
            let peer = source.clone();
            let rings = source.derived_cache_runtime().ring_info().cloned();
            let result = source
                .with_mmff_optimized_with_params(&MmffOptimizationParams {
                    max_iterations: 0,
                    ..Default::default()
                })
                .unwrap();
            assert_eq!(result.needs_more, 1);
            assert_eq!(source, peer);
            assert_eq!(
                result.molecule.properties().prop("_MMFFSanitized"),
                Some(&if existing {
                    cosmolkit_model::PropertyValue::String("already".into())
                } else {
                    cosmolkit_model::PropertyValue::Int(1)
                })
            );
            assert_eq!(
                result
                    .molecule
                    .properties()
                    .is_prop_computed("_MMFFSanitized")
                    .unwrap(),
                !existing
            );
            assert_eq!(
                result.molecule.properties().prop("numArom"),
                source.properties().prop("numArom")
            );
            assert_eq!(
                result.molecule.properties().prop("unchanged"),
                Some(&cosmolkit_model::PropertyValue::String("metadata".into()))
            );
            assert_eq!(
                result.molecule.derived_cache_runtime().ring_info(),
                rings.as_ref()
            );
            assert!(
                result
                    .molecule
                    .derived_cache_runtime()
                    .valid_states()
                    .contains(DerivedState::RINGS.union(DerivedState::AROMATICITY))
            );
            assert!(
                !result
                    .molecule
                    .derived_cache_runtime()
                    .valid_states()
                    .contains(DerivedState::VALENCE)
            );
            assert!(Arc::ptr_eq(
                &source.topology_arc_runtime(),
                &result.molecule.topology_arc_runtime()
            ));
            assert!(Arc::ptr_eq(
                &source.coordinates_arc_runtime(),
                &result.molecule.coordinates_arc_runtime()
            ));
            assert_eq!(
                Arc::ptr_eq(
                    &source.properties_arc_runtime(),
                    &result.molecule.properties_arc_runtime()
                ),
                existing
            );
        }
    }
}

impl MmffOptimizeMoleculeResult {
    pub fn molecule(&self) -> &Molecule {
        &self.molecule
    }
}
impl MmffOptimizeMoleculeConfsResult {
    pub fn molecule(&self) -> &Molecule {
        &self.molecule
    }
    pub fn conformer_results(&self) -> &[MmffOptimizeMoleculeConfResult] {
        &self.conformer_results
    }
}
