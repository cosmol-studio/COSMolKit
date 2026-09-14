//! Compile-only probes for the real operation-body/runtime module boundary.
//!
//! This module is a sibling of `runtime`, exactly like production operation
//! bodies. It is enabled only by the integration compile test.

#![allow(dead_code, unexpected_cfgs)]

#[cfg(any(
    cosmolkit_runtime_privacy_case = "sanitize_allowed",
    cosmolkit_runtime_privacy_case = "sanitize_forbidden"
))]
use crate::SanitizeAccess;
#[cfg(any(
    cosmolkit_runtime_privacy_case = "allowed",
    cosmolkit_runtime_privacy_case = "forbidden"
))]
use crate::{MultiOutputOpParts, WithHydrogensAccess};
use crate::{OpParts, OperationError};

#[cfg(cosmolkit_runtime_privacy_case = "pending_allowed")]
#[cosmolkit_macros::mol_op_body(potential_stereo, context)]
fn pending_result_body() -> Result<
    crate::PotentialStereoResult<crate::PendingMolecule<super::PotentialStereoAccess>>,
    OperationError,
> {
    Ok(crate::PotentialStereoResult {
        stereo: Vec::new(),
        atom_ranks: Vec::new(),
        ring_relations: Vec::new(),
        cleaned_molecule: Some(context.pending_molecule()?),
    })
}

#[cfg(cosmolkit_runtime_privacy_case = "pending_forbidden")]
fn unregistered_pending_capability_is_hidden(parts: &mut OpParts<'_, crate::WithHydrogensAccess>) {
    let _ = parts.pending_molecule();
}

#[cfg(cosmolkit_runtime_privacy_case = "pending_forbidden")]
fn pending_runtime_lifecycle_is_hidden(mut parts: OpParts<'_, super::PotentialStereoAccess>) {
    let _ = parts.pending_molecule_runtime();
    let _ = parts.ensure_unsealed_runtime();
    let _ = parts.finish_result(crate::PotentialStereoResult::<
        crate::PendingMolecule<super::PotentialStereoAccess>,
    > {
        stereo: Vec::new(),
        atom_ranks: Vec::new(),
        ring_relations: Vec::new(),
        cleaned_molecule: None,
    });
}

#[cfg(cosmolkit_runtime_privacy_case = "pending_forbidden")]
fn pending_is_not_a_molecule(pending: crate::PendingMolecule<super::PotentialStereoAccess>) {
    let _ = &pending.identity;
    let _ = &pending.topology;
    let _ = &pending.coordinates;
    let _ = &pending.properties;
    let _ = &pending.derived_cache;
    let _ = pending.clone();
    let _ = pending.num_atoms();
    let _: crate::Molecule = pending;
}

#[cfg(cosmolkit_runtime_privacy_case = "pending_finalizer")]
fn finalizer_cannot_be_constructed() {
    let _ = crate::ResultFinalizer::<super::PotentialStereoAccess> {
        parts: None,
        operation: "forged",
    };
}

#[cfg(cosmolkit_runtime_privacy_case = "pending_wrong_marker")]
fn pending_marker_cannot_be_changed(
    pending: crate::PendingMolecule<super::PotentialStereoAccess>,
) -> crate::PendingMolecule<crate::WithHydrogensAccess> {
    pending
}

#[cfg(cosmolkit_runtime_privacy_case = "pending_reuse")]
fn pending_cannot_be_used_twice(pending: crate::PendingMolecule<super::PotentialStereoAccess>) {
    let first = Some(pending);
    let second = Some(pending);
    drop((first, second));
}

#[cfg(cosmolkit_runtime_privacy_case = "allowed")]
fn declaration_generated_capabilities_are_visible(
    parts: &mut OpParts<'_, WithHydrogensAccess>,
) -> Result<(), OperationError> {
    let topology = parts.checkout_topology()?;
    parts.install_topology(topology)?;
    let coordinates = parts.checkout_coordinates()?;
    parts.install_coordinates(coordinates)?;
    let properties = parts.checkout_properties()?;
    parts.install_properties(properties)?;
    let cache = parts.checkout_derived_cache()?;
    parts.install_derived_cache(cache)?;
    parts.record_topology_edit(crate::TopologyEditKind::Appending)?;
    parts.record_topology_mapping(cosmolkit_model::TopologyMapping::identity(0, 0))?;
    parts.apply_runtime_remap()?;
    parts.clear_cache(crate::DerivedState::RINGS)?;
    parts.apply_cip_policy()
}

#[cfg(cosmolkit_runtime_privacy_case = "forbidden")]
fn unauthorized_generated_method_is_hidden(parts: &OpParts<'_, WithHydrogensAccess>) {
    let _ = parts.properties();
}

#[cfg(cosmolkit_runtime_privacy_case = "forbidden")]
fn unrestricted_single_output_primitives_are_hidden(parts: &mut OpParts<'_, WithHydrogensAccess>) {
    let _ = parts.read_properties_runtime();
    let _ = parts.checkout_topology_runtime();
    let _ = parts.checkout_coordinates_runtime();
    let _ = parts.checkout_properties_runtime();
    let _ = parts.checkout_derived_cache_runtime();
    let _ = OpParts::<WithHydrogensAccess>::install_topology_runtime;
    let _ = OpParts::<WithHydrogensAccess>::install_coordinates_runtime;
    let _ = OpParts::<WithHydrogensAccess>::install_properties_runtime;
    let _ = OpParts::<WithHydrogensAccess>::install_derived_cache_runtime;
}

#[cfg(cosmolkit_runtime_privacy_case = "forbidden")]
fn unrestricted_multiple_output_primitives_are_hidden(
    parts: &mut MultiOutputOpParts<'_, WithHydrogensAccess>,
) {
    let _ = parts.source_topology_runtime();
    let _ = parts.source_coordinates_runtime();
    let _ = parts.source_properties_runtime();
    let _ = MultiOutputOpParts::<WithHydrogensAccess>::emit_all_runtime;
}

#[cfg(cosmolkit_runtime_privacy_case = "forbidden")]
fn runtime_fields_are_hidden(parts: &OpParts<'_, WithHydrogensAccess>) {
    let _ = &parts.spec;
    let _ = &parts.source;
    let _ = &parts.topology;
    let _ = &parts.coordinates;
    let _ = &parts.properties;
    let _ = &parts.derived_cache;
    let _ = &parts.in_place_target;
}

#[cfg(cosmolkit_runtime_privacy_case = "forbidden")]
fn wrapper_owned_lifecycle_is_hidden() {
    let _ = OpParts::<WithHydrogensAccess>::new;
    let _ = OpParts::<WithHydrogensAccess>::new_in_place;
    let _ = OpParts::<WithHydrogensAccess>::finish;
    let _ = OpParts::<WithHydrogensAccess>::abort_in_place;
    let _ = OpParts::<WithHydrogensAccess>::finish_in_place;
    let _ = MultiOutputOpParts::<WithHydrogensAccess>::new;
    let _ = MultiOutputOpParts::<WithHydrogensAccess>::finish;
}

#[cfg(cosmolkit_runtime_privacy_case = "sanitize_allowed")]
fn sanitize_declaration_generated_capabilities_are_visible(
    parts: &mut OpParts<'_, SanitizeAccess>,
) -> Result<(), OperationError> {
    let topology = parts.checkout_topology()?;
    parts.install_topology(topology)?;
    parts.record_topology_edit(crate::TopologyEditKind::Local)?;
    parts.clear_cache(
        crate::DerivedState::RINGS
            .union(crate::DerivedState::RING_FAMILIES)
            .union(crate::DerivedState::VALENCE)
            .union(crate::DerivedState::AROMATICITY)
            .union(crate::DerivedState::STEREO)
            .union(crate::DerivedState::DRAWING)
            .union(crate::DerivedState::FINGERPRINT),
    )?;
    parts.prove_preserved(
        crate::DerivedState::COORDINATES,
        crate::PreservationProof::SanitizeTopologyState,
    )?;
    parts.apply_cip_policy()
}

#[cfg(cosmolkit_runtime_privacy_case = "sanitize_forbidden")]
fn sanitize_unauthorized_generated_methods_are_hidden(parts: &mut OpParts<'_, SanitizeAccess>) {
    let _ = parts.topology();
    let _ = parts.coordinates();
    let _ = parts.properties();
    let _ = parts.derived_cache();
    let _ = parts.checkout_coordinates();
    let _ = parts.install_coordinates(cosmolkit_model::CoordinateBlock::default());
}

#[cfg(cosmolkit_runtime_privacy_case = "sanitize_forbidden")]
fn sanitize_unrestricted_runtime_and_lifecycle_are_hidden(parts: &mut OpParts<'_, SanitizeAccess>) {
    let _ = parts.read_topology_runtime();
    let _ = parts.read_properties_runtime();
    let _ = parts.checkout_topology_runtime();
    let _ = parts.install_topology_runtime(cosmolkit_model::TopologyBlock::default());
    let _ = &parts.spec;
    let _ = &parts.source;
    let _ = &parts.topology;
    let _ = &parts.coordinates;
    let _ = &parts.properties;
    let _ = &parts.derived_cache;
    let _ = &parts.in_place_target;
    let _ = OpParts::<SanitizeAccess>::new;
    let _ = OpParts::<SanitizeAccess>::new_in_place;
    let _ = OpParts::<SanitizeAccess>::finish;
    let _ = OpParts::<SanitizeAccess>::abort_in_place;
    let _ = OpParts::<SanitizeAccess>::finish_in_place;
}
