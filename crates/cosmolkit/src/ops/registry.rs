//! Sole declaration of molecule operations and generated contract matrices.

use cosmolkit_macros::molecule_ops;

use super::{FeatureSpec, SupportStatus};

#[cfg(feature = "hydrogens")]
pub(crate) const HYDROGENS_FEATURE: FeatureSpec = FeatureSpec {
    name: "hydrogens",
    category: "chemistry",
    status: SupportStatus::Unsupported {
        reason: "hydrogen transforms have not yet been ported to their final algorithm owner",
    },
    rdkit_parity_sensitive: true,
    docs: "Explicit-hydrogen topology transformations.",
};

#[cfg(feature = "valence")]
pub(crate) const VALENCE_FEATURE: FeatureSpec = FeatureSpec {
    name: "valence",
    category: "chemistry",
    status: SupportStatus::Unsupported {
        reason: "valence assignment is registered but its live cache integration is not complete",
    },
    rdkit_parity_sensitive: true,
    docs: "Explicit- and implicit-valence assignment over molecule topology.",
};

#[cfg(feature = "radicals")]
pub(crate) const RADICALS_FEATURE: FeatureSpec = FeatureSpec {
    name: "radicals",
    category: "chemistry",
    status: SupportStatus::Unsupported {
        reason: "radical assignment is registered but its live topology integration is not complete",
    },
    rdkit_parity_sensitive: true,
    docs: "Automatic radical-electron assignment over molecule topology.",
};

#[cfg(feature = "rings")]
pub(crate) const RINGS_FEATURE: FeatureSpec = FeatureSpec {
    name: "rings",
    category: "chemistry",
    status: SupportStatus::Unsupported {
        reason: "ring assignments are registered but their live cache integration is not complete",
    },
    rdkit_parity_sensitive: true,
    docs: "RDKit-compatible ring and ring-family assignment over molecule topology.",
};

#[cfg(feature = "stereo")]
pub(crate) const STEREO_FEATURE: FeatureSpec = FeatureSpec {
    name: "stereo",
    category: "chemistry",
    status: SupportStatus::Unsupported {
        reason: "stereochemistry operations are registered but their detached implementations and live operation integrations are not complete",
    },
    rdkit_parity_sensitive: true,
    docs: "RDKit-compatible atom structure-tag assignment, potential-stereochemistry perception, and weak topology-state cleanup.",
};

#[cfg(feature = "transforms")]
pub(crate) const TRANSFORMS_FEATURE: FeatureSpec = FeatureSpec {
    name: "transforms",
    category: "coordinates",
    status: SupportStatus::Unsupported {
        reason: "coordinate transforms are registered but their live operation integration is not complete",
    },
    rdkit_parity_sensitive: true,
    docs: "RDKit-compatible detached coordinate transforms and atom-position replacement.",
};

molecule_ops! {
    #[cfg(feature = "valence")]
    op with_assigned_valence(params: &cosmolkit_core::ValenceParams) {
        method: with_assigned_valence_with_params,
        impl_fn: crate::ops::valence::assign_valence_impl,
        domain: topology,
        kind: weak,
        topology_edit: none,
        access: {
            read: [topology],
            write: [derived_cache],
        },
        may_mutate: [derived_cache],
        auto_remap: [],
        derived_effects: {
            recompute: [valence],
            preserve: [rings, stereo],
            invalidate: [],
            operation_defined: [],
        },
        cip_state: preserve,
        requires_mapping: none,
        feature: crate::ops::registry::VALENCE_FEATURE,
        parity: required_when_supported,
        parity_profile: "assign_valence_rdkit",
        io_roundtrip: false,
        invariant_profile: "weak_valence_cache_assignment",
        default_method: with_assigned_valence,
        default_args: [&cosmolkit_core::ValenceParams::default()],
        inplace: true,
        inplace_method: assign_valence_with_params_,
        default_inplace_method: assign_valence_,
    }

    #[cfg(feature = "radicals")]
    op with_assigned_radicals {
        method: with_assigned_radicals,
        impl_fn: crate::ops::radicals::assign_radicals_impl,
        domain: topology,
        kind: weak,
        topology_edit: none,
        access: {
            read: [],
            write: [topology, properties, derived_cache],
        },
        may_mutate: [topology, properties, derived_cache],
        auto_remap: [],
        derived_effects: {
            recompute: [],
            preserve: [rings, ring_families, coordinates],
            invalidate: [valence, aromaticity, stereo, drawing, fingerprint],
            operation_defined: [],
        },
        cip_state: clear,
        requires_mapping: none,
        feature: crate::ops::registry::RADICALS_FEATURE,
        parity: required_when_supported,
        parity_profile: "assign_radicals_rdkit",
        io_roundtrip: false,
        invariant_profile: "weak_radical_assignment",
        inplace: true,
        inplace_method: assign_radicals_,
    }

    #[cfg(feature = "rings")]
    op with_assigned_rings {
        method: with_assigned_rings,
        impl_fn: crate::ops::rings::assign_rings_impl,
        domain: topology,
        kind: weak,
        topology_edit: none,
        access: {
            read: [topology],
            write: [derived_cache],
        },
        may_mutate: [derived_cache],
        auto_remap: [],
        derived_effects: {
            recompute: [rings],
            preserve: [valence, aromaticity, stereo, coordinates, drawing, fingerprint],
            invalidate: [ring_families],
            operation_defined: [],
        },
        cip_state: preserve,
        requires_mapping: none,
        feature: crate::ops::registry::RINGS_FEATURE,
        parity: required_when_supported,
        parity_profile: "fast_find_rings_rdkit",
        io_roundtrip: false,
        invariant_profile: "weak_ring_cache_assignment",
        inplace: true,
        inplace_method: assign_rings_,
    }

    #[cfg(feature = "rings")]
    op with_assigned_ring_families(params: &cosmolkit_core::RingSearchParams) {
        method: with_assigned_ring_families_with_params,
        impl_fn: crate::ops::rings::assign_ring_families_impl,
        domain: topology,
        kind: weak,
        topology_edit: none,
        access: {
            read: [topology],
            write: [derived_cache],
        },
        may_mutate: [derived_cache],
        auto_remap: [],
        derived_effects: {
            recompute: [ring_families],
            preserve: [rings, valence, aromaticity, stereo, coordinates, drawing, fingerprint],
            invalidate: [],
            operation_defined: [],
        },
        cip_state: preserve,
        requires_mapping: none,
        feature: crate::ops::registry::RINGS_FEATURE,
        parity: required_when_supported,
        parity_profile: "find_ring_families_rdkit",
        io_roundtrip: false,
        invariant_profile: "weak_ring_family_cache_assignment",
        default_method: with_assigned_ring_families,
        default_args: [&cosmolkit_core::RingSearchParams::default()],
        inplace: true,
        inplace_method: assign_ring_families_with_params_,
        default_inplace_method: assign_ring_families_,
    }

    #[cfg(feature = "stereo")]
    op with_chiral_tags_from_structure(params: &cosmolkit_core::StructureTagParams) {
        method: with_chiral_tags_from_structure_with_params,
        impl_fn: crate::ops::structure_tags::assign_chiral_tags_from_structure_impl,
        domain: topology,
        kind: weak,
        topology_edit: none,
        access: {
            read: [coordinates],
            write: [topology, properties, derived_cache],
        },
        may_mutate: [topology, properties, derived_cache],
        auto_remap: [],
        derived_effects: {
            recompute: [],
            preserve: [rings, ring_families, valence, aromaticity, coordinates],
            invalidate: [stereo, drawing, fingerprint],
            operation_defined: [],
        },
        cip_state: clear,
        semantic_preconditions: [trusted_bond_topology, hydrogen_ownership_represented],
        requires_mapping: none,
        feature: crate::ops::registry::STEREO_FEATURE,
        parity: required_when_supported,
        parity_profile: "assign_chiral_tags_from_structure_rdkit",
        io_roundtrip: true,
        invariant_profile: "weak_structure_tag_assignment",
        default_method: with_chiral_tags_from_structure,
        default_args: [&cosmolkit_core::StructureTagParams::default()],
        inplace: true,
        inplace_method: assign_chiral_tags_from_structure_with_params_,
        default_inplace_method: assign_chiral_tags_from_structure_,
    }

    #[cfg(feature = "stereo")]
    op potential_stereo(params: &cosmolkit_core::PotentialStereoParams) {
        method: potential_stereo_with_params,
        impl_fn: crate::ops::potential_stereo::potential_stereo_impl,
        result_type: cosmolkit_core::PotentialStereoAssignment,
        domain: topology,
        kind: weak,
        topology_edit: none,
        access: {
            read: [],
            write: [topology, derived_cache],
        },
        may_mutate: [topology, derived_cache],
        auto_remap: [],
        derived_effects: {
            recompute: [],
            preserve: [rings, ring_families, valence, aromaticity, coordinates],
            invalidate: [stereo, drawing, fingerprint],
            operation_defined: [],
        },
        cip_state: clear,
        semantic_preconditions: [trusted_bond_topology, hydrogen_ownership_represented],
        requires_mapping: none,
        feature: crate::ops::registry::STEREO_FEATURE,
        parity: required_when_supported,
        parity_profile: "find_potential_stereo_rdkit",
        io_roundtrip: true,
        invariant_profile: "weak_potential_stereo_cleanup",
        default_method: potential_stereo,
        default_args: [&cosmolkit_core::PotentialStereoParams::default()],
    }

    #[cfg(feature = "hydrogens")]
    op with_hydrogens {
        method: with_hydrogens,
        impl_fn: crate::ops::hydrogens::add_hydrogens_impl,
        domain: topology,
        kind: strong,
        topology_edit: expanding,
        access: {
            read: [],
            write: [topology, coordinates, properties, derived_cache],
        },
        may_mutate: [topology, coordinates, properties, derived_cache],
        auto_remap: [coordinates, properties],
        derived_effects: {
            recompute: [],
            preserve: [],
            invalidate: [rings, valence, stereo],
            operation_defined: [],
        },
        cip_state: clear,
        requires_mapping: required,
        feature: crate::ops::registry::HYDROGENS_FEATURE,
        parity: required_when_supported,
        parity_profile: "add_hydrogens_rdkit",
        io_roundtrip: false,
        invariant_profile: "strong_topology_with_coordinates",
        inplace: true,
        inplace_method: add_hydrogens_,
    }

    #[cfg(feature = "hydrogens")]
    op without_hydrogens {
        method: without_hydrogens,
        impl_fn: crate::ops::hydrogens::remove_hydrogens_impl,
        domain: topology,
        kind: strong,
        topology_edit: compacting,
        access: {
            read: [],
            write: [topology, coordinates, properties, derived_cache],
        },
        may_mutate: [topology, coordinates, properties, derived_cache],
        auto_remap: [coordinates, properties],
        derived_effects: {
            recompute: [],
            preserve: [],
            invalidate: [rings, stereo],
            operation_defined: [valence],
        },
        cip_state: clear,
        requires_mapping: required,
        feature: crate::ops::registry::HYDROGENS_FEATURE,
        parity: required_when_supported,
        parity_profile: "remove_hydrogens_rdkit",
        io_roundtrip: false,
        invariant_profile: "strong_topology_with_coordinates",
        inplace: true,
        inplace_method: remove_hydrogens_,
    }

    #[cfg(feature = "transforms")]
    op with_atom_position(
        atom: cosmolkit_model::AtomId,
        position: [f64; 3],
        params: &cosmolkit_core::AtomPositionParams,
    ) {
        method: with_atom_position_with_params,
        impl_fn: crate::ops::transforms::with_atom_position_impl,
        domain: coordinate,
        kind: weak,
        topology_edit: none,
        access: {
            read: [topology],
            write: [coordinates, derived_cache],
        },
        may_mutate: [coordinates, derived_cache],
        auto_remap: [],
        derived_effects: {
            recompute: [],
            preserve: [rings, ring_families, valence, aromaticity, fingerprint],
            invalidate: [stereo, drawing],
            operation_defined: [],
        },
        cip_state: clear,
        requires_mapping: none,
        feature: crate::ops::registry::TRANSFORMS_FEATURE,
        parity: required_when_supported,
        parity_profile: "set_atom_position_rdkit",
        io_roundtrip: true,
        invariant_profile: "coordinate_atom_position",
        default_method: with_atom_position,
        default_args: [&cosmolkit_core::AtomPositionParams::default()],
        inplace: true,
        inplace_method: set_atom_position_with_params_,
        default_inplace_method: set_atom_position_,
    }
}
