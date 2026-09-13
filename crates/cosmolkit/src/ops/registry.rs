//! Sole declaration of molecule operations and generated contract matrices.

use cosmolkit_macros::molecule_ops;

use super::{FeatureSpec, SupportStatus};

#[cfg(test)]
pub(crate) const COW_TEST_FEATURE: FeatureSpec = FeatureSpec {
    name: "cow-runtime-test",
    category: "internal-test",
    status: SupportStatus::Experimental,
    rdkit_parity_sensitive: false,
    docs: "Internal registered-operation coverage for block-level COW.",
};

#[cfg(feature = "hydrogens")]
pub(crate) const HYDROGENS_FEATURE: FeatureSpec = FeatureSpec {
    name: "hydrogens",
    category: "chemistry",
    status: SupportStatus::SupportedWithRdkitParity,
    rdkit_parity_sensitive: true,
    docs: "Explicit-hydrogen topology transformations.",
};

#[cfg(feature = "kekulize")]
pub(crate) const KEKULIZE_FEATURE: FeatureSpec = FeatureSpec {
    name: "kekulize",
    category: "chemistry",
    status: SupportStatus::SupportedWithRdkitParity,
    rdkit_parity_sensitive: true,
    docs: "RDKit-compatible Kekule bond assignment over molecule topology.",
};

#[cfg(feature = "aromaticity")]
pub(crate) const AROMATICITY_FEATURE: FeatureSpec = FeatureSpec {
    name: "aromaticity",
    category: "chemistry",
    status: SupportStatus::SupportedWithRdkitParity,
    rdkit_parity_sensitive: true,
    docs: "RDKit-compatible aromatic atom and bond assignment over molecule topology.",
};

#[cfg(feature = "valence")]
pub(crate) const VALENCE_FEATURE: FeatureSpec = FeatureSpec {
    name: "valence",
    category: "chemistry",
    status: SupportStatus::SupportedWithRdkitParity,
    rdkit_parity_sensitive: true,
    docs: "Explicit- and implicit-valence assignment over molecule topology.",
};

#[cfg(feature = "radicals")]
pub(crate) const RADICALS_FEATURE: FeatureSpec = FeatureSpec {
    name: "radicals",
    category: "chemistry",
    status: SupportStatus::SupportedWithRdkitParity,
    rdkit_parity_sensitive: true,
    docs: "Automatic radical-electron assignment over molecule topology.",
};

#[cfg(feature = "rings")]
pub(crate) const RINGS_FEATURE: FeatureSpec = FeatureSpec {
    name: "rings",
    category: "chemistry",
    status: SupportStatus::SupportedWithRdkitParity,
    rdkit_parity_sensitive: true,
    docs: "RDKit-compatible ring and ring-family assignment over molecule topology.",
};

#[cfg(feature = "stereo")]
pub(crate) const STEREO_FEATURE: FeatureSpec = FeatureSpec {
    name: "stereo",
    category: "chemistry",
    status: SupportStatus::SupportedWithRdkitParity,
    rdkit_parity_sensitive: true,
    docs: "RDKit-compatible 3D chiral-tag assignment and potential-stereochemistry perception.",
};

#[cfg(feature = "transforms")]
pub(crate) const TRANSFORMS_FEATURE: FeatureSpec = FeatureSpec {
    name: "transforms",
    category: "coordinates",
    status: SupportStatus::SupportedWithRdkitParity,
    rdkit_parity_sensitive: true,
    docs: "RDKit-compatible detached coordinate transforms and atom-position replacement.",
};

#[cfg(feature = "sanitize")]
pub(crate) const SANITIZE_FEATURE: FeatureSpec = FeatureSpec {
    name: "sanitize",
    category: "chemistry",
    status: SupportStatus::SupportedWithRdkitParity,
    rdkit_parity_sensitive: true,
    docs: "RDKit-compatible molecule sanitization and chemistry-problem detection.",
};

molecule_ops! {
    #[cfg(feature = "sanitize")]
    op sanitize(params: &cosmolkit_core::SanitizeParams) {
        method: sanitize_with_params,
        impl_fn: crate::ops::sanitize::sanitize_impl,
        domain: topology,
        kind: weak,
        topology_edit: local,
        access: {
            read: [],
            write: [topology, properties, derived_cache],
        },
        may_mutate: [topology, properties, derived_cache],
        auto_remap: [],
        derived_effects: {
            recompute: [],
            preserve: [coordinates],
            invalidate: [
                rings,
                ring_families,
                valence,
                aromaticity,
                stereo,
                drawing,
                fingerprint,
            ],
            operation_defined: [],
        },
        cip_state: clear,
        semantic_preconditions: [],
        requires_mapping: none,
        feature: crate::ops::runtime::registry::SANITIZE_FEATURE,
        parity: required_now,
        parity_profile: "sanitize_rdkit",
        io_roundtrip: true,
        invariant_profile: "weak_sanitize_topology_state",
        default_method: sanitize,
        default_args: [&cosmolkit_core::SanitizeParams::default()],
    }

    #[cfg(feature = "kekulize")]
    op with_kekulized_bonds(params: &cosmolkit_core::KekulizeParams) {
        method: with_kekulized_bonds_with_params,
        impl_fn: crate::ops::kekulize::kekulize_bonds_impl,
        domain: topology,
        kind: weak,
        topology_edit: local,
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
        semantic_preconditions: [],
        requires_mapping: none,
        feature: crate::ops::runtime::registry::KEKULIZE_FEATURE,
        parity: required_now,
        parity_profile: "kekulize_rdkit",
        io_roundtrip: false,
        invariant_profile: "weak_kekulize_bond_assignment",
        default_method: with_kekulized_bonds,
        default_args: [&cosmolkit_core::KekulizeParams::default()],
        inplace: true,
        inplace_method: kekulize_bonds_with_params_,
        default_inplace_method: kekulize_bonds_,
    }

    #[cfg(feature = "aromaticity")]
    op with_assigned_aromaticity(params: &cosmolkit_core::AromaticityParams) {
        method: with_assigned_aromaticity_with_params,
        impl_fn: crate::ops::aromaticity::assign_aromaticity_impl,
        domain: topology,
        kind: weak,
        topology_edit: local,
        access: {
            read: [],
            write: [topology, properties, derived_cache],
        },
        may_mutate: [topology, properties, derived_cache],
        auto_remap: [],
        derived_effects: {
            recompute: [aromaticity],
            preserve: [rings, ring_families, coordinates],
            invalidate: [valence, stereo, drawing, fingerprint],
            operation_defined: [],
        },
        cip_state: clear,
        semantic_preconditions: [],
        requires_mapping: none,
        feature: crate::ops::runtime::registry::AROMATICITY_FEATURE,
        parity: required_now,
        parity_profile: "assign_aromaticity_rdkit",
        io_roundtrip: false,
        invariant_profile: "weak_aromaticity_assignment",
        default_method: with_assigned_aromaticity,
        default_args: [&cosmolkit_core::AromaticityParams::default()],
        inplace: true,
        inplace_method: assign_aromaticity_with_params_,
        default_inplace_method: assign_aromaticity_,
    }

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
        feature: crate::ops::runtime::registry::VALENCE_FEATURE,
        parity: required_now,
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
        feature: crate::ops::runtime::registry::RADICALS_FEATURE,
        parity: required_now,
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
        feature: crate::ops::runtime::registry::RINGS_FEATURE,
        parity: required_now,
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
        feature: crate::ops::runtime::registry::RINGS_FEATURE,
        parity: required_now,
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
        semantic_preconditions: [],
        requires_mapping: none,
        feature: crate::ops::runtime::registry::STEREO_FEATURE,
        parity: required_now,
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
        result_type: crate::PotentialStereoResult,
        assemble_fn: crate::ops::potential_stereo::assemble_potential_stereo_result,
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
            preserve: [rings, ring_families, valence, aromaticity, coordinates],
            invalidate: [stereo, drawing, fingerprint],
            operation_defined: [],
        },
        cip_state: clear,
        semantic_preconditions: [],
        requires_mapping: none,
        feature: crate::ops::runtime::registry::STEREO_FEATURE,
        parity: required_now,
        parity_profile: "find_potential_stereo_rdkit",
        io_roundtrip: true,
        invariant_profile: "weak_potential_stereo_cleanup",
        default_method: potential_stereo,
        default_args: [&cosmolkit_core::PotentialStereoParams::default()],
    }

    #[cfg(feature = "hydrogens")]
    op with_hydrogens(params: &cosmolkit_core::AddHsParams) {
        method: with_hydrogens_with_params,
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
            preserve: [rings, ring_families],
            invalidate: [valence, aromaticity, stereo, drawing, fingerprint],
            operation_defined: [],
        },
        cip_state: clear,
        semantic_preconditions: [],
        requires_mapping: required,
        feature: crate::ops::runtime::registry::HYDROGENS_FEATURE,
        parity: required_now,
        parity_profile: "add_hydrogens_rdkit",
        io_roundtrip: false,
        invariant_profile: "strong_topology_with_coordinates",
        default_method: with_hydrogens,
        default_args: [&cosmolkit_core::AddHsParams::default()],
        inplace: true,
        inplace_method: add_hydrogens_with_params_,
        default_inplace_method: add_hydrogens_,
    }

    #[cfg(feature = "hydrogens")]
    op without_hydrogens(params: &cosmolkit_core::RemoveHsParams) {
        method: without_hydrogens_with_params,
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
            invalidate: [rings, ring_families, aromaticity, stereo, drawing, fingerprint],
            operation_defined: [valence],
        },
        cip_state: clear,
        semantic_preconditions: [],
        requires_mapping: required,
        feature: crate::ops::runtime::registry::HYDROGENS_FEATURE,
        parity: required_now,
        parity_profile: "remove_hydrogens_rdkit",
        io_roundtrip: false,
        invariant_profile: "strong_topology_with_coordinates",
        default_method: without_hydrogens,
        default_args: [&cosmolkit_core::RemoveHsParams::default()],
        inplace: true,
        inplace_method: remove_hydrogens_with_params_,
        default_inplace_method: remove_hydrogens_,
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
            read: [],
            write: [topology, coordinates, properties, derived_cache],
        },
        may_mutate: [topology, coordinates, properties, derived_cache],
        auto_remap: [],
        derived_effects: {
            recompute: [],
            preserve: [rings, ring_families, valence, aromaticity, fingerprint],
            invalidate: [stereo, drawing],
            operation_defined: [],
        },
        cip_state: clear,
        requires_mapping: none,
        feature: crate::ops::runtime::registry::TRANSFORMS_FEATURE,
        parity: required_now,
        parity_profile: "set_atom_position_rdkit",
        io_roundtrip: true,
        invariant_profile: "coordinate_atom_position",
        default_method: with_atom_position,
        default_args: [&cosmolkit_core::AtomPositionParams::default()],
        inplace: true,
        inplace_method: set_atom_position_with_params_,
        default_inplace_method: set_atom_position_,
    }

    #[cfg(test)]
    op cow_coordinates_for_test {
        method: cow_coordinates_for_test,
        impl_fn: crate::ops::cow_tests::cow_coordinates_for_test_impl,
        domain: coordinate,
        kind: weak,
        topology_edit: none,
        access: { read: [], write: [coordinates] },
        may_mutate: [coordinates],
        auto_remap: [],
        derived_effects: {
            recompute: [], preserve: [], invalidate: [], operation_defined: [],
        },
        cip_state: preserve,
        semantic_preconditions: [],
        requires_mapping: none,
        feature: crate::ops::runtime::registry::COW_TEST_FEATURE,
        parity: not_applicable,
        io_roundtrip: false,
        invariant_profile: "cow-coordinate-write-test",
    }

    #[cfg(test)]
    op cow_coordinates_failure_for_test {
        method: cow_coordinates_failure_for_test,
        impl_fn: crate::ops::cow_tests::cow_coordinates_failure_for_test_impl,
        domain: coordinate,
        kind: weak,
        topology_edit: none,
        access: { read: [], write: [coordinates] },
        may_mutate: [coordinates],
        auto_remap: [],
        derived_effects: {
            recompute: [], preserve: [], invalidate: [], operation_defined: [],
        },
        cip_state: preserve,
        semantic_preconditions: [],
        requires_mapping: none,
        feature: crate::ops::runtime::registry::COW_TEST_FEATURE,
        parity: not_applicable,
        io_roundtrip: false,
        invariant_profile: "cow-coordinate-failure-test",
        inplace: true,
        inplace_method: cow_coordinates_failure_for_test_,
    }
}
