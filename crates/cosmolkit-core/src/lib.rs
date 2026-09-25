//! Detached chemistry algorithms for the crate migration.
//!
//! The crate is intentionally below the public runtime boundary. It accepts
//! owned values from `cosmolkit-model` and never accepts a live `Molecule`, an
//! operation context, or runtime cache state.

mod aromaticity;
mod atropisomer;
mod attachment_points;
mod bond_dirs;
mod cip_ranks;
mod cleanup;
mod conjugation;
mod double_stereo;
mod hcount;
mod hybridization;
mod hydrogens;
mod kekulize;
mod legacy_stereo;
mod matrices;
mod nontetrahedral_stereo;
mod paths;
mod periodic_table;
mod potential_stereo;
mod query_ops;
mod radicals;
mod rings;
mod sanitize;
mod stereo_order;
mod structure_tags;
mod transforms;
mod valence;

pub use attachment_points::{
    AttachmentExpansionError, AttachmentExpansionResult, AttachmentWarning,
    expand_attachment_points,
};

pub use atropisomer::{
    AtropisomerAssignment, AtropisomerBondUpdate, AtropisomerCarrierEnd, AtropisomerConformer,
    AtropisomerDiagnostic, AtropisomerError, AtropisomerRejectionKind, AtropisomerWedgeAssignment,
    AtropisomerWedgeUpdate, StereoGroupAssignment, atropisomer_carriers,
    cleanup_atropisomer_stereo_groups, detect_atropisomer_chirality,
    does_topology_have_atropisomers, stereo_group_atom_ids, wedge_bonds_from_atropisomers,
};

pub use aromaticity::{
    AromaticityAssignment, AromaticityError, AromaticityModel, AromaticityParams,
    assign_aromaticity, assign_aromaticity_with_query_state,
};

pub use bond_dirs::{BondDirectionStereoError, assign_chiral_types_from_bond_dirs};

pub use hydrogens::{
    AddHsParams, AddHydrogensResult, HydrogenError, HydrogenWarning, RemoveHsParams,
    RemoveHydrogensResult, add_hydrogens_impl, add_hydrogens_topology_with_query_state,
    add_hydrogens_with_params, add_hydrogens_with_query_state, remove_hydrogens_impl,
    remove_hydrogens_with_params, remove_hydrogens_with_query_state,
};

/// Strict-build-only bridge for the detached AddHs topology migration target.
///
/// The topology addition plan is internal algorithm composition state and is
/// not a public facade or binding surface.
#[cfg(feature = "op-contracts-strict")]
#[doc(hidden)]
pub mod __migration_hydrogens {
    pub use crate::attachment_points::{attachment_query_rows, expand_attachment_points};
    pub use crate::hydrogens::{
        AddHydrogensCoordinateResult, AddHydrogensTopologyResult, AddedHydrogen, AddedHydrogenKind,
        HydrogenError, PreparedHydrogenRemoval, add_hydrogen_coordinates, add_hydrogens_topology,
        place_terminal_attachment_coordinates, prepare_hydrogen_removal_stereo,
        remove_hydrogen_candidates, remove_hydrogen_candidates_with_query_state,
    };
}

pub use cip_ranks::{
    CipRankError, assign_atom_cip_ranks, assign_atom_cip_ranks_with_query_state,
    refine_atom_cip_ranks_from_invariants, refine_atom_cip_ranks_from_invariants_with_query_state,
};
pub(crate) use cleanup::{CleanupError, CleanupParams, cleanup};

/// Strict-build-only bridge for detached cleanup migration validation.
///
/// Cleanup remains a private sanitize phase and is not a public facade API.
#[cfg(feature = "op-contracts-strict")]
#[doc(hidden)]
pub mod __migration_cleanup {
    pub use crate::cleanup::{CleanupError, CleanupParams, cleanup};
}
pub(crate) use conjugation::{ConjugationError, assign_conjugation, atom_has_conjugated_bond};

/// Strict-build-only bridge for detached migration validation.
///
/// The conjugation phase remains absent from the default public surface and
/// from the `cosmolkit` facade. This exact bridge exists only so the owning
/// crate's integration target can validate the crate-private sanitize phase.
#[cfg(feature = "op-contracts-strict")]
#[doc(hidden)]
pub mod __migration_conjugation {
    pub use crate::conjugation::{ConjugationError, assign_conjugation, atom_has_conjugated_bond};
}
pub use double_stereo::{
    DoubleBondControl, DoubleBondStereoAssignment, DoubleBondStereoDescriptor,
    DoubleBondStereoError, DoubleBondStereoInfo, DoubleBondStereoSpecified,
    assign_directional_double_bond_stereo, assign_double_bond_stereo_from_directions,
    clear_bond_directions, clear_single_bond_directions, double_bond_stereo_info,
    find_double_bond_stereo_atoms, has_stereo_bond_direction, is_double_bond_stereo_candidate,
    neighboring_directed_bond, opposite_stereo_bond_direction, set_double_bond_neighbor_directions,
    should_detect_double_bond_stereo, translate_ez_to_cis_trans, with_double_bond_stereo_reference,
};

pub use hcount::total_hydrogen_count;
pub(crate) use hybridization::{HybridizationAssignment, HybridizationError, assign_hybridization};
pub use nontetrahedral_stereo::{
    non_tetrahedral_across_ligand, non_tetrahedral_ideal_angle, trigonal_bipyramidal_axial_ligand,
};

/// Strict-build-only bridge for detached hybridization migration validation.
///
/// The hybridization phase remains absent from the default public surface and
/// from the `cosmolkit` facade. This exact bridge exists only so the owning
/// crate's integration target can validate the crate-private sanitize phase.
#[cfg(feature = "op-contracts-strict")]
#[doc(hidden)]
pub mod __migration_hybridization {
    pub use crate::hybridization::{
        HybridizationAssignment, HybridizationError, assign_hybridization,
    };
}

/// Strict-build-only bridge for detached sanitize-substage migration tests.
///
/// These values remain implementation details of sanitization and are not a
/// default facade or binding surface.
#[cfg(feature = "op-contracts-strict")]
#[doc(hidden)]
pub mod __migration_sanitize {
    pub use crate::atropisomer::cleanup_invalid_atropisomers;
    pub use crate::hcount::{AdjustHsAssignment, AdjustHsError, adjust_hs};
    pub use crate::hybridization::HybridizationAssignment;
    pub use crate::sanitize::{
        SanitizeAssignment, SanitizeError, SanitizeOperations, SanitizeParams, SanitizeStage,
        sanitize_topology,
    };
    pub use crate::structure_tags::cleanup_chirality;
}

pub use kekulize::{
    CanonicalRankError, CanonicalRankParams, KekulizeAssignment, KekulizeAttempt, KekulizeError,
    KekulizeParams, kekulize, kekulize_if_possible, kekulize_if_possible_with_query_state,
    kekulize_with_query_state, rank_fragment_atoms, rank_mol_atoms_with_params,
};

pub use legacy_stereo::{
    LegacyStereoError, assign_legacy_stereochemistry, assign_legacy_stereochemistry_for_depiction,
    assign_legacy_stereochemistry_with_flags, assign_legacy_stereochemistry_with_query_state,
};

pub use structure_tags::cleanup_stereo_groups;

pub use matrices::{
    AdjacencyMatrixParams, DenseMatrix, DistanceMatrix3dParams, MatrixError,
    TopologicalDistanceMatrixParams, adjacency_matrix, distance_matrix_3d,
    topological_distance_matrix,
};

pub use paths::{
    AtomEnvironment, AtomEnvironmentParams, ConnectedComponents, DetachedPathSubgraph, GraphPath,
    PathError, PathRepresentation, PathSearchParams, SubgraphSearchParams, SubtopologyParams,
    SubtopologyResult, UniqueSubgraphParams, all_paths_in_range, all_paths_of_length,
    all_subgraphs_in_range, all_subgraphs_of_length, atom_environment, bond_ids_from_atom_path,
    connected_components, shortest_path, subtopology_from_path, unique_subgraphs_of_length,
};

pub use periodic_table::{
    PeriodicTableError, atomic_mass, element_info, isotope_abundance, isotope_mass,
    most_common_isotope, most_common_isotope_mass,
};

pub use potential_stereo::{
    PotentialStereoAssignment, PotentialStereoCenter, PotentialStereoDescriptor,
    PotentialStereoError, PotentialStereoInfo, PotentialStereoParams, PotentialStereoSpecified,
    PotentialStereoType, RingStereoRelation, potential_stereo,
};

pub use radicals::{RadicalAssignment, RadicalDiagnostic, RadicalError, assign_radicals};

pub use rings::{
    RingFindType, RingFindingError, RingInfo, RingSearchParams, fast_find_rings,
    fast_find_rings_from_parts, find_ring_families, find_ring_families_from_parts, find_sssr,
    find_sssr_from_parts, find_sssr_with_options_from_parts, is_atom_bridgehead_from_topology,
    symmetrize_sssr_with_options_from_parts, symmetrized_sssr,
};

pub use sanitize::{
    ChemistryProblem, ChemistryProblemError, ChemistryProblemReport, SanitizeAssignment,
    SanitizeError, SanitizeOperations, SanitizeParams, SanitizeStage, detect_chemistry_problems,
    sanitize_topology, sanitize_topology_with_query_state,
};
pub(crate) use sanitize::{
    PropertyCacheAssignment, PropertyCacheError, PropertyCacheParams, assign_property_cache,
};

/// Strict-build-only bridge for detached property-cache migration validation.
///
/// Property-cache calculation remains a private sanitize phase and this bridge
/// exposes no live cache installation or facade API.
#[cfg(feature = "op-contracts-strict")]
#[doc(hidden)]
pub mod __migration_property_cache {
    pub use crate::sanitize::{
        PropertyCacheAssignment, PropertyCacheError, PropertyCacheParams, assign_property_cache,
    };
}

pub use stereo_order::{
    StereoOrderError, TetrahedralLigand, TetrahedralRemap, atom_nonzero_degree,
    bond_affects_atom_chirality, count_swaps_to_interconvert, incident_tetrahedral_bond_order,
    invert_tetrahedral_tag, remap_tetrahedral_center, tetrahedral_tag_after_order_change,
};

pub use structure_tags::{
    StereoError, StructureTagAssignment, StructureTagParams, assign_chiral_tags_from_structure,
};

pub use transforms::{
    AtomPositionParams, CanonicalTransformParams, CentroidParams, PrincipalAxesAndMoments,
    PrincipalAxesKind, PrincipalAxesParams, Transform3D, TransformError, angle_degrees,
    angle_radians, bond_length, canonical_transform, canonicalize_conformer, centroid,
    dihedral_degrees, dihedral_radians, principal_axes_and_moments, transform_conformer,
    with_angle_degrees, with_angle_radians, with_atom_position, with_bond_length,
    with_dihedral_degrees, with_dihedral_radians,
};

pub use valence::{
    ValenceAssignment, ValenceError, ValenceModel, ValenceParams, ValencePhase,
    assign_explicit_valence_for_atom_from_parts,
    assign_implicit_valence_for_atom_from_parts_with_explicit_valence, assign_valence,
    assign_valence_for_topology, assign_valence_state_for_atom_from_parts,
    assign_valence_with_options_for_topology, assign_valence_with_options_from_parts,
    atom_has_valence_violation_for_topology, atom_has_valence_violation_from_parts,
    bond_type_as_double, bond_valence_contrib, calculate_explicit_valence_for_topology,
    calculate_explicit_valence_from_parts, calculate_implicit_valence_for_topology,
    calculate_implicit_valence_from_parts, can_be_hypervalent, explicit_valence_for_atom,
    get_effective_atomic_num, has_valence_violation, implicit_valence_for_atom,
    periodic_table_more_electronegative, periodic_table_outer_electrons, periodic_table_row,
    rdkit_atomic_number_from_symbol, rdkit_default_valence, rdkit_element_symbol, rdkit_rb0,
    rdkit_valence_list, required_valence_list,
};
