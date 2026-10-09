//! Detached source-backed chemistry algorithms.
//!
//! The crate is intentionally below the public runtime boundary. It accepts
//! owned values from `cosmolkit-model` and never accepts a live `Molecule`, an
//! operation context, or runtime cache state.

mod alignment;
mod aromaticity;
mod atropisomer;
mod attachment_points;
mod bond_dirs;
mod cip_ranks;
mod cleanup;
mod conjugation;
mod coordinate_input;
mod double_stereo;
mod fragments;
mod hcount;
mod hybridization;
mod hydrogens;
mod kekulize;
mod legacy_stereo;
mod matrices;
mod nontetrahedral_stereo;
#[doc(hidden)]
pub mod parser_helpers;
#[doc(hidden)]
pub mod parser_stereo_order;
mod paths;
mod periodic_table;
mod platform_threads;
mod polymer_sgroup;
mod potential_stereo;
mod property_string;
mod query_ops;
mod radicals;
mod random;
mod rings;
mod sanitize;
mod scaffolds;
#[doc(hidden)]
pub mod source_control_c;
mod source_sort;
#[doc(hidden)]
pub mod stereo_graph;
mod stereo_order;
mod structure_tags;
mod transforms;
mod valence;
mod wedge;

#[cfg(test)]
mod tests;

pub use platform_threads::{
    ThreadCountError, cpu_online_count, observe_hardware_threads, rdkit_thread_count,
    rdkit_threads_with_observed_hardware,
};

pub use random::{RdkitRandomEngine, RdkitRandomGenerator, with_rdkit_random_generator};
pub use stereo_order::{invert_atom_chirality, invert_bond_chirality};

/// Detached, source-ordered fragment extraction with explicit topology mappings.
pub use fragments::{
    MoleculeFragment, MoleculeFragmentsError, get_largest_molecule_fragment, get_molecule_fragments,
};
pub use scaffolds::{
    ScaffoldError, ScaffoldResult, murcko_decompose, murcko_scaffold, net_scaffold,
};

/// Narrow borrowed fragment bridge used by sibling algorithm crates.
#[doc(hidden)]
pub use fragments::{
    FragmentCoordinateView, FragmentCoordinateViewError,
    get_molecule_fragments_with_coordinate_view,
};

pub use attachment_points::{
    AttachmentExpansionError, AttachmentExpansionResult, AttachmentWarning,
    expand_attachment_points,
};

pub use atropisomer::{
    AtropisomerAssignment, AtropisomerBondUpdate, AtropisomerCarrierEnd, AtropisomerConformer,
    AtropisomerDiagnostic, AtropisomerError, AtropisomerRejectionKind, AtropisomerWedgeAssignment,
    AtropisomerWedgeUpdate, StereoGroupAssignment, atropisomer_carriers,
    atropisomer_carriers_for_bonds, cleanup_atropisomer_stereo_groups,
    collect_query_stereo_group_atom_ids_source, collect_stereo_group_atom_ids_source,
    detect_atropisomer_chirality, does_topology_have_atropisomers,
    get_all_atom_ids_for_stereo_group, get_all_atom_ids_for_stereo_groups, stereo_group_atom_ids,
    wedge_atropisomer_no_conformer_source, wedge_atropisomer_three_d_source,
    wedge_atropisomer_two_d_source, wedge_bonds_from_atropisomers,
    wedge_bonds_from_atropisomers_projected_source, wedge_bonds_from_atropisomers_source,
};

pub use aromaticity::{
    AromaticityAssignment, AromaticityError, AromaticityModel, AromaticityParams,
    assign_aromaticity, assign_aromaticity_with_query_state, assign_mmff_aromaticity_prepared,
};

pub use bond_dirs::{BondDirectionStereoError, assign_chiral_types_from_bond_dirs};

pub use hydrogens::{
    AddHsParams, AddHydrogensResult, HydrogenError, HydrogenWarning, RemoveHsParams,
    RemoveHydrogensResult, add_hydrogens_impl, add_hydrogens_topology_with_query_state,
    add_hydrogens_with_params, add_hydrogens_with_query_state, add_hydrogens_with_source_valence,
    remove_hydrogens_impl, remove_hydrogens_with_params, remove_hydrogens_with_query_state,
};

pub use cip_ranks::{
    CipRankError, assign_atom_cip_ranks, assign_atom_cip_ranks_with_query_state,
    refine_atom_cip_ranks_from_invariants, refine_atom_cip_ranks_from_invariants_with_query_state,
};
pub(crate) use cleanup::{CleanupError, CleanupParams, cleanup};

pub use conjugation::{ConjugationError, assign_conjugation_flags};
pub(crate) use conjugation::{assign_conjugation, atom_has_conjugated_bond};

pub use double_stereo::{
    DoubleBondControl, DoubleBondStereoAssignment, DoubleBondStereoDescriptor,
    DoubleBondStereoError, DoubleBondStereoInfo, DoubleBondStereoSpecified, DoubleBondStereoUpdate,
    StereoAtomSearch, StereoAtomSearchWarning, assign_directional_double_bond_stereo,
    assign_double_bond_stereo_from_directions, clear_bond_directions, clear_single_bond_directions,
    double_bond_stereo_info, double_bond_stereo_reference_atoms, find_double_bond_stereo_atoms,
    find_double_bond_stereo_atoms_with_rank_reader, has_stereo_bond_direction,
    is_double_bond_stereo_candidate, neighboring_directed_bond,
    neighboring_directed_bond_from_incident, opposite_stereo_bond_direction,
    set_double_bond_neighbor_directions, should_detect_double_bond_stereo,
    translate_ez_to_cis_trans, with_double_bond_stereo_reference,
};

pub use hcount::{total_hydrogen_count, total_hydrogen_count_from_validated};
pub(crate) use hybridization::assign_hybridization;
pub use hybridization::{
    HybridizationAssignment, HybridizationError, assign_hybridization_with_conjugation,
};
pub use nontetrahedral_stereo::{
    non_tetrahedral_across_ligand, non_tetrahedral_ideal_angle, trigonal_bipyramidal_axial_ligand,
};

pub use kekulize::{
    CanonicalRankError, CanonicalRankParams, KekulizeAssignment, KekulizeAttempt, KekulizeError,
    KekulizeParams, kekulize, kekulize_if_possible, kekulize_if_possible_with_query_state,
    kekulize_if_possible_with_query_state_and_ring_info, kekulize_selected_fragment,
    kekulize_with_query_state, kekulize_with_query_state_and_ring_info, rank_fragment_atoms,
    rank_fragment_atoms_with_params, rank_fragment_atoms_with_prepared_state,
    rank_mol_atoms_with_params, source_kekulize_attempt,
};

pub use legacy_stereo::{
    LegacyStereoAssignment, LegacyStereoError, assign_legacy_stereochemistry,
    assign_legacy_stereochemistry_for_depiction, assign_legacy_stereochemistry_source,
    assign_legacy_stereochemistry_with_assignments, assign_legacy_stereochemistry_with_flags,
    assign_legacy_stereochemistry_with_query_state,
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
    connected_components, query_atom_paths_in_range, query_bond_paths_in_range,
    query_connected_components, query_subgraphs_in_range, shortest_path, subtopology_from_path,
    unique_subgraphs_of_length,
};

pub use periodic_table::{
    PeriodicTableError, atomic_mass, element_info, is_rdkit_organic_subset, isotope_abundance,
    isotope_mass, most_common_isotope, most_common_isotope_mass,
};

pub use potential_stereo::{
    PotentialStereoAssignment, PotentialStereoCenter, PotentialStereoDescriptor,
    PotentialStereoError, PotentialStereoInfo, PotentialStereoParams, PotentialStereoSpecified,
    PotentialStereoType, RingStereoRelation, potential_stereo,
    potential_tetrahedral_centers_for_atoms,
};

pub use polymer_sgroup::{
    PolymerSGroupError, finalize_polymer_sgroup, setup_unmarked_polymer_sgroup,
};
pub use property_string::{
    PropertyStringError, RequiredPropertyStringError, int_vector_to_string,
    property_value_to_string, required_property_value_to_string,
};

pub use radicals::{RadicalAssignment, RadicalDiagnostic, RadicalError, assign_radicals};

pub use rings::{
    RingFindType, RingFindingError, RingInfo, RingSearchParams, fast_find_rings,
    fast_find_rings_from_parts, find_ring_families, find_ring_families_from_parts, find_sssr,
    find_sssr_from_parts, find_sssr_with_options_from_parts,
    find_sssr_with_source_outputs_from_parts, is_atom_bridgehead_from_topology,
    ring_info_from_selected_rows, symmetrize_sssr_with_options_from_parts, symmetrized_sssr,
    symmetrized_sssr_with_properties,
};

pub use sanitize::{
    ChemistryProblem, ChemistryProblemError, ChemistryProblemReport, SanitizeAssignment,
    SanitizeError, SanitizeOperations, SanitizeParams, SanitizeStage, detect_chemistry_problems,
    sanitize_topology, sanitize_topology_with_query_state, source_sanitize_tautomer_product,
};
pub(crate) use sanitize::{
    PropertyCacheAssignment, PropertyCacheError, PropertyCacheParams, assign_property_cache,
};

pub use stereo_order::{
    StereoOrderError, TetrahedralLigand, TetrahedralRemap, atom_nonzero_degree,
    atom_perturbation_order, bond_affects_atom_chirality, count_swaps_to_interconvert,
    incident_tetrahedral_bond_order, invert_tetrahedral_tag, remap_tetrahedral_center,
    tetrahedral_tag_after_order_change,
};

pub use structure_tags::{
    StereoError, StructureTagAssignment, StructureTagParams, assign_chiral_tags_from_structure,
    nontetrahedral_enabled, unsigned_dihedral_radians,
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
    AtomMetadata, AtomicSymbolLookupError, SymbolValenceLookupError, ValenceAssignment,
    ValenceError, ValenceModel, ValenceParams, ValencePhase,
    assign_explicit_valence_for_atom_from_parts,
    assign_implicit_valence_for_atom_from_parts_with_explicit_valence, assign_valence,
    assign_valence_for_topology, assign_valence_state_for_atom_from_parts,
    assign_valence_with_options_for_topology, assign_valence_with_options_from_parts,
    atom_has_valence_violation_for_topology, atom_has_valence_violation_from_parts, atom_metadata,
    atom_metadata_from_assignment, atom_valence_cache_needs_update, bond_type_as_double,
    bond_valence_contrib, calculate_explicit_valence_for_topology,
    calculate_explicit_valence_from_parts, calculate_implicit_valence_for_topology,
    calculate_implicit_valence_from_parts, can_be_hypervalent, explicit_valence_for_atom,
    get_effective_atomic_num, has_valence_violation, implicit_valence_for_atom,
    num_pi_electrons_for_topology, periodic_table_more_electronegative,
    periodic_table_outer_electrons, periodic_table_row, rdkit_atomic_number_from_c_symbol,
    rdkit_atomic_number_from_symbol, rdkit_default_valence, rdkit_default_valence_from_c_symbol,
    rdkit_default_valence_from_symbol, rdkit_element_symbol, rdkit_rb0, rdkit_valence_list,
    rdkit_valence_list_from_c_symbol, rdkit_valence_list_from_symbol, required_valence_list,
};

pub use wedge::{
    CrossedBondContext, MolFileBondStereoInfo, WedgeAssignments, WedgeError, WedgeInfo,
    determine_bond_wedge_state, get_molfile_bond_stereo_info, pick_bonds_to_wedge,
    pick_bonds_to_wedge_default_source, pick_bonds_to_wedge_source,
    pick_bonds_to_wedge_with_existing_ring_info, pick_bonds_to_wedge_with_ring_info,
};

pub use coordinate_input::{
    Coordinate2DInputParams, Coordinate3DInputParams, Coordinate3DReadError, CoordinateInputError,
    CoordinateZPolicy, Replace3DCoordinatesParams, append_3d_conformer, clear_3d_conformers,
    coordinates_2d_from_input, coordinates_3d_for_id, coordinates_3d_from_input,
    install_2d_coordinates, install_only_3d_conformer, replace_3d_coordinates,
};

/// Foundational quaternion alignment over detached point rows.
pub use alignment::{align_points, alignment_sum_squared_residual, alignment_transform_point};

pub use periodic_table::{covalent_radius, van_der_waals_radius};

mod property_numeric;
#[doc(hidden)]
pub use property_numeric::{PropertyIntReadError, property_value_to_int};
#[doc(hidden)]
pub use property_numeric::{PropertyUIntReadError, UIntLexicalReadError, property_value_to_uint};

#[doc(hidden)]
pub use property_numeric::{
    DoubleLexicalReadError, DoubleLexicalReadErrorKind, PropertyDoubleReadError,
    SourceDoubleExtraction, SourceNumericRoundingMode, SourceNumericStreamState,
    UnsignedStreamArrayError, property_value_to_double, source_extract_double, source_field_double,
    source_lexical_double, source_lexical_double_with_rounding, source_unsigned_stream_array,
    source_unsigned_stream_read,
};

pub use fragments::{
    FragmentSourceMetadata, FragmentSourceMetadataView,
    assign_molecule_fragments_with_source_outputs, get_molecule_fragments_with_source_outputs,
    get_shared_molecule_fragments_with_source_outputs,
};

#[doc(hidden)]
pub use potential_stereo::potential_tetrahedral_center_from_source;

#[doc(hidden)]
pub use wedge::pick_bond_to_wedge_with_source_properties;

#[doc(hidden)]
pub use valence::update_query_atom_property_cache_source;

#[doc(hidden)]
pub use query_ops::query_bond_has_complex_type_query;

#[doc(hidden)]
pub use atropisomer::query_atropisomer_carriers_source;
#[doc(hidden)]
pub use wedge::get_query_directional_bond_stereo_info_source;

#[doc(hidden)]
pub use atropisomer::wedge_query_bonds_from_atropisomers_source;
#[doc(hidden)]
pub use wedge::pick_query_bonds_to_wedge_source;

#[doc(hidden)]
pub use property_numeric::{PropertyULongReadError, property_value_to_ulong};
