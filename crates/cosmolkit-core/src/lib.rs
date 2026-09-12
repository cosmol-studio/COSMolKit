//! Detached chemistry algorithms for the crate migration.
//!
//! The crate is intentionally below the public runtime boundary. It accepts
//! owned values from `cosmolkit-model` and never accepts a live `Molecule`, an
//! operation context, or runtime cache state.

mod atropisomer;
mod cip_ranks;
mod double_stereo;
mod hcount;
mod hydrogens;
mod matrices;
mod paths;
mod periodic_table;
mod potential_stereo;
mod radicals;
mod rings;
mod sanitize;
mod stereo_order;
mod structure_tags;
mod transforms;
mod valence;

pub use atropisomer::{
    AtropisomerAssignment, AtropisomerBondUpdate, AtropisomerConformer, AtropisomerDiagnostic,
    AtropisomerError, AtropisomerRejectionKind, AtropisomerWedgeAssignment, AtropisomerWedgeUpdate,
    StereoGroupAssignment, cleanup_atropisomer_stereo_groups, detect_atropisomer_chirality,
    does_topology_have_atropisomers, stereo_group_atom_ids, wedge_bonds_from_atropisomers,
};

pub use hydrogens::{
    AddHsParams, CoreOperationError, DetachedBlocks, RemoveHsParams, add_hydrogens_impl,
    add_hydrogens_with_params, remove_hydrogens_impl, remove_hydrogens_with_params,
};

pub use cip_ranks::{CipRankError, assign_atom_cip_ranks, refine_atom_cip_ranks_from_invariants};
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
    PropertyCacheAssignment, SanitizeError, assign_property_cache_for_topology,
    assign_valence_properties_for_topology,
};

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
