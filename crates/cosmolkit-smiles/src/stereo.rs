//! Notation adapter imports the one shared core ordering implementation.
pub(crate) use cosmolkit_core::parser_stereo_order::{
    atom_has_fourth_valence, bond_order_as_double, chiral_atom_needs_tag_inversion,
    count_swaps_to_interconvert, insert_implicit_nontetrahedral_neighbors, invert_tetrahedral_tag,
    nontetrahedral_chiral_permutation,
};

#[cfg(test)]
pub(crate) use cosmolkit_core::parser_stereo_order::nontetrahedral_max_neighbors;
