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
    let atom_count = topology.atoms.len();
    let conformer = cosmolkit_depict::compute_2d_coordinates(topology, params)
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
        coordinates.conformers_2d.clear();
        0
    } else {
        next_2d_conformer_id(&coordinates).map_err(OperationError::Coordinate2D)?
    };
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
