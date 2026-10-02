//! SMILES-specific post-chemistry stereo dispatch over detached values.

use cosmolkit_core::{
    DoubleBondStereoError, LegacyStereoError, RingFindingError, RingSearchParams,
    ValenceAssignment, ValenceError, ValenceModel, assign_legacy_stereochemistry,
    assign_valence_state_for_atom_from_parts, assign_valence_with_options_for_topology,
    clear_single_bond_directions, set_double_bond_neighbor_directions, symmetrized_sssr,
};
use cosmolkit_model::Conformer3D;

use crate::{SmilesParseParams, SmilesRecord};

/// Structured failures from the source SMILES post-parse stereo stage.
#[derive(Debug, thiserror::Error)]
pub enum SmilesStereoError {
    #[error(transparent)]
    Properties(#[from] cosmolkit_model::MoleculePropertyError),
    #[error("prepared valence field {field} has {actual} rows; expected {expected}")]
    ValenceRows {
        field: &'static str,
        actual: usize,
        expected: usize,
    },
    #[error(transparent)]
    Directions(#[from] DoubleBondStereoError),
    #[error(transparent)]
    Assignment(#[from] LegacyStereoError),
    #[error(transparent)]
    Rings(#[from] RingFindingError),
    #[error(transparent)]
    Valence(#[from] ValenceError),
}

/// Complete stereo after the caller has performed the requested source
/// sanitize/RemoveHs stage. Never accepts or constructs a live molecule.
pub fn finalize_smiles_stereo(
    mut record: SmilesRecord,
    params: &SmilesParseParams,
    prepared_valence: &mut Option<ValenceAssignment>,
) -> Result<SmilesRecord, SmilesStereoError> {
    // RDKit SmilesParse.cpp, MolFromSmiles (2026.03.1):
    // RDKit✔️✔️:   if (res && (params.sanitize || params.removeHs)) {
    // The canonical core owners have already run the preceding removeHs or
    // sanitizeMol branch. Continue with their final topology and coordinates.
    // RDKit✔️✔️:     if (res->hasProp(SmilesParseOps::detail::_needsDetectBondStereo)) {
    // RDKit✔️✔️:       // we encountered either wiggly bond in the CXSMILES,
    // RDKit✔️✔️:       // these need to be handled the same way they were in mol files
    // RDKit✔️✔️:       if (conf || conf3d) {
    // RDKit✔️✔️:         MolOps::clearSingleBondDirFlags(*res);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       MolOps::setDoubleBondNeighborDirections(*res, conf ? conf : conf3d);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     res->clearProp(SmilesParseOps::detail::_needsDetectBondStereo);
    // RDKit✔️✔️:     // figure out stereochemistry:
    // RDKit✔️✔️:     bool cleanIt = true, force = true, flagPossible = true;
    // RDKit✔️✔️:     MolOps::assignStereochemistry(*res, cleanIt, force, flagPossible);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     //  we still need to do something about double bond stereochemistry
    // RDKit✔️✔️:     //  (was github issue 337)
    // RDKit✔️✔️:     //  now that atom stereochem has been perceived, the wedging
    // RDKit✔️✔️:     //  information is no longer needed, so we clear
    // RDKit✔️✔️:     //  single bond dir flags:
    // RDKit✔️✔️:     MolOps::clearSingleBondDirFlags(*res, true);
    // RDKit✔️✔️:   }
    // Behavior: retain the marker when BOTH flags are false. Otherwise consume
    // it only after successful direction reconstruction, then run the existing
    // fixed-profile legacy assignment (not merely Cis/Trans -> Z/E relabeling).
    // Complexity: dispatch reuses the unique core algorithms. Borrow an already
    // computed assignment; only a cold call computes all rows. Detached outputs
    // carry no runtime cache authority.
    // The optional XY lift below costs O(V); no extra topology clone is needed
    // at this boundary. Errors discard the owned record without live mutation.
    if !params.sanitize && !params.remove_hydrogens {
        record.topology = clear_single_bond_directions(record.topology, true)?;
        return Ok(record);
    }

    let rings = symmetrized_sssr(&record.topology, &RingSearchParams::default())?;
    if record.properties.prop("_needsDetectBondStereo").is_some() {
        // RDKit✔️✔️:       if (!testConf->is3D()) {
        // RDKit✔️✔️:         if (conf == nullptr) {  // only take the first 2d conf
        // RDKit✔️✔️:           conf = testConf;
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       } else {
        // RDKit✔️✔️:         if (conf3d == nullptr) {  // only take the first 3d conf
        // RDKit✔️✔️:           conf3d = testConf;
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        // RDKit's geometry kernel accepts XYZ also for a 2D conformer. Lift
        // only a borrowed 2D input, prefer it over 3D, and keep stored rows intact.
        let lifted = record.coordinates.conformers_2d.first().map(|conformer| {
            Conformer3D::new(
                conformer.id(),
                conformer
                    .coordinates()
                    .iter()
                    .map(|xy| [xy[0], xy[1], 0.0])
                    .collect(),
                false,
            )
        });
        let conformer = lifted
            .as_ref()
            .or_else(|| record.coordinates.conformers_3d.first());
        if conformer.is_some() {
            record.topology = clear_single_bond_directions(record.topology, false)?;
        }
        record.topology = set_double_bond_neighbor_directions(record.topology, &rings, conformer)?;
    }
    record.properties.clear_prop("_needsDetectBondStereo");
    // assignStereochemistry updates a missing property cache non-strictly.
    // This local assignment is not a claim that unsanitized chemistry passed
    // strict sanitization, nor does it undo CK-VALENCE-001 runtime invalidation.
    if prepared_valence.is_none() {
        *prepared_valence = Some(assign_valence_with_options_for_topology(
            &record.topology,
            ValenceModel::RdkitLike,
            false,
        )?);
    }
    let valence = prepared_valence.as_mut().expect("assigned above");
    // Validate both fields before any core stereo path indexes a prepared row.
    for (field, actual) in [
        ("explicit_valence", valence.explicit_valence.len()),
        ("implicit_hydrogens", valence.implicit_hydrogens.len()),
    ] {
        if actual != record.topology.atoms.len() {
            return Err(SmilesStereoError::ValenceRows {
                field,
                actual,
                expected: record.topology.atoms.len(),
            });
        }
    }
    // Remember only possible H-cleanup rows, not a cloned topology or an
    // additional complete valence assignment.
    let cleanup_candidates = record
        .topology
        .atoms
        .iter()
        .filter(|atom| {
            matches!(
                atom.chiral_tag(),
                cosmolkit_model::ChiralTag::TetrahedralCw
                    | cosmolkit_model::ChiralTag::TetrahedralCcw
            ) && atom.explicit_hydrogens() == 1
                && atom.formal_charge() == 0
                && !atom.is_aromatic()
        })
        .map(|atom| atom.id())
        .collect::<Vec<_>>();
    record.topology = assign_legacy_stereochemistry(record.topology, valence, &rings)?;
    // RDKit✔️✔️:       atom->setNumExplicitHs(0);
    // RDKit✔️✔️:       atom->setNoImplicit(false);
    // RDKit✔️✔️:       atom->calcExplicitValence(false);
    // RDKit✔️✔️:       atom->calcImplicitValence(false);
    // Behavior: legacy cleanup can change H state after sanitization. Refresh
    // exactly those changed rows through the existing atom-valence owner.
    // Complexity: O(V) candidate scan, O(tagged atoms) temporary IDs and O(degree)
    // per changed row; no second whole-graph assignment or topology clone.
    for id in cleanup_candidates {
        if record.topology.atoms[id.index()].explicit_hydrogens() == 0 {
            let (explicit, implicit) = assign_valence_state_for_atom_from_parts(
                &record.topology.atoms,
                &record.topology.bonds,
                &record.topology.adjacency,
                id,
                false,
            )?;
            valence.explicit_valence[id.index()] = explicit;
            valence.implicit_hydrogens[id.index()] = implicit;
        }
    }
    // RDKit Chirality.cpp::assignStereochemistry (pinned legacy profile):
    // RDKit✔️✔️:   mol.setProp(common_properties::_StereochemDone, 1, true);
    // Behavior: topology-only core dispatch cannot write molecule properties.
    // Transport the successful computed marker here, so downstream consumers
    // preserve source hasProp() semantics without repeating stereo assignment.
    // Complexity: one property insertion; no topology copy or perception pass.
    record
        .properties
        .set_computed_prop("_StereochemDone", "1")?;
    Ok(record)
}
