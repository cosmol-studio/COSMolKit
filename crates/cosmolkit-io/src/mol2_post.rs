//! MOL2-specific detached finalization; all chemistry reuses existing core owners.
use crate::{Mol2Record, MolPostError};
use cosmolkit_core::{SanitizeOperations as S, SanitizeParams};

#[derive(Debug, thiserror::Error)]
#[error("MOL2 {stage} failed: {source}")]
pub struct Mol2PostError {
    pub stage: &'static str,
    #[source]
    source: Box<dyn std::error::Error + Send + Sync>,
}
impl Mol2PostError {
    fn stage<E: std::error::Error + Send + Sync + 'static>(stage: &'static str, source: E) -> Self {
        Self {
            stage,
            source: Box::new(source),
        }
    }
}

pub fn finish_mol2_record(
    mut record: Mol2Record,
    sanitize: bool,
    remove_hydrogens: bool,
) -> Result<Mol2Record, Mol2PostError> {
    // RDKit❗❌:   // set chirality prior to sanitization since it happens from 3D and it's not
    // RDKit❗❌:   // possible anymore once the hydrogens are removed
    // RDKit❗❌:   // FIX: for now this is only for the first conformer - need to be changed once
    // RDKit❗❌:   // we use multiconformer files
    // RDKit❗❌:   MolOps::assignChiralTypesFrom3D(*res);
    // RDKit❗❌:
    // RDKit❗❌:   if (res && params.sanitize) {
    // RDKit❗❌:     MolOps::cleanUp(*res);
    // RDKit❗❌:
    // RDKit❗❌:     try {
    // RDKit❗❌:       // when we sanitize for mol2, we skip the cleanup organometallic step since it's
    // RDKit❗❌:       // not really compatible with the semantics of mol2 files
    // RDKit❗❌:       constexpr auto sanitizeFlags = MolOps::SanitizeFlags::SANITIZE_ALL ^
    // RDKit❗❌:                             MolOps::SanitizeFlags::SANITIZE_CLEANUP_ORGANOMETALLICS;
    // RDKit❗❌:       if (params.removeHs) {
    // RDKit❗❌:         // Bond stereo detection must happen before H removal, or
    // RDKit❗❌:         // else we might be removing stereogenic H atoms in double
    // RDKit❗❌:         // bonds (e.g. imines). But before we run stereo detection,
    // RDKit❗❌:         // we need to run mol cleanup so don't have trouble with
    // RDKit❗❌:         // e.g. nitro groups. Sadly, this a;; means we will find
    // RDKit❗❌:         // run both cleanup and ring finding twice (a fast find
    // RDKit❗❌:         // rings in bond stereo detection, and another in
    // RDKit❗❌:         // sanitization's SSSR symmetrization).
    // RDKit❗❌:         unsigned int failedOp = 0;
    // RDKit❗❌:         MolOps::sanitizeMol(*res, failedOp, MolOps::SanitizeFlags::SANITIZE_CLEANUP);
    // RDKit❗❌:         MolOps::detectBondStereochemistry(*res);
    // RDKit❗❌:         MolOps::RemoveHsParameters rhp;
    // RDKit❗❌:         bool sanitize = false;
    // RDKit❗❌:         MolOps::removeHs(*res, rhp, sanitize);
    // RDKit❗❌:         MolOps::sanitizeMol(*res, failedOp, sanitizeFlags);
    // RDKit❗❌:       } else {
    // RDKit❗❌:         unsigned int failedOp;
    // RDKit❗❌:         MolOps::sanitizeMol(*res, failedOp, sanitizeFlags);
    // RDKit❗❌:         MolOps::detectBondStereochemistry(*res);
    // RDKit❗❌:       }
    // RDKit❗❌:
    // RDKit❗❌:     } catch (MolSanitizeException &se) {
    // RDKit❗❌:       BOOST_LOG(rdWarningLog) << "sanitize ";
    // RDKit❗❌:       std::string molName;
    // RDKit❗❌:       res->getProp(common_properties::_Name, molName);
    // RDKit❗❌:       BOOST_LOG(rdWarningLog) << molName << ": ";
    // RDKit❗❌:       throw se;
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     res->updatePropertyCache(false);
    // RDKit❗❌:     MolOps::assignStereochemistry(*res, true, true);
    // RDKit❗❌:   }
    // Behavior: source-ordered 3D tagging, MOL2 organometallic exclusion,
    // cleanup/double-bond perception before removal, and final legacy stereo.
    // Complexity: detached owner assignments allocate validated output blocks;
    // this is more work than the source's in-place RWMol stages (❌).
    let valence = cosmolkit_core::assign_valence(
        &record.topology,
        &cosmolkit_core::ValenceParams {
            model: cosmolkit_core::ValenceModel::RdkitLike,
            strict: false,
        },
    )
    .map_err(|e| Mol2PostError::stage("Valence", e))?;
    record.topology = cosmolkit_core::assign_chiral_tags_from_structure(
        &record.topology,
        &record.coordinates,
        &valence,
        &cosmolkit_core::StructureTagParams::default(),
    )
    .map_err(|e| Mol2PostError::stage("StructureStereo", e))?
    .topology;
    if sanitize {
        // ALL's reserved bits trigger no stages; enumerate the actual named
        // source stages except CLEANUP_ORGANOMETALLICS, retaining every other
        // currently modeled stage. The core's raw-bit validator stays intact.
        let operations = S::CLEANUP
            | S::PROPERTIES
            | S::SYMM_RINGS
            | S::KEKULIZE
            | S::FIND_RADICALS
            | S::SET_AROMATICITY
            | S::SET_CONJUGATION
            | S::SET_HYBRIDIZATION
            | S::CLEANUP_CHIRALITY
            | S::ADJUST_HS
            | S::CLEANUP_ATROPISOMERS;
        record.topology = cosmolkit_core::sanitize_topology(
            &record.topology,
            &SanitizeParams {
                operations: S::CLEANUP,
            },
        )
        .map_err(|e| Mol2PostError::stage("Cleanup", e))?
        .topology;
        if remove_hydrogens {
            record.topology = cosmolkit_core::sanitize_topology(
                &record.topology,
                &SanitizeParams {
                    operations: S::CLEANUP,
                },
            )
            .map_err(|e| Mol2PostError::stage("CleanupBeforeHydrogens", e))?
            .topology;
            record.topology = crate::mol_post::detect_double_bond_stereochemistry(
                record.topology,
                &record.coordinates,
            )
            .map_err(|e: MolPostError| Mol2PostError::stage("BondStereo", e))?;
            let result = cosmolkit_core::remove_hydrogens_with_params(
                record.topology,
                record.coordinates,
                record.properties,
                &cosmolkit_core::RemoveHsParams {
                    sanitize: false,
                    ..Default::default()
                },
            )
            .map_err(|e| Mol2PostError::stage("RemoveHydrogens", e))?;
            record.topology = result.topology;
            record.coordinates = result.coordinates;
            record.properties = result.properties;
            record.topology =
                cosmolkit_core::sanitize_topology(&record.topology, &SanitizeParams { operations })
                    .map_err(|e| Mol2PostError::stage("Sanitize", e))?
                    .topology;
        } else {
            record.topology =
                cosmolkit_core::sanitize_topology(&record.topology, &SanitizeParams { operations })
                    .map_err(|e| Mol2PostError::stage("Sanitize", e))?
                    .topology;
            record.topology = crate::mol_post::detect_double_bond_stereochemistry(
                record.topology,
                &record.coordinates,
            )
            .map_err(|e| Mol2PostError::stage("BondStereo", e))?;
        }
        let valence = cosmolkit_core::assign_valence(
            &record.topology,
            &cosmolkit_core::ValenceParams {
                model: cosmolkit_core::ValenceModel::RdkitLike,
                strict: false,
            },
        )
        .map_err(|e| Mol2PostError::stage("PropertyCache", e))?;
        let rings = cosmolkit_core::symmetrized_sssr(
            &record.topology,
            &cosmolkit_core::RingSearchParams::default(),
        )
        .map_err(|e| Mol2PostError::stage("Rings", e))?;
        record.topology =
            cosmolkit_core::assign_legacy_stereochemistry(record.topology, &valence, &rings)
                .map_err(|e| Mol2PostError::stage("LegacyStereo", e))?;
        // Move the already computed final cache after all H mappings and stereo stages.
        record.post_state = crate::MolPostDerivedState {
            valence: Some(valence),
            rings: Some(rings),
        };
        record
            .properties
            .set_computed_prop("_StereochemDone", "1")
            .map_err(|e| Mol2PostError::stage("StereoMetadata", e))?;
    }
    Ok(record)
}
