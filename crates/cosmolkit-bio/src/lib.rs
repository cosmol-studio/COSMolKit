//! Detached structural-biology value and algorithm boundaries.

mod hierarchy;
mod metadata;
mod protein;
mod relationships;
mod residue;
mod secondary_structure;
mod source_ids;
mod structure_metadata;

pub use hierarchy::{
    AltLocRequest, BioAltLocGroupId, BioAssembly, BioAssemblyGenerator, BioAssemblyId,
    BioAssemblyOperator, BioAssemblySpecialKind, BioAtomId, BioAtomRow, BioCalcFlag, BioChainId,
    BioChainRow, BioCoordinateBlock, BioCoordinateFormat, BioCrystalCell, BioCrystalInfo,
    BioEntityDbRef, BioEntityId, BioEntityRow, BioModelId, BioModelRow, BioNcsOperator,
    BioResidueId, BioResidueRow, BioRowSpan, BioSiftsUnpResidue, BioStructure, BioStructureError,
    BioStructureParts, BioTransform, ChainKind, EntityKind, PolymerKind, ResidueKind,
    altloc_matches, is_same_conformer,
};
pub use metadata::{
    BioBasicRefinementInfo, BioDiffractionInfo, BioExperimentInfo, BioExperimentalCrystalInfo,
    BioMetadata, BioRefinementInfo, BioRefinementRestraint, BioReflectionsInfo,
    BioSoftwareClassification, BioSoftwareItem, BioTlsGroup, BioTlsSelection,
};
pub use protein::{
    Protein, ProteinAtomRef, ProteinChainRef, ProteinProjectionError, ProteinResidueRef,
};
pub use relationships::{
    AtomAddress, BioAsu, BioCisPep, BioConnection, BioConnectionKind, BioModRes, ResidueAddress,
};
pub use residue::{
    ResidueCode, ResidueInfo, ResidueInfoKind, ResidueSequenceError,
    UNKNOWN_TABULATED_RESIDUE_INDEX, expand_one_letter, expand_one_letter_sequence,
    find_residue_info, find_residue_info_index, residue_code, residue_info, residue_info_checked,
};
pub use secondary_structure::{BioHelix, BioHelixClass, BioSheet, BioStrand};
pub use source_ids::{
    AltLocLabel, AtomName, AtomSourceIds, ChainSourceIds, EntitySourceIds, PdbAtomSerial,
    PdbChainId, PdbSeqId, ResidueName, ResidueSourceIds,
};
pub use structure_metadata::BioStructureSourceState;

use cosmolkit_model::TopologyBlock;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BioError {
    Unsupported,
}

pub fn select_residues(topology: &TopologyBlock, query: &str) -> Result<TopologyBlock, BioError> {
    let _ = (topology, query);
    Err(BioError::Unsupported)
}
