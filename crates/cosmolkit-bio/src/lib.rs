//! Detached structural-biology value and algorithm boundaries.

#![cfg_attr(
    not(feature = "test-probes"),
    doc = "Normal builds keep selection test probes private to explicit test configurations.\n\n```compile_fail\nuse cosmolkit_bio::__bio_rows_probe;\n```"
)]

mod hierarchy;
mod metadata;
mod protein;
mod relationships;
mod residue;
mod secondary_structure;
mod selection;
mod source_ids;
mod structure_metadata;

/// Test-only bridge for the owning IO integration target (BIO-ROWS R12,
/// `migration_io_bio_rows.rs`). Doc-hidden, registered in no contract,
/// and NOT public API: it exposes the private selection module's
/// row-integration entry points behind plain-data types because the
/// packet's prescribed IO gate command enables no probe feature. Removal
/// is authorized with the owning tests.
/// Test infrastructure only (BIO-ROWS-PROBE): this module exists solely
/// for the owning IO integration target and is compiled only when the
/// default-disabled `test-probes` feature is explicitly enabled.
/// Doc-hidden is NOT sufficient on its own — without the feature the
/// module does not exist at all in normal builds.
#[cfg(feature = "test-probes")]
#[doc(hidden)]
pub mod __bio_rows_probe {
    use crate::hierarchy::{BioCoordinateFormat, BioStructureData};
    use crate::selection::{
        BioRowTraverseError, BioSelectedAtomIds, Selection, SelectionList, SelectionSequenceId,
        bio_first_rows,
    };

    /// Plain-data selector description mirroring the private `Selection`
    /// field layout (subset used by the integration tests).
    #[derive(Debug, Clone, Default)]
    pub struct ProbeSelection {
        pub mdl: i32,
        pub chain_list: Option<String>,
        pub residue_names: Option<String>,
        pub atom_names: Option<String>,
        pub from_seqid: Option<(i32, u8)>,
        pub to_seqid: Option<(i32, u8)>,
        /// Explicit Gemmi El ordinals (X=0, H=1, ..., Og=118, D=119)
        /// selecting the element mask. BIO-CID C02: the probe no longer
        /// parses CID element lists — the lexical parser moved to the IO
        /// owner and BIO keeps no parser — so callers state ordinals
        /// directly.
        pub element_ordinals: Option<Vec<u8>>,
        pub altloc_list: Option<String>,
        /// Strict occupancy lower bound (property 'q', relation '>').
        pub occ_gt: Option<f64>,
    }

    impl ProbeSelection {
        fn build(&self) -> Selection {
            let mut sel = Selection::default();
            sel.mdl = self.mdl;
            if let Some(list) = &self.chain_list {
                sel.chain_ids = SelectionList {
                    all: false,
                    inverted: false,
                    list: list.clone(),
                };
            }
            if let Some(list) = &self.residue_names {
                sel.residue_names = SelectionList {
                    all: false,
                    inverted: false,
                    list: list.clone(),
                };
            }
            if let Some(list) = &self.atom_names {
                sel.atom_names = SelectionList {
                    all: false,
                    inverted: false,
                    list: list.clone(),
                };
            }
            if let Some((num, icode)) = self.from_seqid {
                sel.from_seqid = SelectionSequenceId { seqnum: num, icode };
            }
            if let Some((num, icode)) = self.to_seqid {
                sel.to_seqid = SelectionSequenceId { seqnum: num, icode };
            }
            if let Some(ordinals) = &self.element_ordinals {
                // Explicit mask construction (BIO-CID C02): no CID parsing
                // in BIO; ordinals address the canonical 120-slot mask.
                let mut mask = [false; crate::selection::GEMMI_EL_END];
                for &ordinal in ordinals {
                    mask[usize::from(ordinal)] = true;
                }
                sel.elements = Some(mask);
            }
            if let Some(list) = &self.altloc_list {
                sel.altlocs = SelectionList {
                    all: false,
                    inverted: false,
                    list: list.clone(),
                };
            }
            if let Some(value) = self.occ_gt {
                sel.atom_inequalities
                    .push(crate::selection::SelectionAtomInequality {
                        property: b'q',
                        relation: 1,
                        value,
                    });
            }
            sel
        }
    }

    /// Plain-data first-match hit (original row table indices).
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub struct ProbeHit {
        pub model: u32,
        pub chain: u32,
        pub residue: u32,
        pub atom: u32,
    }

    /// End-to-end first-match over a detached reader result.
    pub fn probe_first_rows(
        probe: &ProbeSelection,
        data: &BioStructureData,
        residue_flag: u8,
        atom_flag: u8,
    ) -> Result<Option<ProbeHit>, String> {
        let sel = probe.build();
        bio_first_rows(
            &sel,
            data.models(),
            data.chains(),
            data.residues(),
            data.atoms(),
            data.input_format(),
            residue_flag,
            atom_flag,
        )
        .map_err(|e: BioRowTraverseError| e.to_string())
        .map(|hit| {
            hit.map(|h| ProbeHit {
                model: h.model.value(),
                chain: h.cra.chain.value(),
                residue: h.cra.residue.value(),
                atom: h.cra.atom.value(),
            })
        })
    }

    /// End-to-end selected atom-ID iteration (collected for comparison).
    pub fn probe_selected_atom_ids(
        probe: &ProbeSelection,
        data: &BioStructureData,
        residue_flag: u8,
        atom_flag: u8,
    ) -> Result<Vec<u32>, String> {
        let sel = probe.build();
        let mut ids = Vec::new();
        for item in BioSelectedAtomIds::over_structure(
            &sel,
            data.models(),
            data.chains(),
            data.residues(),
            data.atoms(),
            data.input_format(),
            residue_flag,
            atom_flag,
        ) {
            match item {
                Ok(id) => ids.push(id.value()),
                Err(e) => return Err(e.to_string()),
            }
        }
        Ok(ids)
    }

    /// The document provenance of a detached reader result.
    pub fn probe_input_format(data: &BioStructureData) -> BioCoordinateFormat {
        data.input_format()
    }
}

pub use hierarchy::{
    AltLocRequest, BioAltLocGroupId, BioAssembly, BioAssemblyGenerator, BioAssemblyId,
    BioAssemblyOperator, BioAssemblySpecialKind, BioAtomId, BioAtomRow, BioCalcFlag, BioChainId,
    BioChainRow, BioCoordinateBlock, BioCoordinateFormat, BioCrystalCell, BioCrystalInfo,
    BioEntityDbRef, BioEntityId, BioEntityRow, BioModelId, BioModelRow, BioNcsOperator,
    BioNearestImage, BioResidueId, BioResidueRow, BioRowSpan, BioSiftsUnpResidue,
    BioStructureCopySource, BioStructureData, BioStructureError, BioStructureParts, BioTransform,
    ChainKind, EntityKind, PolymerKind, ResidueKind, altloc_matches, atom_name_logical_view,
    find_nearest_image, is_same_conformer, residue_name_logical_view, set_crystal_cell,
    set_crystal_fractional_transform, set_crystal_space_group_hm, set_crystal_z_pdb_if_nonempty,
    setup_cell_images,
};
pub use metadata::{
    BioBasicRefinementInfo, BioDiffractionInfo, BioExperimentInfo, BioExperimentalCrystalInfo,
    BioMetadata, BioRefinementInfo, BioRefinementRestraint, BioReflectionsInfo,
    BioSoftwareClassification, BioSoftwareItem, BioTlsGroup, BioTlsSelection,
};
pub use protein::{
    ProteinAtomIter, ProteinAtomRef, ProteinChainIter, ProteinChainRef, ProteinData,
    ProteinProjectionError, ProteinResidueIter, ProteinResidueRef, ProteinSelectionSummary,
    protein_atoms, protein_chain, protein_chains, protein_residues, protein_selection_summary,
    validate_protein_structure,
};
pub use relationships::{
    AtomAddress, BioAsu, BioCisPep, BioConnection, BioConnectionKind, BioModRes, ResidueAddress,
};
pub use residue::{
    ResidueCode, ResidueCodeParseError, ResidueIdentity, ResidueInfo, ResidueInfoKind,
    ResidueSequenceError, UNKNOWN_TABULATED_RESIDUE_INDEX, expand_one_letter,
    expand_one_letter_sequence, find_residue_info, find_residue_info_index, residue_code,
    residue_info, residue_info_checked,
};
pub use secondary_structure::{BioHelix, BioHelixClass, BioSheet, BioStrand};
pub use selection::{
    BioRowChainError, BioRowModelError, BioRowTraverseError, BioSelectionBlocks,
    BioSelectionCopyCause, BioSelectionCopyError, BioSelectionData, BioSelectionMatchError,
    BioSelectionParts, BioSelectionPartsView, CidListUtf8Error, GEMMI_EL_END, GEMMI_ELEMENT_NAMES,
    GEMMI_IS_METAL, GemmiElementMask, SelectionAtomInequality, SelectionCidList, SelectionList,
    SelectionSequenceId, copy_selection_blocks, gemmi_element_name, gemmi_find_element,
};
pub use source_ids::{
    AltLocLabel, AtomName, AtomSourceIds, ChainSourceIds, EntitySourceIds, PdbAtomSerial,
    PdbChainId, PdbSeqId, ResidueName, ResidueSourceIds,
};
pub use structure_metadata::BioStructureSourceState;

use cosmolkit_model::TopologyBlock;

/// Translate coordinate rows only; metadata and anisotropic tensors are unchanged.
pub fn translate_coordinates(coordinates: &mut BioCoordinateBlock, offset: [f64; 3]) {
    // Project-defined coordinate-only translation. No full Gemmi structure
    // transformation or metadata rewrite is claimed. One pass, O(n), no allocation.
    for position in coordinates.positions_mut() {
        for axis in 0..3 {
            position[axis] += offset[axis];
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BioError {
    Unsupported,
}

pub fn select_residues(topology: &TopologyBlock, query: &str) -> Result<TopologyBlock, BioError> {
    let _ = (topology, query);
    Err(BioError::Unsupported)
}
