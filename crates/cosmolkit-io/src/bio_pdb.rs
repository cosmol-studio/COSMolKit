//! Gemmi-aligned structural PDB input primitives.
//!
//! This is separate from `pdb.rs`, which reads PDB into a detached chemical
//! molecule rather than a `BioStructureData` hierarchy.

use std::collections::HashMap;
use std::hash::{Hash, Hasher};

use cosmolkit_bio::{
    AltLocLabel, AltLocRequest, AtomAddress, AtomName, AtomSourceIds, BioAssembly,
    BioAssemblyGenerator, BioAssemblyOperator, BioAssemblySpecialKind, BioAsu, BioAtomRow,
    BioCalcFlag, BioChainId, BioChainRow, BioCisPep, BioConnection, BioConnectionKind,
    BioCoordinateBlock, BioCoordinateFormat, BioCrystalCell, BioCrystalInfo, BioDiffractionInfo,
    BioEntityDbRef, BioEntityId, BioEntityRow, BioExperimentInfo, BioExperimentalCrystalInfo,
    BioHelix, BioMetadata, BioModRes, BioModelId, BioModelRow, BioNcsOperator, BioRefinementInfo,
    BioRefinementRestraint, BioResidueId, BioResidueRow, BioRowSpan, BioSheet, BioSiftsUnpResidue,
    BioSoftwareClassification, BioSoftwareItem, BioStrand, BioStructureData, BioStructureError,
    BioStructureParts, BioStructureSourceState, BioTlsGroup, BioTlsSelection, BioTransform,
    ChainKind, ChainSourceIds, EntityKind, EntitySourceIds, PdbAtomSerial, PdbChainId, PdbSeqId,
    ResidueAddress, ResidueInfoKind, ResidueName, ResidueSourceIds, find_nearest_image,
    find_residue_info, is_same_conformer, set_crystal_cell, set_crystal_fractional_transform,
    set_crystal_space_group_hm, set_crystal_z_pdb_if_nonempty, setup_cell_images,
};
use cosmolkit_types::Element;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
struct PdbInputOptions {
    // Gemmi's PdbReadOptions default is zero; the stream reader normalizes it
    // to the format's 120-byte maximum before reading.
    max_line_length: i32,
    // Gemmi source: third_party/gemmi/include/gemmi/model.hpp, PdbReadOptions.
    // Gemmi❗✔️:   bool ignore_ter = false; // ignores TER records completely
    ignore_ter: bool,
    // Gemmi❗✔️:   bool split_chain_on_ter = false;
    split_chain_on_ter: bool,
}

impl PdbInputOptions {
    fn effective_max_line_length(self) -> usize {
        // Gemmi source: third_party/gemmi/src/pdb.cpp, populate_structure_from_pdb_stream
        // Gemmi❗✔️:   if (options.max_line_length <= 0 || options.max_line_length > 120)
        // Gemmi❗✔️:     options.max_line_length = 120;
        if self.max_line_length <= 0 || self.max_line_length > 120 {
            120
        } else {
            self.max_line_length as usize
        }
    }
}

/// Source options for reading PDB text into the detached BIO hierarchy.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct BioPdbReadParams {
    /// Values `<= 0` or `> 120` are normalized to Gemmi's 120-byte maximum.
    pub max_line_length: i32,
    /// Record the first line containing a non-ASCII byte without rejecting it.
    pub check_non_ascii: bool,
    /// Ignore TER records completely.
    pub ignore_ter: bool,
    /// Split chain parts at TER without inferring residue entity kinds.
    pub split_chain_on_ter: bool,
    /// Retain raw remarks but skip their later metadata interpretation.
    pub skip_remarks: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PdbReaderStage {
    Stream,
    Record([u8; 4]),
    Finalization,
}

#[derive(Debug)]
enum PdbReaderCause {
    LineNumberOutsideSourceIntRange,
    InputFormat(PdbInputFormatError),
    AtomFields(PdbAtomFieldError),
    ModelTransition(PdbModelTransitionError),
    Remark(PdbRemarkError),
    Conect(PdbConectError),
    ModRes(PdbModResError),
    Helix(PdbHelixError),
    Sheet(PdbSheetError),
    Header(PdbHeaderError),
    DbRef(PdbDbRefError),
    Mtrix(PdbMtrixError),
    Crystal(PdbCrystalError),
    Finalization(PdbStructureFinalizationError),
}

impl std::fmt::Display for PdbReaderCause {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for PdbReaderCause {}

#[derive(Debug)]
struct PdbReaderError {
    line_number: Option<i32>,
    stage: PdbReaderStage,
    cause: PdbReaderCause,
}

impl std::fmt::Display for PdbReaderError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.line_number {
            Some(line_number) => write!(
                formatter,
                "PDB input failed at line {line_number} during {:?}: {}",
                self.stage, self.cause
            ),
            None => write!(
                formatter,
                "PDB input failed during {:?}: {}",
                self.stage, self.cause
            ),
        }
    }
}

impl std::error::Error for PdbReaderError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(&self.cause)
    }
}

/// Observable phase for a structured detached PDB reader failure.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BioPdbReadStage {
    Stream,
    Record,
    Finalization,
}

/// A source-context-preserving error from the detached BIO PDB reader.
#[derive(Debug)]
pub struct BioPdbReadError {
    inner: PdbReaderError,
}

impl BioPdbReadError {
    /// Physical one-based input line, when the failure occurred while handling a row.
    #[must_use]
    pub const fn line_number(&self) -> Option<i32> {
        self.inner.line_number
    }

    /// High-level reader phase that failed.
    #[must_use]
    pub const fn stage(&self) -> BioPdbReadStage {
        match self.inner.stage {
            PdbReaderStage::Stream => BioPdbReadStage::Stream,
            PdbReaderStage::Record(_) => BioPdbReadStage::Record,
            PdbReaderStage::Finalization => BioPdbReadStage::Finalization,
        }
    }

    /// Four-byte record discriminator for record-local failures.
    #[must_use]
    pub const fn record_tag(&self) -> Option<[u8; 4]> {
        match self.inner.stage {
            PdbReaderStage::Record(tag) => Some(tag),
            PdbReaderStage::Stream | PdbReaderStage::Finalization => None,
        }
    }
}

impl std::fmt::Display for BioPdbReadError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Display::fmt(&self.inner, formatter)
    }
}

impl std::error::Error for BioPdbReadError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(&self.inner.cause)
    }
}

struct PdbBioReader<'text, 'source> {
    cursor: PdbLineCursor<'text>,
    source: &'source str,
    options: PdbInputOptions,
    skip_remarks: bool,
    state: PdbReaderState,
    transition: PdbModelTransition,
    connection_records: Vec<Vec<u8>>,
}

impl<'text, 'source> PdbBioReader<'text, 'source> {
    fn new(
        text: &'text str,
        source: &'source str,
        options: PdbInputOptions,
        check_non_ascii: bool,
        skip_remarks: bool,
    ) -> Self {
        // Gemmi source: third_party/gemmi/src/pdb.cpp,
        // populate_structure_from_pdb_stream.
        // Gemmi❗✔️:   Transform matrix;
        // Gemmi❗✔️:   std::vector<std::string> conn_records;
        // Gemmi❗✔️:   Model *model = nullptr;
        // Gemmi❗✔️:   Chain *chain = nullptr;
        // Gemmi❗✔️:   Residue *resi = nullptr;
        // Gemmi❗✔️:   char line[122] = {0};
        // Gemmi❗✔️:   int line_num = 0;
        // Gemmi❗✔️:   bool after_ter = false;
        // Behavior review: parser-local state is split across the existing
        // typed reader, hierarchy transition and fixed-buffer cursor owners;
        // raw connection rows remain encounter-ordered until finalization.
        // The source string is borrowed for format diagnostics rather than
        // copied, and no second structural hierarchy is created.
        // Complexity review: initialization is constant-size plus the
        // existing fixed 122-byte cursor buffer; empty Vec allocation is
        // deferred and the input/source strings remain borrowed.
        Self {
            cursor: PdbLineCursor::new(text, options),
            source,
            options,
            skip_remarks,
            state: PdbReaderState::new(source, check_non_ascii, BioStructureSourceState::default()),
            transition: PdbModelTransition::default(),
            connection_records: Vec::new(),
        }
    }

    fn record_error(&self, record: [u8; 4], cause: PdbReaderCause) -> PdbReaderError {
        // Gemmi source: third_party/gemmi/src/pdb.cpp,
        // populate_structure_from_pdb_stream.
        // Gemmi❗✔️:   auto wrong = [&line_num](const std::string& msg) {
        // Gemmi❗✔️:     fail("Problem in line ", std::to_string(line_num), ": ", msg);
        // Gemmi❗✔️:   };
        // Behavior review: preserve the current one-based line and exact
        // record-stage context while transporting the owner error enum itself;
        // no cause is converted to a string or discarded. The wrapper's
        // Display is presentation only and its `cause` remains matchable.
        // Complexity review: copies the fixed record tag and moves one typed
        // cause; it performs no input scan or source-text allocation.
        PdbReaderError {
            line_number: Some(self.state.line_number),
            stage: PdbReaderStage::Record(record),
            cause,
        }
    }

    fn stream_error(&self, cause: PdbReaderCause) -> PdbReaderError {
        PdbReaderError {
            line_number: None,
            stage: PdbReaderStage::Stream,
            cause,
        }
    }

    fn finalization_error(&self, cause: PdbReaderCause) -> PdbReaderError {
        PdbReaderError {
            line_number: None,
            stage: PdbReaderStage::Finalization,
            cause,
        }
    }

    fn dispatch_hierarchy_record(&mut self, copied_len: usize) -> Result<bool, PdbReaderError> {
        let source_line_buffer = &self.cursor.line_buffer;
        let line = &source_line_buffer[..copied_len];
        let mut record = [0; 4];
        let record_len = line.len().min(record.len());
        record[..record_len].copy_from_slice(&line[..record_len]);

        // Gemmi source: third_party/gemmi/src/pdb.cpp,
        // populate_structure_from_pdb_stream.
        // Gemmi❗✔️:     if (is_record_type4(line, "ATOM") || is_record_type4(line, "HETATM")) {
        // Gemmi❗✔️:       if (len < 55)
        // Gemmi❗✔️:         wrong("The line is too short to be correct:\n" + std::string(line));
        // Gemmi❗✔️:       std::string chain_name = read_string(line+20, 2);
        // Gemmi❗✔️:       ResidueId rid = read_res_id(line+22, line+17);
        // Gemmi❗✔️:       if (len > 72)
        // Gemmi❗✔️:         rid.segment = read_string(line+72, 4);
        // Gemmi❗✔️:       Atom atom;
        // Gemmi❗✔️:       atom.serial = read_serial(line+6);
        // Gemmi❗✔️:       atom.name = read_string(line+12, 4);
        // Gemmi❗✔️:       atom.altloc = read_altloc(line[16]);
        // Gemmi❗✔️:       atom.pos.x = read_double(line+30, 8);
        // Gemmi❗✔️:       atom.pos.y = read_double(line+38, 8);
        // Gemmi❗✔️:       atom.pos.z = read_double(line+46, 8);
        // Gemmi❗✔️:       if (len > 58)
        // Gemmi❗✔️:         atom.occ = (float) read_double(line+54, 6);
        // Gemmi❗✔️:       if (len > 64)
        // Gemmi❗✔️:         atom.b_iso = (float) read_double(line+60, 6);
        // Gemmi❗✔️:       if (len > 76 && (std::isalpha(line[76]) || std::isalpha(line[77])))
        // Gemmi❗✔️:         atom.element = Element(line + 76);
        // Gemmi❗✔️:       else
        // Gemmi❗✔️:         atom.element = infer_element_from_padded_name(line+12);
        // Gemmi❗✔️:       atom.charge = (len > 78 ? read_charge(line[78], line[79]) : 0);
        // Gemmi❗✔️:       resi->atoms.emplace_back(atom);
        // Gemmi❗✔️:     } else if (is_record_type4(line, "ANISOU")) {
        // Gemmi❗✔️:       if (!model || !chain || !resi || resi->atoms.empty())
        // Gemmi❗✔️:         wrong("ANISOU record not directly after ATOM/HETATM.");
        // Gemmi❗✔️:       // We assume that ANISOU refers to the last atom.
        // Gemmi❗✔️:       // Can it not be the case?
        // Gemmi❗✔️:       Atom &atom = resi->atoms.back();
        // Gemmi❗✔️:       if (atom.aniso.u11 != 0.)
        // Gemmi❗✔️:         wrong("Duplicated ANISOU record or not directly after ATOM/HETATM.");
        // Gemmi❗✔️:       atom.aniso.u11 = read_int(line+28, 7) * 1e-4f;
        // Gemmi❗✔️:       atom.aniso.u22 = read_int(line+35, 7) * 1e-4f;
        // Gemmi❗✔️:       atom.aniso.u33 = read_int(line+42, 7) * 1e-4f;
        // Gemmi❗✔️:       atom.aniso.u12 = read_int(line+49, 7) * 1e-4f;
        // Gemmi❗✔️:       atom.aniso.u13 = read_int(line+56, 7) * 1e-4f;
        // Gemmi❗✔️:       atom.aniso.u23 = read_int(line+63, 7) * 1e-4f;
        // Gemmi❗✔️:     } else if (is_record_type3(line, "TER") && !options.ignore_ter) {
        // Gemmi❗✔️:       if (!chain || st.ter_status == 'e')
        // Gemmi❗✔️:         continue;
        // Gemmi❗✔️:       st.ter_status = 'y';
        // Gemmi❗✔️:       if (options.split_chain_on_ter) {
        // Gemmi❗✔️:         chain = nullptr;
        // Gemmi❗✔️:         // split_chain_on_ter is used for AMBER files that can have TER records
        // Gemmi❗✔️:         // in various places. So in such case TER doesn't imply entity_type.
        // Gemmi❗✔️:         continue;
        // Gemmi❗✔️:       }
        // Gemmi❗✔️:       if (after_ter) {
        // Gemmi❗✔️:         st.ter_status = 'e';
        // Gemmi❗✔️:         continue;
        // Gemmi❗✔️:       }
        // Gemmi❗✔️:       for (Residue& res : chain->residues) {
        // Gemmi❗✔️:         res.entity_type = EntityType::Polymer;
        // Gemmi❗✔️:         if GEMMI_UNLIKELY(res.is_water())
        // Gemmi❗✔️:           st.ter_status = 'e';
        // Gemmi❗✔️:       }
        // Gemmi❗✔️:       after_ter = true;
        // Gemmi❗✔️:     } else if (is_record_type4(line, "MODEL")) {
        // Gemmi❗✔️:       if (model && chain)
        // Gemmi❗✔️:         wrong("MODEL without ENDMDL?");
        // Gemmi❗✔️:       int num = read_int(line+6, 8);
        // Gemmi❗✔️:       model = &st.find_or_add_model(num);
        // Gemmi❗✔️:       if (!model->chains.empty())
        // Gemmi❗✔️:         wrong("duplicate MODEL number: " + std::to_string(num));
        // Gemmi❗✔️:       chain = nullptr;
        // Gemmi❗✔️:     } else if (is_record_type4(line, "ENDMDL")) {
        // Gemmi❗✔️:       model = nullptr;
        // Gemmi❗✔️:       chain = nullptr;
        // Gemmi❗✔️:     } else if (is_record_type3(line, "END")) {
        // Gemmi❗✔️:       break;
        // Gemmi❗✔️:     }
        // Behavior review: the hierarchy subset retains source-relative branch
        // precedence (ATOM/HETATM, ANISOU, TER, MODEL, ENDMDL, END); metadata
        // and crystal arms are inserted between these disjoint predicates by
        // the later full dispatcher. Atom fields use the existing fixed-width
        // decoder, and the existing typed chain/residue owners carry all
        // representable source fields. Errors retain record and line context.
        // Complexity review: record tests are constant-count fixed-width
        // comparisons; atom decode is bounded by the PDB line/field widths,
        // chain transition uses Gemmi's first-match scan, and residue lookup
        // uses the existing expected-O(1) index. No line or hierarchy clone.
        if gemmi_record_type4(line, 0, *b"ATOM") || gemmi_record_type4(line, 0, *b"HETA") {
            if line.len() < 55 {
                return Err(self.record_error(
                    record,
                    PdbReaderCause::AtomFields(PdbAtomFieldError::TooShort { len: line.len() }),
                ));
            }
            let raw_chain: [u8; 2] = line[20..22]
                .try_into()
                .expect("the preceding source-length gate checked the atom record");
            let chain_bytes = read_string(&raw_chain);
            let chain_name = PdbChainId::from_ascii(&chain_bytes).ok_or_else(|| {
                self.record_error(
                    record,
                    PdbReaderCause::AtomFields(PdbAtomFieldError::ChainNameNotRepresentable {
                        bytes: raw_chain,
                    }),
                )
            })?;
            let sequence_field: &[u8; 5] = line[22..27]
                .try_into()
                .expect("the preceding source-length gate checked the atom record");
            let residue_name_field: &[u8; 3] = line[17..20]
                .try_into()
                .expect("the preceding source-length gate checked the atom record");
            let residue = read_res_id(sequence_field, residue_name_field).map_err(|cause| {
                self.record_error(
                    record,
                    PdbReaderCause::AtomFields(PdbAtomFieldError::ResidueKey(cause)),
                )
            })?;
            let residue = read_pdb_residue_segment(residue, line).map_err(|cause| {
                self.record_error(
                    record,
                    PdbReaderCause::AtomFields(PdbAtomFieldError::ResidueKey(cause)),
                )
            })?;
            let fields = decode_pdb_atom_fields(line, source_line_buffer)
                .map_err(|cause| self.record_error(record, PdbReaderCause::AtomFields(cause)))?;
            self.transition
                .add_atom_record(chain_name, residue, line[0], fields)
                .map_err(|cause| {
                    self.record_error(record, PdbReaderCause::ModelTransition(cause))
                })?;
            return Ok(true);
        }
        if gemmi_record_type4(line, 0, *b"ANIS") {
            self.transition
                .anisou_record(source_line_buffer)
                .map_err(|cause| {
                    self.record_error(record, PdbReaderCause::ModelTransition(cause))
                })?;
            return Ok(true);
        }
        if gemmi_record_type3(line, 0, *b"TER\0") && !self.options.ignore_ter {
            self.transition
                .ter_record(line, self.options, &mut self.state.source_state)
                .map_err(|cause| {
                    self.record_error(record, PdbReaderCause::ModelTransition(cause))
                })?;
            return Ok(true);
        }
        if gemmi_record_type4(line, 0, *b"MODE") {
            self.transition
                .model_record(source_line_buffer)
                .map_err(|cause| {
                    self.record_error(record, PdbReaderCause::ModelTransition(cause))
                })?;
            return Ok(true);
        }
        if gemmi_record_type4(line, 0, *b"ENDM") {
            self.transition.end_model_record();
            return Ok(true);
        }
        if gemmi_record_type3(line, 0, *b"END\0") {
            self.transition.end_record();
            return Ok(true);
        }
        Ok(false)
    }

    fn dispatch_metadata_record(
        &mut self,
        copied_len: usize,
    ) -> Result<Option<PdbDbRefAction>, PdbReaderError> {
        let source_line_buffer = &self.cursor.line_buffer;
        let line = &source_line_buffer[..copied_len];
        let mut record = [0; 4];
        let record_len = line.len().min(record.len());
        record[..record_len].copy_from_slice(&line[..record_len]);

        // Gemmi source: third_party/gemmi/src/pdb.cpp,
        // populate_structure_from_pdb_stream. These verbatim branch entries
        // preserve the pinned record order; each delegated owner carries the
        // complete corresponding field-processing anchor in its own body.
        // Gemmi❗✔️:     } else if (is_record_type4(line, "REMARK")) {
        // Gemmi❗✔️:     } else if (is_record_type4(line, "CONECT")) {
        // Gemmi❗✔️:     } else if (is_record_type4(line, "SEQRES")) {
        // Gemmi❗✔️:     } else if (is_record_type4(line, "HELIX")) {
        // Gemmi❗✔️:     } else if (is_record_type4(line, "SHEET")) {
        // Gemmi❗✔️:     } else if (is_record_type4(line, "SSBOND") ||
        // Gemmi❗✔️:                is_record_type4(line, "LINK") ||
        // Gemmi❗✔️:                is_record_type4(line, "CISPEP")) {
        // Gemmi❗✔️:       conn_records.emplace_back(line);
        // Gemmi❗✔️:     } else if (is_record_type4(line, "MODRES")) {
        // Gemmi❗✔️:     } else if (is_record_type4(line, "HETNAM")) {
        // Gemmi❗✔️:     } else if (is_record_type4(line, "DBREF")) { // DBREF or DBREF1 or DBREF2
        // Gemmi❗✔️:       else if (ent.dbrefs.empty()) // DBREF2 without DBREF1?
        // Gemmi❗✔️:         break;
        // Gemmi❗✔️:     } else if (is_record_type4(line, "HEADER")) {
        // Gemmi❗✔️:     } else if (is_record_type4(line, "TITLE")) {
        // Gemmi❗✔️:     } else if (is_record_type4(line, "KEYWDS")) {
        // Gemmi❗✔️:     } else if (is_record_type4(line, "EXPDTA")) {
        // Gemmi❗✔️:     } else if (is_record_type4(line, "AUTHOR") && len > 10) {
        // Gemmi❗✔️:     }
        // Behavior review: each raw PDB row is routed to its existing
        // source-backed state owner; record matching and source order follow
        // the pinned chain. A matched row returns Some(Continue), an orphan
        // DBREF2 returns Some(Stop) so the complete reader can stop scanning
        // but still finalize, and None leaves the row for a later dispatcher.
        // REMARK capture is deliberately unconditional here; interpretation
        // remains deferred to the finalizer's skip_remarks branch. Connection
        // records retain encounter order and exact visible source bytes.
        // Complexity review: a fixed sequence of record checks invokes at
        // most one bounded existing parser; only deferred connection rows
        // allocate a source-line copy, matching Gemmi's per-record string
        // construction. No hierarchy scan or parser-state clone is added.
        if self
            .state
            .remark_record(source_line_buffer, copied_len)
            .map_err(|cause| self.record_error(record, PdbReaderCause::Remark(cause)))?
        {
            return Ok(Some(PdbDbRefAction::Continue));
        }
        if self
            .state
            .conect_record(source_line_buffer, copied_len)
            .map_err(|cause| self.record_error(record, PdbReaderCause::Conect(cause)))?
        {
            return Ok(Some(PdbDbRefAction::Continue));
        }
        if self.state.seqres_record(line) {
            return Ok(Some(PdbDbRefAction::Continue));
        }
        if self
            .state
            .helix_record(source_line_buffer, copied_len)
            .map_err(|cause| self.record_error(record, PdbReaderCause::Helix(cause)))?
        {
            return Ok(Some(PdbDbRefAction::Continue));
        }
        if self
            .state
            .sheet_record(source_line_buffer, copied_len)
            .map_err(|cause| self.record_error(record, PdbReaderCause::Sheet(cause)))?
        {
            return Ok(Some(PdbDbRefAction::Continue));
        }
        if gemmi_record_type4(line, 0, *b"SSBO")
            || gemmi_record_type4(line, 0, *b"LINK")
            || gemmi_record_type4(line, 0, *b"CISP")
        {
            self.connection_records.push(line.to_vec());
            return Ok(Some(PdbDbRefAction::Continue));
        }
        if self
            .state
            .modres_record(source_line_buffer, copied_len)
            .map_err(|cause| self.record_error(record, PdbReaderCause::ModRes(cause)))?
        {
            return Ok(Some(PdbDbRefAction::Continue));
        }
        if self.state.hetnam_record(source_line_buffer, copied_len) {
            return Ok(Some(PdbDbRefAction::Continue));
        }

        let dbref_action = self
            .state
            .dbref_record(source_line_buffer)
            .map_err(|cause| self.record_error(record, PdbReaderCause::DbRef(cause)))?;
        if dbref_action != PdbDbRefAction::NotDbRef {
            return Ok(Some(dbref_action));
        }
        if self
            .state
            .header_record(source_line_buffer, copied_len)
            .map_err(|cause| self.record_error(record, PdbReaderCause::Header(cause)))?
        {
            return Ok(Some(PdbDbRefAction::Continue));
        }
        Ok(None)
    }

    fn dispatch_crystal_format_record(
        &mut self,
        copied_len: usize,
    ) -> Result<bool, PdbReaderError> {
        let source_line_buffer = &self.cursor.line_buffer;
        let line = &source_line_buffer[..copied_len];
        let mut record = [0; 4];
        let record_len = line.len().min(record.len());
        record[..record_len].copy_from_slice(&line[..record_len]);

        // Gemmi source: third_party/gemmi/src/pdb.cpp,
        // populate_structure_from_pdb_stream. The record-order entries below
        // are copied from the source dispatcher; the invoked row helpers keep
        // their complete field and state-transition anchors in their owners.
        // Gemmi❗✔️:     } else if (is_record_type4(line, "SCALEn")) {
        // Gemmi❗✔️:       if (read_matrix(matrix, line, len) == 3) {
        // Gemmi❗✔️:         st.cell.set_matrices_from_fract(matrix);
        // Gemmi❗✔️:         matrix.set_identity();
        // Gemmi❗✔️:       }
        // Gemmi❗✔️:     } else if (is_record_type4(line, "ORIGX")) {
        // Gemmi❗✔️:       st.has_origx = true;
        // Gemmi❗✔️:       read_matrix(st.origx, line, len);
        // Gemmi❗✔️:     } else if (is_record_type4(line, "CRYST1")) {
        // Gemmi❗✔️:       if (len > 54)
        // Gemmi❗✔️:         st.cell.set(read_double(line+6, 9),
        // Gemmi❗✔️:                     read_double(line+15, 9),
        // Gemmi❗✔️:                     read_double(line+24, 9),
        // Gemmi❗✔️:                     read_double(line+33, 7),
        // Gemmi❗✔️:                     read_double(line+40, 7),
        // Gemmi❗✔️:                     read_double(line+47, 7));
        // Gemmi❗✔️:       if (len > 56)
        // Gemmi❗✔️:         st.spacegroup_hm = read_string(line+55, 11);
        // Gemmi❗✔️:       if (len > 67) {
        // Gemmi❗✔️:         std::string z = read_string(line+66, 4);
        // Gemmi❗✔️:         if (!z.empty())
        // Gemmi❗✔️:           st.info["_cell.Z_PDB"] = z;
        // Gemmi❗✔️:       }
        // Gemmi❗✔️:     } else if (is_record_type4(line, "MTRIXn")) {
        // Gemmi❗✔️:       if (read_matrix(matrix, line, len) == 3) {
        // Gemmi❗✔️:         std::string id = read_string(line+7, 3);
        // Gemmi❗✔️:         if (matrix.is_identity()) {
        // Gemmi❗✔️:           st.info["_struct_ncs_oper.id"] = id;
        // Gemmi❗✔️:         } else {
        // Gemmi❗✔️:           bool given = len > 59 && line[59] == '1';
        // Gemmi❗✔️:           st.ncs.push_back({id, given, matrix});
        // Gemmi❗✔️:           matrix.set_identity();
        // Gemmi❗✔️:         }
        // Gemmi❗✔️:       }
        // Gemmi❗✔️:     } else if (is_record_type4(line, "data")) {
        // Gemmi❗✔️:       if (line[4] == '_' && !model)
        // Gemmi❗✔️:         fail("Incorrect file format (perhaps it is cif not pdb?): " + source);
        // Gemmi❗✔️:     } else if (is_record_type4(line, "{\"da")) {
        // Gemmi❗✔️:       if (ialpha3_id(line+4) == ialpha3_id("ta_") && !model)
        // Gemmi❗✔️:         fail("Incorrect file format (perhaps it is mmJSON not pdb?): " + source);
        // Behavior review: route matrix rows in the exact source order through
        // the existing shared SCALE/MTRIX staging value and source-specific
        // ORIGX/CRYST1 owners. Only the model-active state suppresses the two
        // input-format errors. These checks inspect the reusable 122-byte
        // source buffer (including suffix bytes after a C-string NUL), matching
        // the source's direct indexed reads; a nonmatching record is an
        // intentional no-op for a later dispatcher/default branch.
        // Complexity review: four fixed record checks plus one constant-size
        // format discriminator; matrix parsing and crystal updates are bounded
        // fixed-field work in their existing owners. No row scans, clones, or
        // allocations are added on the ordinary unrecognized-record path.
        if gemmi_record_type4(line, 0, *b"SCAL") {
            self.state.scale_record(line);
            return Ok(true);
        }
        if gemmi_record_type4(line, 0, *b"ORIG") {
            self.state.origx_record(line);
            return Ok(true);
        }
        if gemmi_record_type4(line, 0, *b"CRYS") {
            self.state
                .cryst1_record(line)
                .map_err(|cause| self.record_error(record, PdbReaderCause::Crystal(cause)))?;
            return Ok(true);
        }
        if gemmi_record_type4(line, 0, *b"MTRI") {
            self.state
                .mtrix_record(line)
                .map_err(|cause| self.record_error(record, PdbReaderCause::Mtrix(cause)))?;
            return Ok(true);
        }

        if let Some(cause) =
            self.transition
                .wrong_input_format(&self.state, source_line_buffer, self.source)
        {
            return Err(self.record_error(record, PdbReaderCause::InputFormat(cause)));
        }
        Ok(false)
    }

    fn read_structure(mut self) -> Result<BioStructureData, PdbReaderError> {
        // Gemmi source: third_party/gemmi/src/pdb.cpp,
        // populate_structure_from_pdb_stream.
        // Gemmi❗✔️:   while (size_t len = line_reader.copy_line(line, options.max_line_length+1)) {
        // Gemmi❗✔️:     ++line_num;
        // Gemmi❗✔️:     if (options.check_non_ascii && st.non_ascii_line == 0)
        // Gemmi❗✔️:       for (size_t i = 0; i < len; ++i)
        // Gemmi❗✔️:         if (static_cast<unsigned char>(line[i]) >= 0x80) {
        // Gemmi❗✔️:           st.non_ascii_line = line_num;
        // Gemmi❗✔️:           break;
        // Gemmi❗✔️:         }
        // Gemmi❗✔️:   }
        // Gemmi❗✔️:   } else if (is_record_type3(line, "END")) {
        // Gemmi❗✔️:     break;
        // Gemmi❗✔️:   } else if (is_record_type4(line, "DBREF")) {
        // Gemmi❗✔️:     else if (ent.dbrefs.empty()) // DBREF2 without DBREF1?
        // Gemmi❗✔️:       break;
        // Behavior review: the already source-backed record owners run once
        // per copied physical line in hierarchy, metadata, then crystal/format
        // order. END and orphan DBREF2 stop this loop but still reach the one
        // shared finalizer below. The cursor's normalized fixed-buffer and
        // suffix behavior is preserved; line counting and first-non-ASCII
        // recording occur before dispatch. Unknown rows do not mutate state.
        // Complexity review: one O(line length) bounded copy/scan per physical
        // line plus a fixed sequence of existing record predicates and owners;
        // the same fixed cursor is reused and the driver adds no per-line
        // allocation or hierarchy scan.
        while let Some(copied_len) = self.cursor.copy_line().map(|line| line.len()) {
            let source_line = &self.cursor.line_buffer[..copied_len];
            self.state
                .record_line_checked(source_line)
                .map_err(|cause| self.stream_error(cause))?;

            if self.dispatch_hierarchy_record(copied_len)? {
                if self.transition.stopped_at_end_record() {
                    break;
                }
                continue;
            }

            if let Some(action) = self.dispatch_metadata_record(copied_len)? {
                if action == PdbDbRefAction::Stop {
                    break;
                }
                continue;
            }

            self.dispatch_crystal_format_record(copied_len)?;
        }

        // Gemmi❗✔️:   if (st.models.empty())
        // Gemmi❗✔️:     st.models.emplace_back(1);
        // Gemmi❗✔️:   st.setup_cell_images();
        // Gemmi❗✔️:   process_conn(st, conn_records);
        // Gemmi❗✔️:   if (!options.skip_remarks)
        // Gemmi❗✔️:     read_metadata_from_remarks(st);
        // Gemmi❗✔️:   restore_full_ccd_codes(st);
        // The finalizer owns this complete source-ordered closure and is called
        // exactly once after EOF, END, or DBREF2's source break. Preserve the
        // nested finalization cause and distinguish it from record/stream
        // failures without formatting it into a string.
        let Self {
            state,
            transition,
            connection_records,
            skip_remarks,
            ..
        } = self;
        finalize_pdb_structure(state, transition, &connection_records, skip_remarks).map_err(
            |cause| PdbReaderError {
                line_number: None,
                stage: PdbReaderStage::Finalization,
                cause: PdbReaderCause::Finalization(cause),
            },
        )
    }
}

/// Read PDB text as a detached canonical BIO structure.
pub fn read_pdb_bio_structure(
    text: &str,
    source_name: &str,
    params: &BioPdbReadParams,
) -> Result<BioStructureData, BioPdbReadError> {
    // Gemmi source: third_party/gemmi/include/gemmi/pdb.hpp::read_pdb_string.
    // Gemmi❗✔️: inline Structure read_pdb_string(const std::string& str,
    // Gemmi❗✔️:                                  const std::string& name,
    // Gemmi❗✔️:                                  PdbReadOptions options={}) {
    // Gemmi❗✔️:   return read_pdb_from_memory(str.c_str(), str.length(), name, options);
    // Gemmi❗✔️: }
    // Gemmi❗✔️: inline Structure read_pdb_from_memory(const char* data, size_t size,
    // Gemmi❗✔️:                                       const std::string& name,
    // Gemmi❗✔️:                                       PdbReadOptions options={}) {
    // Gemmi❗✔️:   MemoryStream stream{data, size};
    // Gemmi❗✔️:   return read_pdb_from_stream(stream, name, options);
    // Gemmi❗✔️: }
    // Gemmi❗✔️: inline Structure read_pdb_from_stream(AnyStream& line_reader,
    // Gemmi❗✔️:                                       const std::string& source,
    // Gemmi❗✔️:                                       PdbReadOptions options) {
    // Gemmi❗✔️:   gemmi::Structure st;
    // Gemmi❗✔️:   populate_structure_from_pdb_stream(line_reader, source, st, options);
    // Gemmi❗✔️:   return st;
    // Gemmi❗✔️: };
    // Behavior review: the public domain entry forwards the original text,
    // source name and all five source options into the sole fixed-buffer
    // reader/finalizer path; it does not expose the parser, staging hierarchy,
    // or mutable state. Nested typed causes and line/phase context remain in
    // BioPdbReadError rather than being flattened into a message.
    // Complexity review: constructing the five-field private options value
    // and wrapping a returned error are O(1); text, structure rows and causes
    // are not copied by this forwarding boundary.
    let options = PdbInputOptions {
        max_line_length: params.max_line_length,
        ignore_ter: params.ignore_ter,
        split_chain_on_ter: params.split_chain_on_ter,
    };
    PdbBioReader::new(
        text,
        source_name,
        options,
        params.check_non_ascii,
        params.skip_remarks,
    )
    .read_structure()
    .map_err(|inner| BioPdbReadError { inner })
}

#[derive(Debug, Clone, PartialEq)]
struct PdbReaderState {
    input_format: BioCoordinateFormat,
    source_state: BioStructureSourceState,
    metadata: BioMetadata,
    assemblies: Vec<BioAssembly>,
    ncs_operators: Vec<BioNcsOperator>,
    remark350_matrix: BioTransform,
    crystal_info: BioCrystalInfo,
    matrix: BioTransform,
    helices: Vec<BioHelix>,
    sheets: Vec<BioSheet>,
    mod_residues: Vec<BioModRes>,
    shortened_ccd_codes: Vec<PdbCcdAlias>,
    entities: PdbEntityState,
    remark3_continuation: Option<PdbRemark3Continuation>,
    remark200_continuation: Option<usize>,
    line_number: i32,
    check_non_ascii: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct PdbRemark3Continuation {
    refinement_index: usize,
    tls_group_index: usize,
    selection_index: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct PdbCcdAlias {
    pub(super) full_code: Vec<u8>,
    pub(super) short_code: Vec<u8>,
}

impl PdbReaderState {
    fn new(source: &str, check_non_ascii: bool, mut source_state: BioStructureSourceState) -> Self {
        // Gemmi source: third_party/gemmi/src/pdb.cpp,
        // populate_structure_from_pdb_stream.
        // Gemmi❗✔️:   st.input_format = CoorFormat::Pdb;
        // Gemmi❗✔️:   st.name = path_basename(source, {".gz", ".pdb"});
        // Gemmi❗✔️:   int line_num = 0;
        // Gemmi❗✔️:   Transform matrix;
        // Gemmi❗✔️:   UnitCell cell;
        // Gemmi❗✔️:   std::string spacegroup_hm;
        // Gemmi❗✔️:   double a = 1.0, b = 1.0, c = 1.0;
        // Gemmi❗✔️:   double alpha = 90.0, beta = 90.0, gamma = 90.0;
        // Gemmi❗✔️:   double volume = 1.0;
        // Gemmi❗✔️:   double ar = 1.0, br = 1.0, cr = 1.0;
        // Gemmi❗✔️:   double cos_alphar = 0.0, cos_betar = 0.0, cos_gammar = 0.0;
        // Gemmi❗✔️:   bool explicit_matrices = false;
        // Gemmi❗✔️:   short cs_count = 0;
        // This state owns canonical format/name and line diagnostics; parser-
        // local hierarchy pointers and `after_ter` belong to later stages.
        // The one staged matrix and canonical crystal value below preserve
        // Gemmi's default transform/unit-cell state for subsequent records.
        source_state.name = gemmi_pdb_path_basename(source);
        Self {
            input_format: BioCoordinateFormat::Pdb,
            source_state,
            metadata: BioMetadata::default(),
            assemblies: Vec::new(),
            ncs_operators: Vec::new(),
            remark350_matrix: BioTransform::identity(),
            crystal_info: BioCrystalInfo::new(
                BioCrystalCell::default(),
                None,
                None,
                BioTransform::identity(),
                BioTransform::identity(),
                false,
                0,
                Vec::new(),
            ),
            matrix: BioTransform::identity(),
            helices: Vec::new(),
            sheets: Vec::new(),
            mod_residues: Vec::new(),
            shortened_ccd_codes: Vec::new(),
            entities: PdbEntityState::default(),
            remark3_continuation: None,
            remark200_continuation: None,
            line_number: 0,
            check_non_ascii,
        }
    }

    fn scale_record(&mut self, line: &[u8]) {
        // Gemmi source: third_party/gemmi/src/pdb.cpp,
        // populate_structure_from_pdb_stream.
        // Gemmi❗✔️:     } else if (is_record_type4(line, "SCALEn")) {
        // Gemmi❗✔️:       if (read_matrix(matrix, line, len) == 3) {
        // Gemmi❗✔️:         st.cell.set_matrices_from_fract(matrix);
        // Gemmi❗✔️:         matrix.set_identity();
        // Gemmi❗✔️:       }
        // Behavior review: the caller owns record-prefix dispatch. Each row
        // mutates the single parser scratch transform; only row 3 installs it
        // through the canonical BIO transition, then resets the same scratch
        // transform for subsequent SCALE/MTRIX records.
        // Complexity review: one fixed-width row parse, constant-size matrix
        // transition and identity reset; no allocations or repeated scans.
        if read_matrix(&mut self.matrix, line) == 3 {
            set_crystal_fractional_transform(&mut self.crystal_info, self.matrix);
            self.matrix = BioTransform::identity();
        }
    }

    fn origx_record(&mut self, line: &[u8]) {
        // Gemmi source: third_party/gemmi/src/pdb.cpp,
        // populate_structure_from_pdb_stream.
        // Gemmi❗✔️:     } else if (is_record_type4(line, "ORIGX")) {
        // Gemmi❗✔️:       st.has_origx = true;
        // Gemmi❗✔️:       read_matrix(st.origx, line, len);
        // Behavior review: record dispatch is owned by the caller. Presence
        // is set before row decoding, including short or invalid-row records;
        // read_matrix then updates only a valid complete row and leaves other
        // rows and the default identity/zero components unchanged.
        // Complexity review: one presence assignment and one fixed-width
        // matrix-row parse with constant-size state; no allocation or scan.
        self.source_state.has_origx = true;
        read_matrix(&mut self.source_state.origx, line);
    }

    fn mtrix_record(&mut self, line: &[u8]) -> Result<(), PdbMtrixError> {
        // Gemmi source: third_party/gemmi/src/pdb.cpp,
        // populate_structure_from_pdb_stream.
        // Gemmi❗✔️:     } else if (is_record_type4(line, "MTRIXn")) {
        // Gemmi❗✔️:       if (read_matrix(matrix, line, len) == 3) {
        // Gemmi❗✔️:         std::string id = read_string(line+7, 3);
        // Gemmi❗✔️:         if (matrix.is_identity()) {
        // Gemmi❗✔️:           // store only ID that will be used when writing to file
        // Gemmi❗✔️:           st.info["_struct_ncs_oper.id"] = id;
        // Gemmi❗✔️:         } else {
        // Gemmi❗✔️:           bool given = len > 59 && line[59] == '1';
        // Gemmi❗✔️:           st.ncs.push_back({id, given, matrix});
        // Gemmi❗✔️:           matrix.set_identity();
        // Gemmi❗✔️:         }
        // Gemmi❗✔️:       }
        // Gemmi❗✔️:     }
        // Gemmi source: third_party/gemmi/include/gemmi/math.hpp,
        // Mat33::is_identity.
        // Gemmi❗✔️:   bool is_identity() const {
        // Gemmi❗✔️:     return a[0][0] == 1 && a[0][1] == 0 && a[0][2] == 0 &&
        // Gemmi❗✔️:            a[1][0] == 0 && a[1][1] == 1 && a[1][2] == 0 &&
        // Gemmi❗✔️:            a[2][0] == 0 && a[2][1] == 0 && a[2][2] == 1;
        // Gemmi❗✔️:   }
        // Gemmi source: third_party/gemmi/include/gemmi/math.hpp,
        // Transform::is_identity.
        // Gemmi❗✔️:   bool is_identity() const {
        // Gemmi❗✔️:     return mat.is_identity() && vec.x == 0. && vec.y == 0. && vec.z == 0.;
        // Gemmi❗✔️:   }
        // Behavior review: only suffix-3 rows finalize an operation; its own
        // ID and length-gated column 60 determine the result. Exact source
        // floating equality classifies identity (including signed zero); an
        // identity operation overwrites only the source info key, while a
        // nonidentity operation appends in source order then resets scratch.
        // Complexity review: each row uses the existing fixed-width parser;
        // completion performs a fixed 12 equality checks and one map insertion
        // or amortized vector append, with no scan over existing operators.
        if read_matrix(&mut self.matrix, line) != 3 {
            return Ok(());
        }

        let id_bytes = remark3_fixed_field(line, 7, 3);
        let id = String::from_utf8(id_bytes.to_vec()).map_err(|error| {
            PdbMtrixError::TextFieldNotUtf8 {
                field_offset: 7,
                valid_up_to: error.utf8_error().valid_up_to(),
            }
        })?;

        let matrix = self.matrix.matrix();
        let translation = self.matrix.translation();
        let is_identity = matrix[0][0] == 1.0
            && matrix[0][1] == 0.0
            && matrix[0][2] == 0.0
            && matrix[1][0] == 0.0
            && matrix[1][1] == 1.0
            && matrix[1][2] == 0.0
            && matrix[2][0] == 0.0
            && matrix[2][1] == 0.0
            && matrix[2][2] == 1.0
            && translation[0] == 0.0
            && translation[1] == 0.0
            && translation[2] == 0.0;

        if is_identity {
            self.source_state
                .info
                .insert("_struct_ncs_oper.id".to_owned(), id);
        } else {
            let given = line.len() > 59 && line[59] == b'1';
            self.ncs_operators
                .push(BioNcsOperator::new(id, given, self.matrix));
            self.matrix = BioTransform::identity();
        }

        Ok(())
    }

    fn cryst1_record(&mut self, line: &[u8]) -> Result<(), PdbCrystalError> {
        // Gemmi source: third_party/gemmi/src/pdb.cpp,
        // populate_structure_from_pdb_stream.
        // Gemmi❗✔️:     } else if (is_record_type4(line, "CRYST1")) {
        // Gemmi❗✔️:       if (len > 54)
        // Gemmi❗✔️:         st.cell.set(read_double(line+6, 9),
        // Gemmi❗✔️:                     read_double(line+15, 9),
        // Gemmi❗✔️:                     read_double(line+24, 9),
        // Gemmi❗✔️:                     read_double(line+33, 7),
        // Gemmi❗✔️:                     read_double(line+40, 7),
        // Gemmi❗✔️:                     read_double(line+47, 7));
        // Gemmi❗✔️:       if (len > 56)
        // Gemmi❗✔️:         st.spacegroup_hm = read_string(line+55, 11);
        // Gemmi❗✔️:       if (len > 67) {
        // Gemmi❗✔️:         std::string z = read_string(line+66, 4);
        // Gemmi❗✔️:         if (!z.empty())
        // Gemmi❗✔️:           st.info["_cell.Z_PDB"] = z;
        // Gemmi❗✔️:       }
        // Behavior review: the exact strict length gates and cell -> group ->
        // Z order are retained. Cell calculation errors propagate before later
        // assignments; group replaces even with decoded empty text; empty Z
        // preserves its prior value. The BIO value owns all resulting state.
        // Complexity review: six fixed-width numeric reads and at most two
        // bounded string conversions; no line rescans beyond fixed fields.
        let read_fixed_double = |offset: usize, width: usize| {
            let suffix = line.get(offset..).unwrap_or_default();
            read_double(&suffix[..suffix.len().min(width)])
        };

        if line.len() > 54 {
            set_crystal_cell(
                &mut self.crystal_info,
                BioCrystalCell {
                    a: read_fixed_double(6, 9),
                    b: read_fixed_double(15, 9),
                    c: read_fixed_double(24, 9),
                    alpha: read_fixed_double(33, 7),
                    beta: read_fixed_double(40, 7),
                    gamma: read_fixed_double(47, 7),
                },
            )?;
        }
        if line.len() > 56 {
            let value = std::str::from_utf8(&remark3_fixed_field(line, 55, 11))
                .map(str::to_owned)
                .map_err(|error| PdbCrystalError::TextFieldNotUtf8 {
                    field_offset: 55,
                    valid_up_to: error.valid_up_to(),
                })?;
            set_crystal_space_group_hm(&mut self.crystal_info, value);
        }
        if line.len() > 67 {
            let value = std::str::from_utf8(&remark3_fixed_field(line, 66, 4))
                .map(str::to_owned)
                .map_err(|error| PdbCrystalError::TextFieldNotUtf8 {
                    field_offset: 66,
                    valid_up_to: error.valid_up_to(),
                })?;
            set_crystal_z_pdb_if_nonempty(&mut self.crystal_info, value);
        }
        Ok(())
    }

    fn seqres_record(&mut self, line: &[u8]) -> bool {
        self.entities.seqres_record(line)
    }

    fn remark_record(
        &mut self,
        source_line_buffer: &[u8; 122],
        line_len: usize,
    ) -> Result<bool, PdbRemarkError> {
        // Gemmi source: third_party/gemmi/src/pdb.cpp,
        // populate_structure_from_pdb_stream.
        // Gemmi❗✔️:     } else if (is_record_type4(line, "REMARK")) {
        // Gemmi❗✔️:       if (line[len-1] == '\n')
        // Gemmi❗✔️:         --len;
        // Gemmi❗✔️:       if (line[len-1] == '\r')
        // Gemmi❗✔️:         --len;
        // Gemmi❗✔️:       st.raw_remarks.emplace_back(line, line+len);
        // The reader invokes this capture branch independently of the later
        // `skip_remarks`-guarded metadata interpretation. Remove at most one
        // terminal LF and then at most one terminal CR, retaining all other
        // bytes and row order. The canonical destination is the existing
        // `BioStructureSourceState::raw_remarks`; invalid UTF-8 is reported explicitly
        // because its String value cannot retain arbitrary source bytes.
        // Behavior review: only a matching four-byte record prefix is stored;
        // trailing byte removal is ordered and single-occurrence, not trim().
        // Complexity review: one bounded source-line scan for UTF-8 validity
        // and one O(line_len) owned copy/appending String, matching Gemmi's
        // string copy plus amortized vector append.
        if line_len > source_line_buffer.len() {
            return Err(PdbRemarkError::LineLengthOutsideBuffer { line_len });
        }
        if !gemmi_record_type4(source_line_buffer, 0, *b"REMA") {
            return Ok(false);
        }

        let mut content_end = line_len;
        if content_end > 0 && source_line_buffer[content_end - 1] == b'\n' {
            content_end -= 1;
        }
        if content_end > 0 && source_line_buffer[content_end - 1] == b'\r' {
            content_end -= 1;
        }
        let remark =
            String::from_utf8(source_line_buffer[..content_end].to_vec()).map_err(|error| {
                PdbRemarkError::TextFieldNotUtf8 {
                    valid_up_to: error.utf8_error().valid_up_to(),
                }
            })?;
        self.source_state.raw_remarks.push(remark);
        Ok(true)
    }

    fn remark3_record(&mut self, line: &str) -> Result<(), PdbRemark3Error> {
        // Gemmi source: third_party/gemmi/src/pdb.cpp, read_remark3_line.
        // Gemmi❗✔️:   const char* key_start = skip_blank(line + 10);
        // Gemmi❗✔️:   const char* colon = std::strchr(key_start, ':');
        // Gemmi❗✔️:   const char* key_end = rtrim_cstr(key_start, colon);
        // Gemmi❗✔️:   std::string key(key_start, key_end);
        // Gemmi❗✔️:   if (possibly_unfinished_remark3) {
        // Gemmi❗✔️:     if (key_start > line + 17) {
        // Gemmi❗✔️:       *possibly_unfinished_remark3 += ' ';
        // Gemmi❗✔️:       possibly_unfinished_remark3->append(key);
        // Gemmi❗✔️:       return;
        // Gemmi❗✔️:     }
        // Gemmi❗✔️:     possibly_unfinished_remark3 = nullptr;
        // Gemmi❗✔️:   }
        // Gemmi❗✔️:   if (colon) {
        // Gemmi❗✔️:     const char* value = skip_blank(colon + 1);
        // Gemmi❗✔️:     const char* end = rtrim_cstr(value);
        // Gemmi❗✔️:     if (end - value == 4 && std::strncmp(value, "NULL", 4) == 0)
        // Gemmi❗✔️:       return;
        // `line` is the already-captured REMARK 3 text with source LF/CR
        // removal applied. The source's `skip_blank` is space/tab only;
        // `rtrim_cstr` uses C-locale `isspace`. Keep these two byte predicates
        // distinct. Continuation targets are represented by stable vector
        // indices rather than a pointer into a `Vec` element; P27 installs the
        // target when the source creates a TLS selection.
        let bytes = line.as_bytes();
        let c_string_end = bytes
            .iter()
            .position(|byte| *byte == 0)
            .unwrap_or(bytes.len());
        let bytes = &bytes[..c_string_end];
        if bytes.len() <= 11 {
            return Ok(());
        }

        let mut key_start = 10;
        while key_start < bytes.len() && matches!(bytes[key_start], b' ' | b'\t') {
            key_start += 1;
        }
        let colon = bytes[key_start..]
            .iter()
            .position(|byte| *byte == b':')
            .map(|offset| key_start + offset);
        let key_end_bound = colon.unwrap_or(bytes.len());
        let key_end = (key_start..key_end_bound)
            .rfind(|index| !gemmi_is_space(bytes[*index]))
            .map_or(key_start, |index| index + 1);
        let key = &bytes[key_start..key_end];

        if let Some(target) = self.remark3_continuation {
            if key_start > 17 {
                let selection = self
                    .metadata
                    .refinement
                    .get_mut(target.refinement_index)
                    .and_then(|refinement| refinement.tls_groups.get_mut(target.tls_group_index))
                    .and_then(|group| group.selections.get_mut(target.selection_index))
                    .ok_or(PdbRemark3Error::MissingContinuationTarget)?;
                selection.details.push(' ');
                selection
                    .details
                    .push_str(std::str::from_utf8(key).map_err(|error| {
                        PdbRemark3Error::TextNotUtf8 {
                            valid_up_to: error.valid_up_to(),
                        }
                    })?);
                self.update_remark3_resolution();
                return Ok(());
            }
            self.remark3_continuation = None;
        }

        let Some(colon) = colon else {
            if key == b"DATA USED IN REFINEMENT." {
                let mut refinement = BioRefinementInfo::default();
                refinement.id = (self.metadata.refinement.len() + 1).to_string();
                self.metadata.refinement.push(refinement);
            } else if key == b"FIT IN THE HIGHEST RESOLUTION BIN." {
                if let Some(refinement) = self.metadata.refinement.last_mut() {
                    refinement.bins.push(Default::default());
                }
            }
            self.update_remark3_resolution();
            return Ok(());
        };

        let mut value_start = colon + 1;
        while value_start < bytes.len() && matches!(bytes[value_start], b' ' | b'\t') {
            value_start += 1;
        }
        let value_end = (value_start..bytes.len())
            .rfind(|index| !gemmi_is_space(bytes[*index]))
            .map_or(value_start, |index| index + 1);
        let value = &bytes[value_start..value_end];
        if value == b"NULL" {
            self.update_remark3_resolution();
            return Ok(());
        }

        // Gemmi source helper: third_party/gemmi/src/pdb.cpp, read_remark3_line.
        // Gemmi❗✔️:     if (same_str(key, "PROGRAM"))
        // Gemmi❗✔️:       add_software(meta, SoftwareItem::Refinement, std::string(value, end));
        if key == b"PROGRAM" {
            self.add_software(BioSoftwareClassification::Refinement, value)?;
        }
        let Some(refinement_index) = self.metadata.refinement.len().checked_sub(1) else {
            self.update_remark3_resolution();
            return Ok(());
        };
        let refinement = &mut self.metadata.refinement[refinement_index];

        // Gemmi❗✔️:     if (same_str(key, "RESOLUTION RANGE HIGH (ANGSTROMS)")) {
        // Gemmi❗✔️:       ref_info.resolution_high = fast_atof(value);
        if key == b"RESOLUTION RANGE HIGH (ANGSTROMS)" {
            refinement.basic.resolution_high = read_double(value);
        // Gemmi❗✔️:     } else if (same_str(key, "RESOLUTION RANGE LOW  (ANGSTROMS)")) {
        // Gemmi❗✔️:       ref_info.resolution_low = fast_atof(value);
        } else if key == b"RESOLUTION RANGE LOW  (ANGSTROMS)" {
            refinement.basic.resolution_low = read_double(value);
        // Gemmi❗✔️:     } else if (same_str(key, "COMPLETENESS FOR RANGE        (%)")) {
        // Gemmi❗✔️:       ref_info.completeness = fast_atof(value);
        } else if key == b"COMPLETENESS FOR RANGE        (%)" {
            refinement.basic.completeness = read_double(value);
        // Gemmi❗✔️:     } else if (same_str(key, "NUMBER OF REFLECTIONS")) {
        // Gemmi❗✔️:       ref_info.reflection_count = std::atoi(value);
        } else if key == b"NUMBER OF REFLECTIONS" {
            refinement.basic.reflection_count = remark3_int(value)?;
        // Gemmi❗✔️:     } else if (same_str(key, "CROSS-VALIDATION METHOD")) {
        // Gemmi❗✔️:       ref_info.cross_validation_method = std::string(value, end);
        } else if key == b"CROSS-VALIDATION METHOD" {
            refinement.cross_validation_method = remark3_text(value)?;
        // Gemmi❗✔️:     } else if (same_str(key, "FREE R VALUE TEST SET SELECTION")) {
        // Gemmi❗✔️:       ref_info.rfree_selection_method = std::string(value, end);
        } else if key == b"FREE R VALUE TEST SET SELECTION" {
            refinement.rfree_selection_method = remark3_text(value)?;
        // Gemmi❗✔️:     } else if (same_str(key, "R VALUE     (WORKING + TEST SET)")) {
        // Gemmi❗✔️:       ref_info.r_all = fast_atof(value);
        } else if key == b"R VALUE     (WORKING + TEST SET)" {
            refinement.basic.r_all = read_double(value);
        // Gemmi❗✔️:     } else if (same_str(key, "R VALUE            (WORKING SET)")) {
        // Gemmi❗✔️:       ref_info.r_work = fast_atof(value);
        } else if key == b"R VALUE            (WORKING SET)" {
            refinement.basic.r_work = read_double(value);
        // Gemmi❗✔️:     } else if (same_str(key, "FREE R VALUE")) {
        // Gemmi❗✔️:       ref_info.r_free = fast_atof(value);
        } else if key == b"FREE R VALUE" {
            refinement.basic.r_free = read_double(value);
        // Gemmi❗✔️:     } else if (same_str(key, "FREE R VALUE TEST SET COUNT")) {
        // Gemmi❗✔️:       ref_info.rfree_set_count = atoi(value);
        } else if key == b"FREE R VALUE TEST SET COUNT" {
            refinement.basic.rfree_set_count = remark3_int(value)?;
        // Gemmi❗✔️:     } else if (same_str(key, "TOTAL NUMBER OF BINS USED")) {
        // Gemmi❗✔️:       ref_info.bin_count = std::atoi(value);
        } else if key == b"TOTAL NUMBER OF BINS USED" {
            refinement.bin_count = remark3_int(value)?;
        // Gemmi❗✔️:     } else if (same_str(key, "BIN RESOLUTION RANGE HIGH       (A)")) {
        // Gemmi❗✔️:       if (!ref_info.bins.empty())
        // Gemmi❗✔️:         ref_info.bins.back().resolution_high = fast_atof(value);
        } else if key == b"BIN RESOLUTION RANGE HIGH       (A)" {
            if let Some(bin) = refinement.bins.last_mut() {
                bin.resolution_high = read_double(value);
            }
        // Gemmi❗✔️:     } else if (same_str(key, "BIN RESOLUTION RANGE LOW        (A)")) {
        // Gemmi❗✔️:       if (!ref_info.bins.empty())
        // Gemmi❗✔️:         ref_info.bins.back().resolution_low = fast_atof(value);
        } else if key == b"BIN RESOLUTION RANGE LOW        (A)" {
            if let Some(bin) = refinement.bins.last_mut() {
                bin.resolution_low = read_double(value);
            }
        // Gemmi❗✔️:     } else if (same_str(key, "BIN COMPLETENESS (WORKING+TEST) (%)")) {
        // Gemmi❗✔️:       if (!ref_info.bins.empty())
        // Gemmi❗✔️:         ref_info.bins.back().completeness = fast_atof(value);
        } else if key == b"BIN COMPLETENESS (WORKING+TEST) (%)" {
            if let Some(bin) = refinement.bins.last_mut() {
                bin.completeness = read_double(value);
            }
        // Gemmi❗✔️:     } else if (same_str(key, "REFLECTIONS IN BIN   (WORKING+TEST)")) {
        // Gemmi❗✔️:       if (!ref_info.bins.empty())
        // Gemmi❗✔️:         ref_info.bins.back().reflection_count = std::atoi(value);
        } else if key == b"REFLECTIONS IN BIN   (WORKING+TEST)" {
            if let Some(bin) = refinement.bins.last_mut() {
                bin.reflection_count = remark3_int(value)?;
            }
        // Gemmi❗✔️:     } else if (same_str(key, "BIN R VALUE          (WORKING+TEST)")) {
        // Gemmi❗✔️:       if (!ref_info.bins.empty())
        // Gemmi❗✔️:         ref_info.bins.back().r_all = fast_atof(value);
        } else if key == b"BIN R VALUE          (WORKING+TEST)" {
            if let Some(bin) = refinement.bins.last_mut() {
                bin.r_all = read_double(value);
            }
        // Gemmi❗✔️:     } else if (same_str(key, "REFLECTIONS IN BIN    (WORKING SET)")) {
        // Gemmi❗✔️:       if (!ref_info.bins.empty())
        // Gemmi❗✔️:         ref_info.bins.back().work_set_count = std::atoi(value);
        } else if key == b"REFLECTIONS IN BIN    (WORKING SET)" {
            if let Some(bin) = refinement.bins.last_mut() {
                bin.work_set_count = remark3_int(value)?;
            }
        // Gemmi❗✔️:     } else if (same_str(key, "BIN R VALUE           (WORKING SET)")) {
        // Gemmi❗✔️:       if (!ref_info.bins.empty())
        // Gemmi❗✔️:         ref_info.bins.back().r_work = fast_atof(value);
        } else if key == b"BIN R VALUE           (WORKING SET)" {
            if let Some(bin) = refinement.bins.last_mut() {
                bin.r_work = read_double(value);
            }
        // Gemmi❗✔️:     } else if (same_str(key, "BIN FREE R VALUE")) {
        // Gemmi❗✔️:       if (!ref_info.bins.empty())
        // Gemmi❗✔️:         ref_info.bins.back().r_free = fast_atof(value);
        } else if key == b"BIN FREE R VALUE" {
            if let Some(bin) = refinement.bins.last_mut() {
                bin.r_free = read_double(value);
            }
        // Gemmi❗✔️:     } else if (same_str(key, "BIN FREE R VALUE TEST SET COUNT")) {
        // Gemmi❗✔️:       if (!ref_info.bins.empty())
        // Gemmi❗✔️:         ref_info.bins.back().rfree_set_count = std::atoi(value);
        } else if key == b"BIN FREE R VALUE TEST SET COUNT" {
            if let Some(bin) = refinement.bins.last_mut() {
                bin.rfree_set_count = remark3_int(value)?;
            }
        // Gemmi❗✔️:     } else if (same_str(key, "FROM WILSON PLOT           (A**2)")) {
        // Gemmi❗✔️:       // TODO
        // Gemmi❗✔️:       // exper.b_wilson = fast_atof(value);
        } else if key == b"FROM WILSON PLOT           (A**2)" {
            // The pinned source branch is itself an explicit TODO and performs
            // no assignment; retain the existing canonical state unchanged.
            // Gemmi❗✔️:     } else if (same_str(key, "MEAN B VALUE      (OVERALL, A**2)")) {
            // Gemmi❗✔️:       ref_info.mean_b = fast_atof(value);
        } else if key == b"MEAN B VALUE      (OVERALL, A**2)" {
            refinement.mean_b = read_double(value);
        // Gemmi❗✔️:     } else if (same_str(key, "B11 (A**2)")) {
        // Gemmi❗✔️:       ref_info.aniso_b.u11 = fast_atof(value);
        } else if key == b"B11 (A**2)" {
            refinement.aniso_b[0] = read_double(value);
        // Gemmi❗✔️:     } else if (same_str(key, "B22 (A**2)")) {
        // Gemmi❗✔️:       ref_info.aniso_b.u22 = fast_atof(value);
        } else if key == b"B22 (A**2)" {
            refinement.aniso_b[1] = read_double(value);
        // Gemmi❗✔️:     } else if (same_str(key, "B33 (A**2)")) {
        // Gemmi❗✔️:       ref_info.aniso_b.u33 = fast_atof(value);
        } else if key == b"B33 (A**2)" {
            refinement.aniso_b[2] = read_double(value);
        // Gemmi❗✔️:     } else if (same_str(key, "B12 (A**2)")) {
        // Gemmi❗✔️:       ref_info.aniso_b.u12 = fast_atof(value);
        } else if key == b"B12 (A**2)" {
            refinement.aniso_b[3] = read_double(value);
        // Gemmi❗✔️:     } else if (same_str(key, "B13 (A**2)")) {
        // Gemmi❗✔️:       ref_info.aniso_b.u13 = fast_atof(value);
        } else if key == b"B13 (A**2)" {
            refinement.aniso_b[4] = read_double(value);
        // Gemmi❗✔️:     } else if (same_str(key, "B23 (A**2)")) {
        // Gemmi❗✔️:       ref_info.aniso_b.u23 = fast_atof(value);
        } else if key == b"B23 (A**2)" {
            refinement.aniso_b[5] = read_double(value);
        // Gemmi❗✔️:     } else if (same_str(key, "ESD FROM LUZZATI PLOT                    (A)")) {
        // Gemmi❗✔️:       ref_info.luzzati_error = fast_atof(value);
        } else if key == b"ESD FROM LUZZATI PLOT                    (A)" {
            refinement.luzzati_error = read_double(value);
        // Gemmi❗✔️:     } else if (same_str(key, "DPI (BLOW EQ-10) BASED ON R VALUE        (A)")) {
        // Gemmi❗✔️:       ref_info.dpi_blow_r = fast_atof(value);
        } else if key == b"DPI (BLOW EQ-10) BASED ON R VALUE        (A)" {
            refinement.dpi_blow_r = read_double(value);
        // Gemmi❗✔️:     } else if (same_str(key, "DPI (BLOW EQ-9) BASED ON FREE R VALUE    (A)")) {
        // Gemmi❗✔️:       ref_info.dpi_blow_rfree = fast_atof(value);
        } else if key == b"DPI (BLOW EQ-9) BASED ON FREE R VALUE    (A)" {
            refinement.dpi_blow_rfree = read_double(value);
        // Gemmi❗✔️:     } else if (same_str(key, "DPI (CRUICKSHANK) BASED ON R VALUE       (A)")) {
        // Gemmi❗✔️:       ref_info.dpi_cruickshank_r = fast_atof(value);
        } else if key == b"DPI (CRUICKSHANK) BASED ON R VALUE       (A)" {
            refinement.dpi_cruickshank_r = read_double(value);
        // Gemmi❗✔️:     } else if (same_str(key, "DPI (CRUICKSHANK) BASED ON FREE R VALUE  (A)")) {
        // Gemmi❗✔️:       ref_info.dpi_cruickshank_rfree = fast_atof(value);
        } else if key == b"DPI (CRUICKSHANK) BASED ON FREE R VALUE  (A)" {
            refinement.dpi_cruickshank_rfree = read_double(value);
        // Gemmi❗✔️:     } else if (same_str(key, "CORRELATION COEFFICIENT FO-FC")) {
        // Gemmi❗✔️:       ref_info.cc_fo_fc_work = fast_atof(value);
        } else if key == b"CORRELATION COEFFICIENT FO-FC" {
            refinement.basic.cc_fo_fc_work = read_double(value);
        // Gemmi❗✔️:     } else if (same_str(key, "CORRELATION COEFFICIENT FO-FC FREE")) {
        // Gemmi❗✔️:       ref_info.cc_fo_fc_free = fast_atof(value);
        } else if key == b"CORRELATION COEFFICIENT FO-FC FREE" {
            refinement.basic.cc_fo_fc_free = read_double(value);
        // Gemmi source: third_party/gemmi/src/pdb.cpp, read_remark3_line.
        // Gemmi❗✔️:     } else if (same_str(key, "BOND LENGTHS")) {
        // Gemmi❗✔️:       add_restraint_count_weight(ref_info, "t_bond_d", value);
        // Gemmi❗✔️:     } else if (same_str(key, "BOND ANGLES")) {
        // Gemmi❗✔️:       add_restraint_count_weight(ref_info, "t_angle_deg", value);
        // Gemmi❗✔️:     } else if (same_str(key, "TORSION ANGLES")) {
        // Gemmi❗✔️:       add_restraint_count_weight(ref_info, "t_dihedral_angle_d", value);
        // Gemmi❗✔️:     } else if (same_str(key, "TRIGONAL CARBON PLANES")) {
        // Gemmi❗✔️:       add_restraint_count_weight(ref_info, "t_trig_c_planes", value);
        // Gemmi❗✔️:     } else if (same_str(key, "GENERAL PLANES")) {
        // Gemmi❗✔️:       add_restraint_count_weight(ref_info, "t_gen_planes", value);
        // Gemmi❗✔️:     } else if (same_str(key, "ISOTROPIC THERMAL FACTORS")) {
        // Gemmi❗✔️:       add_restraint_count_weight(ref_info, "t_it", value);
        // Gemmi❗✔️:     } else if (same_str(key, "BAD NON-BONDED CONTACTS")) {
        // Gemmi❗✔️:       add_restraint_count_weight(ref_info, "t_nbd", value);
        // Gemmi❗✔️:     } else if (same_str(key, "IMPROPER TORSIONS")) {
        // Gemmi❗✔️:       add_restraint_count_weight(ref_info, "t_improper_torsion", value);
        // Gemmi❗✔️:     } else if (same_str(key, "CHIRAL IMPROPER TORSION")) {
        // Gemmi❗✔️:       add_restraint_count_weight(ref_info, "t_chiral_improper_torsion", value);
        // Gemmi❗✔️:     } else if (same_str(key, "SUM OF OCCUPANCIES")) {
        // Gemmi❗✔️:       add_restraint_count_weight(ref_info, "t_sum_occupancies", value);
        // Gemmi❗✔️:     } else if (same_str(key, "UTILITY DISTANCES")) {
        // Gemmi❗✔️:       add_restraint_count_weight(ref_info, "t_utility_distance", value);
        // Gemmi❗✔️:     } else if (same_str(key, "UTILITY ANGLES")) {
        // Gemmi❗✔️:       add_restraint_count_weight(ref_info, "t_utility_angle", value);
        // Gemmi❗✔️:     } else if (same_str(key, "UTILITY TORSION")) {
        // Gemmi❗✔️:       add_restraint_count_weight(ref_info, "t_utility_torsion", value);
        // Gemmi❗✔️:     } else if (same_str(key, "IDEAL-DIST CONTACT TERM")) {
        // Gemmi❗✔️:       add_restraint_count_weight(ref_info, "t_ideal_dist_contact", value);
        } else if key == b"BOND LENGTHS" {
            gemmi_add_restraint_count_weight(refinement, "t_bond_d", value)?;
        } else if key == b"BOND ANGLES" {
            gemmi_add_restraint_count_weight(refinement, "t_angle_deg", value)?;
        } else if key == b"TORSION ANGLES" {
            gemmi_add_restraint_count_weight(refinement, "t_dihedral_angle_d", value)?;
        } else if key == b"TRIGONAL CARBON PLANES" {
            gemmi_add_restraint_count_weight(refinement, "t_trig_c_planes", value)?;
        } else if key == b"GENERAL PLANES" {
            gemmi_add_restraint_count_weight(refinement, "t_gen_planes", value)?;
        } else if key == b"ISOTROPIC THERMAL FACTORS" {
            gemmi_add_restraint_count_weight(refinement, "t_it", value)?;
        } else if key == b"BAD NON-BONDED CONTACTS" {
            gemmi_add_restraint_count_weight(refinement, "t_nbd", value)?;
        } else if key == b"IMPROPER TORSIONS" {
            gemmi_add_restraint_count_weight(refinement, "t_improper_torsion", value)?;
        } else if key == b"CHIRAL IMPROPER TORSION" {
            gemmi_add_restraint_count_weight(refinement, "t_chiral_improper_torsion", value)?;
        } else if key == b"SUM OF OCCUPANCIES" {
            gemmi_add_restraint_count_weight(refinement, "t_sum_occupancies", value)?;
        } else if key == b"UTILITY DISTANCES" {
            gemmi_add_restraint_count_weight(refinement, "t_utility_distance", value)?;
        } else if key == b"UTILITY ANGLES" {
            gemmi_add_restraint_count_weight(refinement, "t_utility_angle", value)?;
        } else if key == b"UTILITY TORSION" {
            gemmi_add_restraint_count_weight(refinement, "t_utility_torsion", value)?;
        } else if key == b"IDEAL-DIST CONTACT TERM" {
            gemmi_add_restraint_count_weight(refinement, "t_ideal_dist_contact", value)?;
        // Gemmi source: third_party/gemmi/src/pdb.cpp, read_remark3_line.
        // Gemmi❗✔️:     } else if (same_str(key, "TLS GROUP")) {
        // Gemmi❗✔️:       ref_info.tls_groups.emplace_back();
        // Gemmi❗✔️:       TlsGroup& tls_group = ref_info.tls_groups.back();
        // Gemmi❗✔️:       tls_group.id = std::string(value, end);
        // Gemmi❗✔️:       tls_group.num_id = (short) no_sign_atoi(tls_group.id.c_str());
        // Gemmi❗✔️:     } else if (same_str(key, "SET") ||
        // Gemmi❗✔️:                // "REMARK   3    SELECTION:"            -> TLS
        // Gemmi❗✔️:                // "REMARK   3     SELECTION          :" -> NCS
        // Gemmi❗✔️:                (same_str(key, "SELECTION") && colon == line + 23)) {
        // Gemmi❗✔️:       if (!ref_info.tls_groups.empty()) {
        // Gemmi❗✔️:         TlsGroup& group = ref_info.tls_groups.back();
        // Gemmi❗✔️:         group.selections.emplace_back();
        // Gemmi❗✔️:         group.selections.back().details = std::string(value, end);
        // Gemmi❗✔️:         possibly_unfinished_remark3 = &group.selections.back().details;
        // Gemmi❗✔️:       }
        // Gemmi❗✔️:     } else if (same_str(key, "RESIDUE RANGE")) {
        // Gemmi❗✔️:       if (!ref_info.tls_groups.empty() && end > colon+21) {
        // Gemmi❗✔️:         TlsGroup& group = ref_info.tls_groups.back();
        // Gemmi❗✔️:         group.selections.emplace_back();
        // Gemmi❗✔️:         TlsGroup::Selection& sel = group.selections.back();
        // Gemmi❗✔️:         sel.chain = read_string(colon+1, 5);
        // Gemmi❗✔️:         if (sel.chain == read_string(colon+16, 5)) {
        // Gemmi❗✔️:           try {
        // Gemmi❗✔️:             sel.res_begin = SeqId(read_string(colon+6, 6));
        // Gemmi❗✔️:             sel.res_end = SeqId(read_string(colon+21, 6));
        // Gemmi❗✔️:           } catch (std::invalid_argument&) {
        // Gemmi❗✔️:             group.selections.pop_back();
        // Gemmi❗✔️:           }
        // Gemmi❗✔️:         } else {  // unexpected -- TLS group should be in one chain
        // Gemmi❗✔️:           group.selections.pop_back();
        // Gemmi❗✔️:         }
        // Gemmi❗✔️:       }
        // Gemmi❗✔️:     } else if (same_str(key, "ORIGIN FOR THE GROUP (A)")) {
        // Gemmi❗✔️:       std::vector<std::string> xyz = split_str_multi(std::string(value, end));
        // Gemmi❗✔️:       if (ref_info.tls_groups.empty() || xyz.size() != 3)
        // Gemmi❗✔️:         return;
        // Gemmi❗✔️:       Position& origin = ref_info.tls_groups.back().origin;
        // Gemmi❗✔️:       origin.x = fast_atof(xyz[0].c_str());
        // Gemmi❗✔️:       origin.y = fast_atof(xyz[1].c_str());
        // Gemmi❗✔️:       origin.z = fast_atof(xyz[2].c_str());
        // Gemmi❗✔️:     } else if (is_tls_item(key)) {
        // Gemmi❗✔️:       if (ref_info.tls_groups.empty())
        // Gemmi❗✔️:         return;
        // Gemmi❗✔️:       TlsGroup& tls = ref_info.tls_groups.back();
        // Gemmi❗✔️:       std::vector<std::string> tokens = split_str_multi(key_start);
        // Gemmi❗✔️:       for (size_t i = 0; i + 1 < tokens.size(); i += 2) {
        // Gemmi❗✔️:         std::string& k = tokens[i];
        // Gemmi❗✔️:         if (k.size() == 4 && k[3] == ':')
        // Gemmi❗✔️:           k.resize(3);
        // Gemmi❗✔️:         if (is_tls_item(k)) {
        // Gemmi❗✔️:           int x = k[1] - '1';
        // Gemmi❗✔️:           int y = k[2] - '1';
        // Gemmi❗✔️:           double v = fast_atof(tokens[i+1].c_str());
        // Gemmi❗✔️:           if (k[0] == 'S') {
        // Gemmi❗✔️:             tls.S[x][y] = v;
        // Gemmi❗✔️:           } else {
        // Gemmi❗✔️:             SMat33<double>& tensor = k[0] == 'T' ? tls.T : tls.L;
        // Gemmi❗✔️:             tensor.unchecked_ref(x, y) = v;
        // Gemmi❗✔️:           }
        // Gemmi❗✔️:         }
        // Gemmi❗✔️:       }
        // Gemmi❗✔️:     }
        // Behavior review: the source's exact key/column gates, last-refinement
        // and last-group ownership, fixed range offsets, continuation target,
        // no-sign IDs, token-pair order, and T/L/S assignments are represented.
        // Source-undefined integer overflow and the approved chain-ID width
        // boundary return structured errors. A short final fixed-width suffix
        // follows the source line buffer's zero termination instead.
        // Complexity review: the parser scans each record once for its key and
        // once for selected tokenization; token slices borrow source bytes and
        // do not clone topology or whole metadata state.
        } else if key == b"TLS GROUP" {
            let id = remark3_text(value)?;
            let (numeric_id, _) = gemmi_no_sign_atoi(id.as_bytes())
                .ok_or(PdbRemark3Error::IntegerOutsideSourceDefinedRange)?;
            refinement.tls_groups.push(BioTlsGroup {
                id,
                num_id: numeric_id as i16,
                ..BioTlsGroup::default()
            });
        } else if key == b"SET" || (key == b"SELECTION" && colon == 23) {
            if let Some(group_index) = refinement.tls_groups.len().checked_sub(1) {
                let details = remark3_text(value)?;
                let group = &mut refinement.tls_groups[group_index];
                group.selections.push(BioTlsSelection {
                    details,
                    ..BioTlsSelection::default()
                });
                self.remark3_continuation = Some(PdbRemark3Continuation {
                    refinement_index,
                    tls_group_index: group_index,
                    selection_index: group.selections.len() - 1,
                });
            }
        } else if key == b"RESIDUE RANGE" {
            if !refinement.tls_groups.is_empty() && value_end > colon + 21 {
                let chain_begin = remark3_fixed_field(bytes, colon + 1, 5);
                let sequence_begin = remark3_fixed_field(bytes, colon + 6, 6);
                let chain_end = remark3_fixed_field(bytes, colon + 16, 5);
                let sequence_end = remark3_fixed_field(bytes, colon + 21, 6);
                if chain_begin == chain_end {
                    if let (Some(res_begin), Some(res_end)) = (
                        gemmi_tls_seq_id(sequence_begin)?,
                        gemmi_tls_seq_id(sequence_end)?,
                    ) {
                        let chain = PdbChainId::from_ascii(chain_begin).ok_or(
                            PdbRemark3Error::TlsChainIdNotRepresentable {
                                width: chain_begin.len(),
                            },
                        )?;
                        refinement
                            .tls_groups
                            .last_mut()
                            .expect("checked above")
                            .selections
                            .push(BioTlsSelection {
                                chain,
                                res_begin,
                                res_end,
                                ..BioTlsSelection::default()
                            });
                    }
                }
            }
        } else if key == b"ORIGIN FOR THE GROUP (A)" {
            if let Some(group) = refinement.tls_groups.last_mut() {
                let xyz = gemmi_split_str_multi(value);
                if let [x, y, z] = xyz.as_slice() {
                    group.origin = [read_double(x), read_double(y), read_double(z)];
                }
            }
        } else if gemmi_is_tls_item(key) {
            if let Some(group) = refinement.tls_groups.last_mut() {
                for pair in gemmi_split_str_multi(&bytes[key_start..]).chunks_exact(2) {
                    let mut item_key = pair[0];
                    if item_key.len() == 4 && item_key[3] == b':' {
                        item_key = &item_key[..3];
                    }
                    if !gemmi_is_tls_item(item_key) {
                        continue;
                    }
                    let matrix = item_key[0];
                    let row = usize::from(item_key[1] - b'1');
                    let column = usize::from(item_key[2] - b'1');
                    let tensor_value = read_double(pair[1]);
                    match matrix {
                        b'S' => group.s[row][column] = tensor_value,
                        b'T' | b'L' => {
                            let component = gemmi_tls_symmetric_component_index(row, column);
                            if matrix == b'T' {
                                group.t[component] = tensor_value;
                            } else {
                                group.l[component] = tensor_value;
                            }
                        }
                        _ => unreachable!("TLS key was validated above"),
                    }
                }
            }
        }

        self.update_remark3_resolution();
        Ok(())
    }

    fn remark_200_230_240_record(&mut self, line: &str) -> Result<(), PdbRemark200Error> {
        // Gemmi source: third_party/gemmi/src/pdb.cpp,
        // read_remark_200_230_240.
        // Gemmi❗✔️: void read_remark_200_230_240(const char* line, Metadata& meta, std::string*& cryst_desc) {
        // Gemmi❗✔️:   // multi-line continuation requires special handling
        // Gemmi❗✔️:   if (cryst_desc) {
        // Gemmi❗✔️:     if (line[10] == ' ' && line[11] == ' ') {
        // Gemmi❗✔️:       const char* start = line + 11;
        // Gemmi❗✔️:       cryst_desc->append(start, rtrim_cstr(start) - start);
        // Gemmi❗✔️:       return;
        // Gemmi❗✔️:     }
        // Gemmi❗✔️:     cryst_desc = nullptr;
        // Gemmi❗✔️:   }
        // Gemmi❗✔️:
        // Gemmi❗✔️:   const char* key_start = skip_blank(line + 10);
        // Gemmi❗✔️:   const char* colon = std::strchr(key_start, ':');
        // Gemmi❗✔️:   const char* key_end = rtrim_cstr(key_start, colon);
        // Gemmi❗✔️:   std::string key(key_start, key_end);
        // Gemmi❗✔️:   if (colon) {
        // Gemmi❗✔️:     const char* value = skip_blank(colon + 1);
        // Gemmi❗✔️:     const char* end = rtrim_cstr(value);
        // Gemmi❗✔️:     if (end - value == 4 && std::strncmp(value, "NULL", 4) == 0)
        // Gemmi❗✔️:       return;
        // Gemmi❗✔️:     if (same_str(key, "INTENSITY-INTEGRATION SOFTWARE")) {
        // Gemmi❗✔️:       add_software(meta, SoftwareItem::DataReduction, std::string(value, end));
        // Gemmi❗✔️:     } else if (same_str(key, "DATA SCALING SOFTWARE")) {
        // Gemmi❗✔️:       add_software(meta, SoftwareItem::DataScaling, std::string(value, end));
        // Gemmi❗✔️:     } else if (same_str(key, "SOFTWARE USED")) {
        // Gemmi❗✔️:       add_software(meta, SoftwareItem::Phasing, std::string(value, end));
        // Gemmi❗✔️:     } else if (same_str(key, "METHOD USED TO DETERMINE THE STRUCTURE")) {
        // Gemmi❗✔️:       meta.solved_by = std::string(value, end);
        // Gemmi❗✔️:     } else if (same_str(key, "STARTING MODEL")) {
        // Gemmi❗✔️:       meta.starting_model = std::string(value, end);
        // Gemmi❗✔️:     } else if (!meta.experiments.empty()) {
        // Gemmi❗✔️:       ExperimentInfo& exper = meta.experiments.back();
        // Gemmi❗✔️:       DiffractionInfo& diffr = meta.crystals.back().diffractions[0];
        // Gemmi❗✔️:       if (same_str(key, "EXPERIMENT TYPE")) {
        // Gemmi❗✔️:         exper.method = std::string(value, end);
        // Gemmi❗✔️:       } else if (same_str(key, "NUMBER OF CRYSTALS USED")) {
        // Gemmi❗✔️:         exper.number_of_crystals = std::atoi(value);
        // Gemmi❗✔️:       } else if (same_str(key, "PH")) {
        // Gemmi❗✔️:         if (is_double(value))
        // Gemmi❗✔️:           meta.crystals.back().ph = fast_atof(value);
        // Gemmi❗✔️:         else
        // Gemmi❗✔️:           meta.crystals.back().ph_range = std::string(value, end);
        // Gemmi❗✔️:       } else if (same_str(key, "DATE OF DATA COLLECTION")) {
        // Gemmi❗✔️:         diffr.collection_date = pdb_date_format_to_iso(std::string(value, end));
        // Gemmi❗✔️:       } else if (same_str(key, "TEMPERATURE           (KELVIN)")) {
        // Gemmi❗✔️:         diffr.temperature = fast_atof(value);
        // Gemmi❗✔️:       } else if (same_str(key, "SYNCHROTRON              (Y/N)")) {
        // Gemmi❗✔️:         if (*value == 'Y')
        // Gemmi❗✔️:           diffr.source = "SYNCHROTRON";
        // Gemmi❗✔️:       } else if (same_str(key, "RADIATION SOURCE")) {
        // Gemmi❗✔️:         if (same_str(diffr.source, "SYNCHROTRON"))
        // Gemmi❗✔️:           diffr.synchrotron = std::string(value, end);
        // Gemmi❗✔️:         else
        // Gemmi❗✔️:           diffr.source = std::string(value, end);
        // Gemmi❗✔️:       } else if (same_str(key, "NEUTRON SOURCE")) {
        // Gemmi❗✔️:         diffr.source = std::string(value, end);
        // Gemmi❗✔️:       } else if (same_str(key, "BEAMLINE")) {
        // Gemmi❗✔️:         diffr.beamline = std::string(value, end);
        // Gemmi❗✔️:         if (!diffr.synchrotron.empty() && diffr.source_type.empty())
        // Gemmi❗✔️:           diffr.source_type = diffr.synchrotron + " BEAMLINE " + diffr.beamline;
        // Gemmi❗✔️:       } else if (same_str(key, "X-RAY GENERATOR MODEL")) {
        // Gemmi❗✔️:         diffr.source_type = std::string(value, end);
        // Gemmi❗✔️:       } else if (same_str(key, "MONOCHROMATIC OR LAUE    (M/L)")) {
        // Gemmi❗✔️:         diffr.mono_or_laue = *value;
        // Gemmi❗✔️:       } else if (same_str(key, "WAVELENGTH OR RANGE        (A)")) {
        // Gemmi❗✔️:         diffr.wavelengths = std::string(value, end);
        // Gemmi❗✔️:       } else if (same_str(key, "MONOCHROMATOR")) {
        // Gemmi❗✔️:         diffr.monochromator = std::string(value, end);
        // Gemmi❗✔️:       } else if (same_str(key, "OPTICS")) {
        // Gemmi❗✔️:         diffr.optics = std::string(value, end);
        // Gemmi❗✔️:       } else if (same_str(key, "DETECTOR TYPE")) {
        // Gemmi❗✔️:         diffr.detector = std::string(value, end);
        // Gemmi❗✔️:       } else if (same_str(key, "DETECTOR MANUFACTURER")) {
        // Gemmi❗✔️:         diffr.detector_make = std::string(value, end);
        // Gemmi❗✔️:       } else if (same_str(key, "NUMBER OF UNIQUE REFLECTIONS")) {
        // Gemmi❗✔️:         exper.unique_reflections = std::atoi(value);
        // Gemmi❗✔️:       } else if (same_str(key, "RESOLUTION RANGE HIGH      (A)")) {
        // Gemmi❗✔️:         exper.reflections.resolution_high = fast_atof(value);
        // Gemmi❗✔️:       } else if (same_str(key, "RESOLUTION RANGE LOW       (A)")) {
        // Gemmi❗✔️:         exper.reflections.resolution_low = fast_atof(value);
        // Gemmi❗✔️:       } else if (same_str(key, "COMPLETENESS FOR RANGE     (%)")) {
        // Gemmi❗✔️:         exper.reflections.completeness = fast_atof(value);
        // Gemmi❗✔️:       } else if (same_str(key, "DATA REDUNDANCY")) {
        // Gemmi❗✔️:         exper.reflections.redundancy = fast_atof(value);
        // Gemmi❗✔️:       } else if (same_str(key, "R MERGE                    (I)")) {
        // Gemmi❗✔️:         exper.reflections.r_merge = fast_atof(value);
        // Gemmi❗✔️:       } else if (same_str(key, "R SYM                      (I)")) {
        // Gemmi❗✔️:         exper.reflections.r_sym = fast_atof(value);
        // Gemmi❗✔️:       } else if (same_str(key, "<I/SIGMA(I)> FOR THE DATA SET")) {
        // Gemmi❗✔️:         exper.reflections.mean_I_over_sigma = fast_atof(value);
        // Gemmi❗✔️:       } else if (same_str(key, "REMARK")) {
        // Gemmi❗✔️:         cryst_desc = &meta.crystals.back().description;
        // Gemmi❗✔️:         *cryst_desc = std::string(value, end);
        // Gemmi❗✔️:       } else if (!exper.shells.empty()) {
        // Gemmi❗✔️:         if (same_str(key, "HIGHEST RESOLUTION SHELL, RANGE HIGH (A)")) {
        // Gemmi❗✔️:           exper.shells.back().resolution_high = fast_atof(value);
        // Gemmi❗✔️:         } else if (same_str(key, "HIGHEST RESOLUTION SHELL, RANGE LOW  (A)")) {
        // Gemmi❗✔️:           exper.shells.back().resolution_low = fast_atof(value);
        // Gemmi❗✔️:         } else if (same_str(key, "COMPLETENESS FOR SHELL     (%)")) {
        // Gemmi❗✔️:           exper.shells.back().completeness = fast_atof(value);
        // Gemmi❗✔️:         } else if (same_str(key, "DATA REDUNDANCY IN SHELL")) {
        // Gemmi❗✔️:           exper.shells.back().redundancy = fast_atof(value);
        // Gemmi❗✔️:         } else if (same_str(key, "R MERGE FOR SHELL          (I)")) {
        // Gemmi❗✔️:           exper.shells.back().r_merge = fast_atof(value);
        // Gemmi❗✔️:         } else if (same_str(key, "R SYM FOR SHELL            (I)")) {
        // Gemmi❗✔️:           exper.shells.back().r_sym = fast_atof(value);
        // Gemmi❗✔️:         } else if (same_str(key, "<I/SIGMA(I)> FOR SHELL")) {
        // Gemmi❗✔️:           exper.shells.back().mean_I_over_sigma = fast_atof(value);
        // Gemmi❗✔️:         }
        // Gemmi❗✔️:       }
        // Gemmi❗✔️:     }
        // Gemmi❗✔️:   } else {
        // Gemmi❗✔️:     if (same_str(key, "EXPERIMENTAL DETAILS")) {
        // Gemmi❗✔️:       meta.crystals.emplace_back();
        // Gemmi❗✔️:       CrystalInfo& c = meta.crystals.back();
        // Gemmi❗✔️:       c.id = std::to_string(meta.crystals.size());
        // Gemmi❗✔️:       c.diffractions.emplace_back();
        // Gemmi❗✔️:       c.diffractions[0].id = c.id;
        // Gemmi❗✔️:       meta.experiments.emplace_back();
        // Gemmi❗✔️:       meta.experiments.back().diffraction_ids.push_back(c.id);
        // Gemmi❗✔️:       if (line[8] == '0' && line[9] == '0')
        // Gemmi❗✔️:         c.diffractions[0].scattering_type = "x-ray";
        // Gemmi❗✔️:       else if (line[8] == '3' && line[9] == '0')
        // Gemmi❗✔️:         c.diffractions[0].scattering_type = "neutron";
        // Gemmi❗✔️:       else if (line[8] == '4' && line[9] == '0')
        // Gemmi❗✔️:         c.diffractions[0].scattering_type = "electron";
        // Gemmi❗✔️:     }
        // Gemmi❗✔️:     if (same_str(key, "IN THE HIGHEST RESOLUTION SHELL.")) {
        // Gemmi❗✔️:       if (!meta.experiments.empty())
        // Gemmi❗✔️:         meta.experiments.back().shells.emplace_back();
        // Gemmi❗✔️:     }
        // Gemmi❗✔️:   }
        // Gemmi❗✔️: }
        //
        // Behavior review: preserve continuation-before-dispatch, exact key/value
        // trimming, the source branch order, source-created experiment/crystal/
        // diffraction row alignment, and ordered metadata mutation. The pinned
        // `is_double` loop has a one-past-terminator path for some decimal inputs;
        // `gemmi_pdb_is_double` reports that undefined source case structurally
        // instead of reading beyond the string. C++ `atoi` overflow is likewise
        // outside source-defined behavior and maps to an explicit error.
        // Complexity review: fixed-column remarks are bounded by the PDB reader's
        // 120-byte line limit; parsing uses several O(n) scans and the same
        // ordered vector appends/string allocations as Gemmi. No graph traversal
        // or per-row nested lookup is added.
        let source_bytes = line.as_bytes();
        let c_end = source_bytes
            .iter()
            .position(|byte| *byte == 0)
            .unwrap_or(source_bytes.len());
        let line = &source_bytes[..c_end];

        if let Some(crystal_index) = self.remark200_continuation {
            if line.get(10).copied().unwrap_or(0) == b' '
                && line.get(11).copied().unwrap_or(0) == b' '
            {
                let end = gemmi_pdb_rtrim_cstr_end(line, 11, None);
                let continuation = std::str::from_utf8(&line[11..end]).map_err(|error| {
                    PdbRemark200Error::TextNotUtf8 {
                        valid_up_to: error.valid_up_to(),
                    }
                })?;
                let crystal = self
                    .metadata
                    .crystals
                    .get_mut(crystal_index)
                    .ok_or(PdbRemark200Error::MissingExperimentalRows)?;
                crystal.description.push_str(continuation);
                return Ok(());
            }
            self.remark200_continuation = None;
        }

        let key_start = gemmi_pdb_skip_blank(line, 10);
        let colon = line[key_start..]
            .iter()
            .position(|byte| *byte == b':')
            .map(|offset| key_start + offset);
        let key_end = gemmi_pdb_rtrim_cstr_end(line, key_start, colon);
        let key = &line[key_start..key_end];

        if let Some(colon) = colon {
            let value_start = gemmi_pdb_skip_blank(line, colon + 1);
            let value_end = gemmi_pdb_rtrim_cstr_end(line, value_start, None);
            let value = &line[value_start..value_end];
            let raw_value = &line[value_start..];
            if value == b"NULL" {
                return Ok(());
            }

            if key == b"INTENSITY-INTEGRATION SOFTWARE" {
                self.add_software(BioSoftwareClassification::DataReduction, value)
                    .map_err(PdbRemark200Error::Software)?;
            } else if key == b"DATA SCALING SOFTWARE" {
                self.add_software(BioSoftwareClassification::DataScaling, value)
                    .map_err(PdbRemark200Error::Software)?;
            } else if key == b"SOFTWARE USED" {
                self.add_software(BioSoftwareClassification::Phasing, value)
                    .map_err(PdbRemark200Error::Software)?;
            } else if key == b"METHOD USED TO DETERMINE THE STRUCTURE" {
                self.metadata.solved_by = std::str::from_utf8(value)
                    .map_err(|error| PdbRemark200Error::TextNotUtf8 {
                        valid_up_to: error.valid_up_to(),
                    })?
                    .to_owned();
            } else if key == b"STARTING MODEL" {
                self.metadata.starting_model = std::str::from_utf8(value)
                    .map_err(|error| PdbRemark200Error::TextNotUtf8 {
                        valid_up_to: error.valid_up_to(),
                    })?
                    .to_owned();
            } else if !self.metadata.experiments.is_empty() {
                let experiment_index = self.metadata.experiments.len() - 1;
                let crystal_index = self
                    .metadata
                    .crystals
                    .len()
                    .checked_sub(1)
                    .ok_or(PdbRemark200Error::MissingExperimentalRows)?;
                let metadata = &mut self.metadata;
                let experiment = metadata
                    .experiments
                    .get_mut(experiment_index)
                    .ok_or(PdbRemark200Error::MissingExperimentalRows)?;
                let crystal = metadata
                    .crystals
                    .get_mut(crystal_index)
                    .ok_or(PdbRemark200Error::MissingExperimentalRows)?;
                let diffraction = crystal
                    .diffractions
                    .get_mut(0)
                    .ok_or(PdbRemark200Error::MissingExperimentalRows)?;

                if key == b"EXPERIMENT TYPE" {
                    experiment.method = std::str::from_utf8(value)
                        .map_err(|error| PdbRemark200Error::TextNotUtf8 {
                            valid_up_to: error.valid_up_to(),
                        })?
                        .to_owned();
                } else if key == b"NUMBER OF CRYSTALS USED" {
                    experiment.number_of_crystals = remark3_int(raw_value)
                        .map_err(|_| PdbRemark200Error::IntegerOutsideSourceDefinedRange)?;
                } else if key == b"PH" {
                    if gemmi_pdb_is_double(raw_value)? {
                        crystal.ph = read_double(raw_value);
                    } else {
                        crystal.ph_range = std::str::from_utf8(value)
                            .map_err(|error| PdbRemark200Error::TextNotUtf8 {
                                valid_up_to: error.valid_up_to(),
                            })?
                            .to_owned();
                    }
                } else if key == b"DATE OF DATA COLLECTION" {
                    diffraction.collection_date = pdb_date_format_to_iso(value);
                } else if key == b"TEMPERATURE           (KELVIN)" {
                    diffraction.temperature = read_double(raw_value);
                } else if key == b"SYNCHROTRON              (Y/N)" {
                    if raw_value.first().copied().unwrap_or(0) == b'Y' {
                        diffraction.source = "SYNCHROTRON".to_owned();
                    }
                } else if key == b"RADIATION SOURCE" {
                    let source = std::str::from_utf8(value)
                        .map_err(|error| PdbRemark200Error::TextNotUtf8 {
                            valid_up_to: error.valid_up_to(),
                        })?
                        .to_owned();
                    if diffraction.source == "SYNCHROTRON" {
                        diffraction.synchrotron = source;
                    } else {
                        diffraction.source = source;
                    }
                } else if key == b"NEUTRON SOURCE" {
                    diffraction.source = std::str::from_utf8(value)
                        .map_err(|error| PdbRemark200Error::TextNotUtf8 {
                            valid_up_to: error.valid_up_to(),
                        })?
                        .to_owned();
                } else if key == b"BEAMLINE" {
                    diffraction.beamline = std::str::from_utf8(value)
                        .map_err(|error| PdbRemark200Error::TextNotUtf8 {
                            valid_up_to: error.valid_up_to(),
                        })?
                        .to_owned();
                    if !diffraction.synchrotron.is_empty() && diffraction.source_type.is_empty() {
                        diffraction.source_type = format!(
                            "{} BEAMLINE {}",
                            diffraction.synchrotron, diffraction.beamline
                        );
                    }
                } else if key == b"X-RAY GENERATOR MODEL" {
                    diffraction.source_type = std::str::from_utf8(value)
                        .map_err(|error| PdbRemark200Error::TextNotUtf8 {
                            valid_up_to: error.valid_up_to(),
                        })?
                        .to_owned();
                } else if key == b"MONOCHROMATIC OR LAUE    (M/L)" {
                    diffraction.mono_or_laue = raw_value.first().copied().unwrap_or(0);
                } else if key == b"WAVELENGTH OR RANGE        (A)" {
                    diffraction.wavelengths = std::str::from_utf8(value)
                        .map_err(|error| PdbRemark200Error::TextNotUtf8 {
                            valid_up_to: error.valid_up_to(),
                        })?
                        .to_owned();
                } else if key == b"MONOCHROMATOR" {
                    diffraction.monochromator = std::str::from_utf8(value)
                        .map_err(|error| PdbRemark200Error::TextNotUtf8 {
                            valid_up_to: error.valid_up_to(),
                        })?
                        .to_owned();
                } else if key == b"OPTICS" {
                    diffraction.optics = std::str::from_utf8(value)
                        .map_err(|error| PdbRemark200Error::TextNotUtf8 {
                            valid_up_to: error.valid_up_to(),
                        })?
                        .to_owned();
                } else if key == b"DETECTOR TYPE" {
                    diffraction.detector = std::str::from_utf8(value)
                        .map_err(|error| PdbRemark200Error::TextNotUtf8 {
                            valid_up_to: error.valid_up_to(),
                        })?
                        .to_owned();
                } else if key == b"DETECTOR MANUFACTURER" {
                    diffraction.detector_make = std::str::from_utf8(value)
                        .map_err(|error| PdbRemark200Error::TextNotUtf8 {
                            valid_up_to: error.valid_up_to(),
                        })?
                        .to_owned();
                } else if key == b"NUMBER OF UNIQUE REFLECTIONS" {
                    experiment.unique_reflections = remark3_int(raw_value)
                        .map_err(|_| PdbRemark200Error::IntegerOutsideSourceDefinedRange)?;
                } else if key == b"RESOLUTION RANGE HIGH      (A)" {
                    experiment.reflections.resolution_high = read_double(raw_value);
                } else if key == b"RESOLUTION RANGE LOW       (A)" {
                    experiment.reflections.resolution_low = read_double(raw_value);
                } else if key == b"COMPLETENESS FOR RANGE     (%)" {
                    experiment.reflections.completeness = read_double(raw_value);
                } else if key == b"DATA REDUNDANCY" {
                    experiment.reflections.redundancy = read_double(raw_value);
                } else if key == b"R MERGE                    (I)" {
                    experiment.reflections.r_merge = read_double(raw_value);
                } else if key == b"R SYM                      (I)" {
                    experiment.reflections.r_sym = read_double(raw_value);
                } else if key == b"<I/SIGMA(I)> FOR THE DATA SET" {
                    experiment.reflections.mean_i_over_sigma = read_double(raw_value);
                } else if key == b"REMARK" {
                    crystal.description = std::str::from_utf8(value)
                        .map_err(|error| PdbRemark200Error::TextNotUtf8 {
                            valid_up_to: error.valid_up_to(),
                        })?
                        .to_owned();
                    self.remark200_continuation = Some(crystal_index);
                } else if let Some(shell) = experiment.shells.last_mut() {
                    if key == b"HIGHEST RESOLUTION SHELL, RANGE HIGH (A)" {
                        shell.resolution_high = read_double(raw_value);
                    } else if key == b"HIGHEST RESOLUTION SHELL, RANGE LOW  (A)" {
                        shell.resolution_low = read_double(raw_value);
                    } else if key == b"COMPLETENESS FOR SHELL     (%)" {
                        shell.completeness = read_double(raw_value);
                    } else if key == b"DATA REDUNDANCY IN SHELL" {
                        shell.redundancy = read_double(raw_value);
                    } else if key == b"R MERGE FOR SHELL          (I)" {
                        shell.r_merge = read_double(raw_value);
                    } else if key == b"R SYM FOR SHELL            (I)" {
                        shell.r_sym = read_double(raw_value);
                    } else if key == b"<I/SIGMA(I)> FOR SHELL" {
                        shell.mean_i_over_sigma = read_double(raw_value);
                    }
                }
            }
        } else {
            if key == b"EXPERIMENTAL DETAILS" {
                let crystal_id = (self.metadata.crystals.len() + 1).to_string();
                let mut crystal = BioExperimentalCrystalInfo::default();
                crystal.id = crystal_id.clone();
                let mut diffraction = BioDiffractionInfo::default();
                diffraction.id = crystal_id.clone();
                let scattering_type = match (
                    line.get(8).copied().unwrap_or(0),
                    line.get(9).copied().unwrap_or(0),
                ) {
                    (b'0', b'0') => Some("x-ray"),
                    (b'3', b'0') => Some("neutron"),
                    (b'4', b'0') => Some("electron"),
                    _ => None,
                };
                if let Some(scattering_type) = scattering_type {
                    diffraction.scattering_type = scattering_type.to_owned();
                }
                crystal.diffractions.push(diffraction);
                let mut experiment = BioExperimentInfo::default();
                experiment.diffraction_ids.push(crystal_id);
                self.metadata.crystals.push(crystal);
                self.metadata.experiments.push(experiment);
            }
            if key == b"IN THE HIGHEST RESOLUTION SHELL." {
                if let Some(experiment) = self.metadata.experiments.last_mut() {
                    experiment.shells.push(Default::default());
                }
            }
        }

        Ok(())
    }

    fn remark_2_300_record(&mut self, remark: &str) -> Result<(), PdbRemarkMetadataError> {
        // Gemmi source: third_party/gemmi/src/pdb.cpp,
        // read_metadata_from_remarks.
        // Gemmi❗✔️:     if (remark.size() <= 11)
        // Gemmi❗✔️:       continue;
        // Gemmi❗✔️:     const char* line = remark.c_str();
        // Gemmi❗✔️:     int num = read_int(line + 7, 3);
        // Gemmi❗✔️:     switch (num) {
        // Gemmi❗✔️:       case 2:
        // Gemmi❗✔️:         if (st.resolution == 0.0 && std::strstr(line, "ANGSTROM"))
        // Gemmi❗✔️:           st.resolution = read_double(line + 23, 7);
        // Gemmi❗✔️:         break;
        // Gemmi❗✔️:       case 300:
        // Gemmi❗✔️:         if (!st.meta.remark_300_detail.empty()) {
        // Gemmi❗✔️:           st.meta.remark_300_detail += '\n';
        // Gemmi❗✔️:           st.meta.remark_300_detail += rtrim_str(remark.substr(11));
        // Gemmi❗✔️:         } else if (remark.compare(11, 7, "REMARK:") == 0) {
        // Gemmi❗✔️:           st.meta.remark_300_detail = trim_str(remark.substr(18));
        // Gemmi❗✔️:         }
        // Gemmi❗✔️:         break;
        // Gemmi❗✔️:     }
        // Behavior review: this method is called for each retained raw remark
        // in source order. The source byte-length gate and three-byte number
        // conversion precede dispatch. REMARK 2 searches only the C-string
        // prefix and reads its seven-byte field only when the existing
        // resolution is zero and `ANGSTROM` is present. A physically short
        // triggered field is source-undefined; it returns a typed boundary
        // error instead of reading beyond Rust storage. REMARK 300 preserves
        // the distinct first-row exact `REMARK:` test and later-row append
        // rule, trimming only the four bytes used by Gemmi `trim_str` and
        // `rtrim_str`. The refinement fallback runs after this row via the
        // existing source-anchored `update_remark3_resolution` helper.
        // Complexity review: remark-number parsing is fixed-width; the
        // `ANGSTROM` search and detail trim/copy are linear in this remark's
        // length, matching Gemmi's search/substr/append work. No scan over the
        // hierarchy or unrelated metadata is added.
        let bytes = remark.as_bytes();
        if bytes.len() <= 11 {
            return Ok(());
        }

        let number = read_int(&bytes[7..10])
            .ok_or(PdbRemarkMetadataError::IntegerOutsideSourceDefinedRange)?;
        match number {
            2 => {
                let c_end = bytes
                    .iter()
                    .position(|byte| *byte == 0)
                    .unwrap_or(bytes.len());
                let c_line = &bytes[..c_end];
                if self.source_state.resolution == 0.0
                    && c_line
                        .windows(b"ANGSTROM".len())
                        .any(|window| window == b"ANGSTROM")
                {
                    let field = bytes.get(23..30).ok_or(
                        PdbRemarkMetadataError::UndefinedResolutionField {
                            remark_length: bytes.len(),
                        },
                    )?;
                    self.source_state.resolution = read_double(field);
                }
            }
            300 => {
                let details = &mut self.metadata.remark_300_detail;
                if !details.is_empty() {
                    let text = remark
                        .get(11..)
                        .ok_or(PdbRemarkMetadataError::TextBoundaryNotUtf8 { offset: 11 })?;
                    details.push('\n');
                    details.push_str(
                        text.trim_end_matches(|ch| matches!(ch, ' ' | '\r' | '\n' | '\t')),
                    );
                } else if bytes.get(11..18) == Some(b"REMARK:".as_slice()) {
                    let text = remark
                        .get(18..)
                        .ok_or(PdbRemarkMetadataError::TextBoundaryNotUtf8 { offset: 18 })?;
                    *details = text
                        .trim_matches(|ch| matches!(ch, ' ' | '\r' | '\n' | '\t'))
                        .to_owned();
                }
            }
            _ => {}
        }

        self.update_remark3_resolution();
        Ok(())
    }

    fn remark_350_record(&mut self, remark: &str) -> Result<bool, PdbRemark350Error> {
        // Gemmi source: third_party/gemmi/src/pdb.cpp,
        // read_metadata_from_remarks.
        // Gemmi❗✔️:       case 350: {
        // Gemmi❗✔️:         const char* colon = std::strchr(line+11, ':');
        // Gemmi❗✔️:         if (colon == line+22 && starts_with(line+11, "BIOMOLECULE")) {
        // Gemmi❗✔️:           st.assemblies.emplace_back(read_string(line+23, 20));
        // Gemmi❗✔️:           continue;
        // Gemmi❗✔️:         }
        // Gemmi❗✔️:         if (st.assemblies.empty())
        // Gemmi❗✔️:           continue;
        // Gemmi❗✔️:         Assembly& assembly = st.assemblies.back();
        // Gemmi❗✔️:         auto r350_key = [&](int cpos, const char* text) {
        // Gemmi❗✔️:           return colon == line + cpos && starts_with(line+11, text);
        // Gemmi❗✔️:         };
        // Gemmi❗✔️:         if (starts_with(line+11, "  BIOMT")) {
        // Gemmi❗✔️:           if (read_matrix(matrix, line+13, remark.size()-13) == 3)
        // Gemmi❗✔️:             if (!assembly.generators.empty()) {
        // Gemmi❗✔️:               auto& opers = assembly.generators.back().operators;
        // Gemmi❗✔️:               opers.emplace_back();
        // Gemmi❗✔️:               opers.back().name = read_string(line+20, 3);
        // Gemmi❗✔️:               opers.back().transform = matrix;
        // Gemmi❗✔️:               matrix.set_identity();
        // Gemmi❗✔️:             }
        // Gemmi❗✔️:         } else if (r350_key(44, "AUTHOR DETERMINED")) {
        // Gemmi❗✔️:           assembly.author_determined = true;
        // Gemmi❗✔️:           assembly.oligomeric_details = read_string(line+45, 35);
        // Gemmi❗✔️:         } else if (r350_key(51, "SOFTWARE DETERMINED")) {
        // Gemmi❗✔️:           assembly.software_determined = true;
        // Gemmi❗✔️:           assembly.oligomeric_details = read_string(line+52, 28);
        // Gemmi❗✔️:         } else if (r350_key(24, "SOFTWARE USED")) {
        // Gemmi❗✔️:           assembly.software_name = read_string(line+25, 55);
        // Gemmi❗✔️:         } else if (r350_key(36, "TOTAL BURIED SURFACE AREA")) {
        // Gemmi❗✔️:           assembly.absa = read_double(line+37, 12);
        // Gemmi❗✔️:         } else if (r350_key(38, "SURFACE AREA OF THE COMPLEX")) {
        // Gemmi❗✔️:           assembly.ssa = read_double(line+39, 12);
        // Gemmi❗✔️:         } else if (r350_key(40, "CHANGE IN SOLVENT FREE ENERGY")) {
        // Gemmi❗✔️:           assembly.more = read_double(line+41, 12);
        // Gemmi❗✔️:         } else if (r350_key(40, "APPLY THE FOLLOWING TO CHAINS") ||
        // Gemmi❗✔️:                    r350_key(40, "                   AND CHAINS")) {
        // Gemmi❗✔️:           if (line[11] == 'A') // first line - APPLY ...
        // Gemmi❗✔️:             assembly.generators.emplace_back();
        // Gemmi❗✔️:           else if (assembly.generators.empty())
        // Gemmi❗✔️:             continue;
        // Gemmi❗✔️:           split_str_into_multi(read_string(line+41, 39), ", ",
        // Gemmi❗✔️:                                assembly.generators.back().chains);
        // Gemmi❗✔️:         }
        // Gemmi❗✔️:       }
        // Behavior review: this consumes one raw remark at a time while
        // retaining the source loop's `continue` result as `true` (skip the
        // shared post-row fallback) for short rows, assembly creation,
        // missing current assembly, and continuation rows without a generator.
        // BIOMT matrix rows accumulate independently in source order; only a
        // row numbered 3 with a current generator emits an operator and resets
        // the scratch transform. PDB assembly metadata keeps Gemmi's empty,
        // false, zero, `SpecialKind::NA`, and NaN constructor defaults. Chain
        // references and duplicates are retained verbatim for canonical
        // BioStructureData validation; no chemistry or chain inference is done.
        // Source-undefined fixed-width floating reads and unrepresentable
        // UTF-8 text return typed errors before committing row state.
        // Complexity review: prefix/colon scans are O(record length), each
        // fixed matrix row performs four bounded numeric reads, and chain
        // splitting is linear in the fixed 39-byte field with allocations
        // proportional to retained chain tokens, matching Gemmi's scans and
        // vector/string work. No hierarchy-wide search is added here.

        let bytes = remark.as_bytes();
        if bytes.len() <= 11 {
            return Ok(true);
        }

        let number =
            read_int(&bytes[7..10]).ok_or(PdbRemark350Error::IntegerOutsideSourceDefinedRange)?;
        if number != 350 {
            return Ok(false);
        }

        let tail = &bytes[11..];
        let c_tail_end = tail
            .iter()
            .position(|byte| *byte == 0)
            .unwrap_or(tail.len());
        let c_tail = &tail[..c_tail_end];
        let colon = c_tail
            .iter()
            .position(|byte| *byte == b':')
            .map(|offset| offset + 11);

        // The source forms `line + 22` for the first pointer comparison. The
        // C-string allocation includes its terminal NUL, so an offset through
        // `len + 1` is valid; larger offsets have no source-defined result.
        if 22 > bytes.len().saturating_add(1) {
            return Err(PdbRemark350Error::SourceUndefinedAccess {
                offset: 22,
                width: 0,
                record_length: bytes.len(),
            });
        }

        if colon == Some(22) && c_tail.starts_with(b"BIOMOLECULE") {
            let name_offset = 23;
            if name_offset > bytes.len() {
                return Err(PdbRemark350Error::SourceUndefinedAccess {
                    offset: name_offset,
                    width: 20,
                    record_length: bytes.len(),
                });
            }
            let name_end = (name_offset + 20).min(bytes.len());
            let name = String::from_utf8(read_string(&bytes[name_offset..name_end]).to_vec())
                .map_err(|error| PdbRemark350Error::TextNotUtf8 {
                    field_offset: name_offset,
                    valid_up_to: error.utf8_error().valid_up_to(),
                })?;
            self.assemblies.push(BioAssembly::new(
                name,
                false,
                false,
                BioAssemblySpecialKind::NotApplicable,
                0,
                String::new(),
                String::new(),
                f64::NAN,
                f64::NAN,
                f64::NAN,
                Vec::new(),
            ));
            return Ok(true);
        }

        if self.assemblies.is_empty() {
            return Ok(true);
        }

        let r350_key = |column: usize, text: &[u8]| {
            if column > bytes.len().saturating_add(1) {
                return Err(PdbRemark350Error::SourceUndefinedAccess {
                    offset: column,
                    width: 0,
                    record_length: bytes.len(),
                });
            }
            Ok(colon == Some(column) && c_tail.starts_with(text))
        };

        if c_tail.starts_with(b"  BIOMT") {
            let matrix_input_length = bytes.len() - 13;
            if matrix_input_length >= 46 {
                let row_number = i32::from(bytes[18]) - i32::from(b'0');
                if (1..=3).contains(&row_number) && matrix_input_length < 55 {
                    return Err(PdbRemark350Error::SourceUndefinedAccess {
                        offset: 58,
                        width: 10,
                        record_length: bytes.len(),
                    });
                }
            }
            let mut next_matrix = self.remark350_matrix;
            let row = read_matrix(&mut next_matrix, &bytes[13..]);
            let has_generator = self
                .assemblies
                .last()
                .is_some_and(|assembly| !assembly.generators.is_empty());
            if row == 3 && has_generator {
                let operator_name_offset = 20;
                let operator_name = String::from_utf8(
                    read_string(&bytes[operator_name_offset..operator_name_offset + 3]).to_vec(),
                )
                .map_err(|error| PdbRemark350Error::TextNotUtf8 {
                    field_offset: operator_name_offset,
                    valid_up_to: error.utf8_error().valid_up_to(),
                })?;
                let assembly = self
                    .assemblies
                    .last_mut()
                    .ok_or(PdbRemark350Error::StateInvariantLost)?;
                let generator = assembly
                    .generators
                    .last_mut()
                    .ok_or(PdbRemark350Error::StateInvariantLost)?;
                generator.operators.push(BioAssemblyOperator::new(
                    Some(operator_name),
                    None,
                    next_matrix,
                ));
                self.remark350_matrix = BioTransform::identity();
            } else {
                self.remark350_matrix = next_matrix;
            }
            return Ok(false);
        }

        if r350_key(44, b"AUTHOR DETERMINED")? {
            let details_offset = 45;
            if details_offset > bytes.len() {
                return Err(PdbRemark350Error::SourceUndefinedAccess {
                    offset: details_offset,
                    width: 35,
                    record_length: bytes.len(),
                });
            }
            let details_end = (details_offset + 35).min(bytes.len());
            let details =
                String::from_utf8(read_string(&bytes[details_offset..details_end]).to_vec())
                    .map_err(|error| PdbRemark350Error::TextNotUtf8 {
                        field_offset: details_offset,
                        valid_up_to: error.utf8_error().valid_up_to(),
                    })?;
            let assembly = self
                .assemblies
                .last_mut()
                .ok_or(PdbRemark350Error::StateInvariantLost)?;
            assembly.author_determined = true;
            assembly.oligomeric_details = details;
        } else if r350_key(51, b"SOFTWARE DETERMINED")? {
            let details_offset = 52;
            if details_offset > bytes.len() {
                return Err(PdbRemark350Error::SourceUndefinedAccess {
                    offset: details_offset,
                    width: 28,
                    record_length: bytes.len(),
                });
            }
            let details_end = (details_offset + 28).min(bytes.len());
            let details =
                String::from_utf8(read_string(&bytes[details_offset..details_end]).to_vec())
                    .map_err(|error| PdbRemark350Error::TextNotUtf8 {
                        field_offset: details_offset,
                        valid_up_to: error.utf8_error().valid_up_to(),
                    })?;
            let assembly = self
                .assemblies
                .last_mut()
                .ok_or(PdbRemark350Error::StateInvariantLost)?;
            assembly.software_determined = true;
            assembly.oligomeric_details = details;
        } else if r350_key(24, b"SOFTWARE USED")? {
            let name_offset = 25;
            if name_offset > bytes.len() {
                return Err(PdbRemark350Error::SourceUndefinedAccess {
                    offset: name_offset,
                    width: 55,
                    record_length: bytes.len(),
                });
            }
            let name_end = (name_offset + 55).min(bytes.len());
            let software_name = String::from_utf8(
                read_string(&bytes[name_offset..name_end]).to_vec(),
            )
            .map_err(|error| PdbRemark350Error::TextNotUtf8 {
                field_offset: name_offset,
                valid_up_to: error.utf8_error().valid_up_to(),
            })?;
            self.assemblies
                .last_mut()
                .ok_or(PdbRemark350Error::StateInvariantLost)?
                .software_name = software_name;
        } else if r350_key(36, b"TOTAL BURIED SURFACE AREA")? {
            let field = bytes
                .get(37..49)
                .ok_or(PdbRemark350Error::SourceUndefinedAccess {
                    offset: 37,
                    width: 12,
                    record_length: bytes.len(),
                })?;
            self.assemblies
                .last_mut()
                .ok_or(PdbRemark350Error::StateInvariantLost)?
                .buried_surface_area = read_double(field);
        } else if r350_key(38, b"SURFACE AREA OF THE COMPLEX")? {
            let field = bytes
                .get(39..51)
                .ok_or(PdbRemark350Error::SourceUndefinedAccess {
                    offset: 39,
                    width: 12,
                    record_length: bytes.len(),
                })?;
            self.assemblies
                .last_mut()
                .ok_or(PdbRemark350Error::StateInvariantLost)?
                .surface_area = read_double(field);
        } else if r350_key(40, b"CHANGE IN SOLVENT FREE ENERGY")? {
            let field = bytes
                .get(41..53)
                .ok_or(PdbRemark350Error::SourceUndefinedAccess {
                    offset: 41,
                    width: 12,
                    record_length: bytes.len(),
                })?;
            self.assemblies
                .last_mut()
                .ok_or(PdbRemark350Error::StateInvariantLost)?
                .solvent_free_energy_change = read_double(field);
        } else if r350_key(40, b"APPLY THE FOLLOWING TO CHAINS")?
            || r350_key(40, b"                   AND CHAINS")?
        {
            let apply = bytes[11] == b'A';
            if !apply
                && self
                    .assemblies
                    .last()
                    .is_some_and(|assembly| assembly.generators.is_empty())
            {
                return Ok(true);
            }
            let chains_offset = 41;
            if chains_offset > bytes.len() {
                return Err(PdbRemark350Error::SourceUndefinedAccess {
                    offset: chains_offset,
                    width: 39,
                    record_length: bytes.len(),
                });
            }
            let chains_end = (chains_offset + 39).min(bytes.len());
            let chain_text =
                String::from_utf8(read_string(&bytes[chains_offset..chains_end]).to_vec())
                    .map_err(|error| PdbRemark350Error::TextNotUtf8 {
                        field_offset: chains_offset,
                        valid_up_to: error.utf8_error().valid_up_to(),
                    })?;
            let chains: Vec<String> = chain_text
                .split(|character| character == ',' || character == ' ')
                .filter(|chain| !chain.is_empty())
                .map(str::to_owned)
                .collect();
            let assembly = self
                .assemblies
                .last_mut()
                .ok_or(PdbRemark350Error::StateInvariantLost)?;
            if apply {
                assembly.generators.push(BioAssemblyGenerator::default());
            }
            assembly
                .generators
                .last_mut()
                .ok_or(PdbRemark350Error::StateInvariantLost)?
                .chains
                .extend(chains);
        }

        Ok(false)
    }

    fn add_software(
        &mut self,
        classification: BioSoftwareClassification,
        name: &[u8],
    ) -> Result<(), PdbRemark3Error> {
        // Gemmi source: third_party/gemmi/src/pdb.cpp, add_software.
        // Gemmi❗✔️:   for (size_t start = 0, end = 0; end != std::string::npos; start = end + 1) {
        // Gemmi❗✔️:     end = name.find(',', start);
        // Gemmi❗✔️:     while (end != std::string::npos &&
        // Gemmi❗✔️:            name[end+1] == ' ' && is_digit(name[end+2]))
        // Gemmi❗✔️:       end = name.find(',', end + 1);
        // Gemmi❗✔️:     meta.software.emplace_back();
        // Gemmi❗✔️:     SoftwareItem& item = meta.software.back();
        // Gemmi❗✔️:     item.name = trim_str(name.substr(start, end - start));
        // Gemmi❗✔️:     size_t sep = item.name.find(' ');
        // Gemmi❗✔️:     if (sep != std::string::npos) {
        // Gemmi❗✔️:       size_t ver_start = item.name.find_first_not_of(" (", sep + 1);
        // Gemmi❗✔️:       if (ver_start == std::string::npos) {
        // Gemmi❗✔️:         item.name.resize(sep);
        // Gemmi❗✔️:         item.classification = type;
        // Gemmi❗✔️:         continue;
        // Gemmi❗✔️:       }
        // Gemmi❗✔️:       item.version = item.name.substr(ver_start);
        // Gemmi❗✔️:       item.name.resize(sep);
        // Gemmi❗✔️:       if (!item.version.empty() && item.version.back() == ')') {
        // Gemmi❗✔️:         size_t open_br = item.version.find('(');
        // Gemmi❗✔️:         if (open_br == std::string::npos) {
        // Gemmi❗✔️:           item.version.pop_back();
        // Gemmi❗✔️:         } else if (open_br + 11 == item.version.size() ||
        // Gemmi❗✔️:                    open_br + 13 == item.version.size()) {
        // Gemmi❗✔️:           item.date = pdb_date_format_to_iso(item.version.substr(open_br + 1));
        // Gemmi❗✔️:           if (item.date.size() == 10 && item.date[5] != 'x') {
        // Gemmi❗✔️:             size_t last = item.version.find_last_not_of(' ', open_br - 1);
        // Gemmi❗✔️:             item.version.resize(last + 1);
        // Gemmi❗✔️:           } else {
        // Gemmi❗✔️:             item.date.clear();
        // Gemmi❗✔️:           }
        // Gemmi❗✔️:         }
        // Gemmi❗✔️:       }
        // Gemmi❗✔️:       if (istarts_with(item.version, "version "))
        // Gemmi❗✔️:         item.version.erase(0, 8);
        // Gemmi❗✔️:     }
        // Gemmi❗✔️:     item.classification = type;
        // Gemmi❗✔️:   }
        // The source helper's comma rule, ASCII field trimming, version/date
        // extraction and case-insensitive `version ` removal are retained for
        // this refinement PROGRAM call. String construction and vector append
        // remain source-ordered; no guessed software-name normalization is used.
        let name = std::str::from_utf8(name).map_err(|error| PdbRemark3Error::TextNotUtf8 {
            valid_up_to: error.valid_up_to(),
        })?;
        let mut start = 0;
        loop {
            let mut end = name[start..].find(',').map(|offset| start + offset);
            while let Some(comma) = end {
                let bytes = name.as_bytes();
                if bytes.get(comma + 1) == Some(&b' ')
                    && bytes.get(comma + 2).is_some_and(u8::is_ascii_digit)
                {
                    end = name[comma + 1..].find(',').map(|offset| comma + 1 + offset);
                } else {
                    break;
                }
            }
            let segment_end = end.unwrap_or(name.len());
            let parsed_name = gemmi_pdb_trim_text(name[start..segment_end].as_bytes());
            let parsed_name =
                std::str::from_utf8(parsed_name).map_err(|error| PdbRemark3Error::TextNotUtf8 {
                    valid_up_to: error.valid_up_to(),
                })?;
            let mut item = BioSoftwareItem {
                name: parsed_name.to_owned(),
                version: String::new(),
                date: String::new(),
                description: String::new(),
                contact_author: String::new(),
                contact_author_email: String::new(),
                classification,
            };
            if let Some(separator) = item.name.find(' ') {
                let version_start = item.name[separator + 1..]
                    .find(|character| character != ' ' && character != '(')
                    .map(|offset| separator + 1 + offset);
                let Some(version_start) = version_start else {
                    item.name.truncate(separator);
                    self.metadata.software.push(item);
                    if let Some(comma) = end {
                        start = comma + 1;
                        continue;
                    }
                    return Ok(());
                };
                item.version = item.name[version_start..].to_owned();
                item.name.truncate(separator);
                if item.version.ends_with(')') {
                    if let Some(open_bracket) = item.version.find('(') {
                        if open_bracket + 11 == item.version.len()
                            || open_bracket + 13 == item.version.len()
                        {
                            item.date =
                                pdb_date_format_to_iso(item.version[open_bracket + 1..].as_bytes());
                            if item.date.len() == 10 && item.date.as_bytes()[5] != b'x' {
                                let end = item.version[..open_bracket]
                                    .rfind(|character| character != ' ')
                                    .map_or(0, |index| index + 1);
                                item.version.truncate(end);
                            } else {
                                item.date.clear();
                            }
                        }
                    } else {
                        item.version.pop();
                    }
                }
                if item
                    .version
                    .get(..8)
                    .is_some_and(|prefix| prefix.eq_ignore_ascii_case("version "))
                {
                    item.version.drain(..8);
                }
            }
            self.metadata.software.push(item);
            let Some(comma) = end else {
                return Ok(());
            };
            start = comma + 1;
        }
    }

    fn update_remark3_resolution(&mut self) {
        // Gemmi source: third_party/gemmi/src/pdb.cpp,
        // read_metadata_from_remarks.
        // Gemmi❗✔️:     if (st.resolution == 0.0) {
        // Gemmi❗✔️:       for (const RefinementInfo& ref_info : st.meta.refinement)
        // Gemmi❗✔️:         if (!std::isnan(ref_info.resolution_high) && ref_info.resolution_high != 0.) {
        // Gemmi❗✔️:           st.resolution = ref_info.resolution_high;
        // Gemmi❗✔️:           break;
        // Gemmi❗✔️:         }
        // Gemmi❗✔️:     }
        // Preserve an existing REMARK 2 result, otherwise choose the first
        // source-ordered refinement value that is neither NaN nor zero.
        if self.source_state.resolution == 0.0 {
            if let Some(resolution) = self
                .metadata
                .refinement
                .iter()
                .map(|refinement| refinement.basic.resolution_high)
                .find(|resolution| !resolution.is_nan() && *resolution != 0.0)
            {
                self.source_state.resolution = resolution;
            }
        }
    }

    fn modres_record(
        &mut self,
        source_line_buffer: &[u8; 122],
        line_len: usize,
    ) -> Result<bool, PdbModResError> {
        // Gemmi source: third_party/gemmi/src/pdb.cpp,
        // populate_structure_from_pdb_stream.
        // Gemmi❗✔️:     } else if (is_record_type4(line, "MODRES")) {
        // Gemmi❗✔️:       ModRes modres;
        // Gemmi❗✔️:       modres.chain_name = read_string(line + 15, 2);
        // Gemmi❗✔️:       modres.res_id = read_res_id(line + 18, line + 12);
        // Gemmi❗✔️:       modres.parent_comp_id = read_string(line + 24, 3);
        // Gemmi❗✔️:       if (len >= 30)
        // Gemmi❗✔️:         // this field is named comment in PDB spec, but details in mmCIF
        // Gemmi❗✔️:         modres.details = read_string(line + 29, 41);
        // Gemmi❗✔️:       // Refmac's extension: 73-80 mod_id
        // Gemmi❗✔️:       // Check for spaces to make sure it's not an overflowed comment
        // Gemmi❗✔️:       if (len >= 73 && line[70] == ' ' && line[71] == ' ')
        // Gemmi❗✔️:         modres.mod_id = read_string(line + 72, 8);
        // Gemmi❗✔️:       st.mod_residues.push_back(modres);
        // Gemmi❗✔️:     }
        // `read_string` and `read_res_id` are the existing unique fixed-field
        // owners. The record's optional fields use physical `len` gates even
        // though the backing buffer is NUL-padded, so short records retain
        // default-empty details/mod_id without reading beyond the source.
        // Behavior review: field offsets, optional gates, source order and
        // duplicate retention follow Gemmi. The canonical chain/residue types
        // enforce their already-approved ASCII boundary; invalid UTF-8 in the
        // three canonical String fields is returned as a typed error before
        // appending, so a failed record cannot partially mutate this state.
        // Complexity review: five bounded field reads and one amortized Vec
        // append per matched row; no scans over prior rows or hierarchy state.
        if !gemmi_record_type4(source_line_buffer, 0, *b"MODR") {
            return Ok(false);
        }

        let raw_chain_name = read_string(&source_line_buffer[15..17]);
        let chain_name = PdbChainId::from_ascii(raw_chain_name).ok_or_else(|| {
            PdbModResError::ChainNameNotRepresentable {
                width: raw_chain_name.len(),
                non_ascii: raw_chain_name.iter().copied().find(|byte| !byte.is_ascii()),
            }
        })?;

        let sequence_field: &[u8; 5] = source_line_buffer[18..23]
            .try_into()
            .expect("fixed PDB line buffer contains the MODRES sequence field");
        let residue_name_field: &[u8; 3] = source_line_buffer[12..15]
            .try_into()
            .expect("fixed PDB line buffer contains the MODRES residue name");
        let res_id = read_res_id(sequence_field, residue_name_field)
            .map_err(PdbModResError::ResidueAddress)?;
        let parent_comp_id = String::from_utf8(read_string(&source_line_buffer[24..27]).to_vec())
            .map_err(|error| PdbModResError::TextFieldNotUtf8 {
            field_offset: 24,
            valid_up_to: error.utf8_error().valid_up_to(),
        })?;
        let details = if line_len >= 30 {
            String::from_utf8(read_string(&source_line_buffer[29..70]).to_vec()).map_err(
                |error| PdbModResError::TextFieldNotUtf8 {
                    field_offset: 29,
                    valid_up_to: error.utf8_error().valid_up_to(),
                },
            )?
        } else {
            String::new()
        };
        let mod_id =
            if line_len >= 73 && source_line_buffer[70] == b' ' && source_line_buffer[71] == b' ' {
                String::from_utf8(read_string(&source_line_buffer[72..80]).to_vec()).map_err(
                    |error| PdbModResError::TextFieldNotUtf8 {
                        field_offset: 72,
                        valid_up_to: error.utf8_error().valid_up_to(),
                    },
                )?
            } else {
                String::new()
            };

        self.mod_residues.push(BioModRes {
            chain_name,
            res_id,
            parent_comp_id,
            mod_id,
            details,
        });
        Ok(true)
    }

    fn helix_record(
        &mut self,
        source_line_buffer: &[u8; 122],
        line_len: usize,
    ) -> Result<bool, PdbHelixError> {
        // Gemmi❗✔️:     } else if (is_record_type4(line, "HELIX")) {
        // Gemmi❗✔️:       if (len < 40)
        // Gemmi❗✔️:         continue;
        // Gemmi❗✔️:       Helix helix;
        // Gemmi❗✔️:       helix.start.chain_name = read_string(line+18, 2);
        // Gemmi❗✔️:       helix.start.res_id = read_res_id(line+21, line+15);
        // Gemmi❗✔️:       helix.end.chain_name = read_string(line+30, 2);
        // Gemmi❗✔️:       helix.end.res_id = read_res_id(line+33, line+27);
        // Gemmi❗✔️:       helix.set_helix_class_as_int(read_int(line+38, 2));
        // Gemmi❗✔️:       if (len > 72)
        // Gemmi❗✔️:         helix.length = read_int(line+72, 5);
        // Gemmi❗✔️:       st.helices.emplace_back(helix);
        // Gemmi❗✔️:     }
        // Behavior review: preserve the four-byte record dispatch and physical
        // length gates, source fixed-column offsets, trimmed source helpers,
        // default Unknown/-1 state, and setter semantics. The canonical BIO
        // chain/residue representation rejects bytes it cannot represent,
        // rather than truncating Gemmi's unrestricted strings; such rows
        // return a typed boundary error before appending any helix.
        // Complexity review: the branch performs a fixed number of bounded
        // field reads, constant-time canonical value construction, and one
        // amortized Vec append; it does not scan previously parsed helices.
        if !gemmi_record_type4(source_line_buffer, 0, *b"HELI") {
            return Ok(false);
        }
        if line_len < 40 {
            return Ok(true);
        }

        let start_chain_bytes = read_string(&source_line_buffer[18..20]);
        let start_chain = PdbChainId::from_ascii(start_chain_bytes).ok_or_else(|| {
            PdbHelixError::ChainNameNotRepresentable {
                field_offset: 18,
                width: start_chain_bytes.len(),
                non_ascii: start_chain_bytes
                    .iter()
                    .copied()
                    .find(|byte| !byte.is_ascii()),
            }
        })?;
        let start_sequence_field: &[u8; 5] = source_line_buffer[21..26]
            .try_into()
            .expect("fixed PDB line buffer contains the HELIX start sequence field");
        let start_residue_name: &[u8; 3] = source_line_buffer[15..18]
            .try_into()
            .expect("fixed PDB line buffer contains the HELIX start residue name");
        let start_residue = read_res_id(start_sequence_field, start_residue_name)
            .map_err(PdbHelixError::ResidueAddress)?;

        let end_chain_bytes = read_string(&source_line_buffer[30..32]);
        let end_chain = PdbChainId::from_ascii(end_chain_bytes).ok_or_else(|| {
            PdbHelixError::ChainNameNotRepresentable {
                field_offset: 30,
                width: end_chain_bytes.len(),
                non_ascii: end_chain_bytes
                    .iter()
                    .copied()
                    .find(|byte| !byte.is_ascii()),
            }
        })?;
        let end_sequence_field: &[u8; 5] = source_line_buffer[33..38]
            .try_into()
            .expect("fixed PDB line buffer contains the HELIX end sequence field");
        let end_residue_name: &[u8; 3] = source_line_buffer[27..30]
            .try_into()
            .expect("fixed PDB line buffer contains the HELIX end residue name");
        let end_residue = read_res_id(end_sequence_field, end_residue_name)
            .map_err(PdbHelixError::ResidueAddress)?;

        let class = read_int(&source_line_buffer[38..40]).ok_or(
            PdbHelixError::IntegerOutsideSourceDefinedRange {
                field_offset: 38,
                width: 2,
            },
        )?;
        let length = if line_len > 72 {
            read_int(&source_line_buffer[72..77]).ok_or(
                PdbHelixError::IntegerOutsideSourceDefinedRange {
                    field_offset: 72,
                    width: 5,
                },
            )?
        } else {
            -1
        };

        let mut helix = BioHelix::default();
        helix.start = cosmolkit_bio::AtomAddress::new(start_chain, start_residue, "", None);
        helix.end = cosmolkit_bio::AtomAddress::new(end_chain, end_residue, "", None);
        helix.set_helix_class_as_int(class);
        helix.length = length;
        self.helices.push(helix);
        Ok(true)
    }

    fn sheet_record(
        &mut self,
        source_line_buffer: &[u8; 122],
        line_len: usize,
    ) -> Result<bool, PdbSheetError> {
        // Gemmi❗✔️:     } else if (is_record_type4(line, "SHEET")) {
        // Gemmi❗✔️:       if (len < 40)
        // Gemmi❗✔️:         continue;
        // Gemmi❗✔️:       std::string sheet_id = read_string(line+11, 3);
        // Gemmi❗✔️:       Sheet& sheet = impl::find_or_add(st.sheets, sheet_id);
        // Gemmi❗✔️:       sheet.strands.emplace_back();
        // Gemmi❗✔️:       Sheet::Strand& strand = sheet.strands.back();
        // Gemmi❗✔️:       strand.start.chain_name = read_string(line+20, 2);
        // Gemmi❗✔️:       strand.start.res_id = read_res_id(line+22, line+17);
        // Gemmi❗✔️:       strand.end.chain_name = read_string(line+31, 2);
        // Gemmi❗✔️:       strand.end.res_id = read_res_id(line+33, line+28);
        // Gemmi❗✔️:       strand.sense = read_int(line+38, 2);
        // Gemmi❗✔️:       if (len > 67) {
        // Gemmi❗✔️:         // the SHEET record has no altloc for atoms of hydrogen bond
        // Gemmi❗✔️:         strand.hbond_atom2.atom_name = read_string(line+41, 4);
        // Gemmi❗✔️:         strand.hbond_atom2.chain_name = read_string(line+48, 2);
        // Gemmi❗✔️:         strand.hbond_atom2.res_id = read_res_id(line+50, line+45);
        // Gemmi❗✔️:         strand.hbond_atom1.atom_name = read_string(line+56, 4);
        // Gemmi❗✔️:         strand.hbond_atom1.chain_name = read_string(line+63, 2);
        // Gemmi❗✔️:         strand.hbond_atom1.res_id = read_res_id(line+65, line+60);
        // Gemmi❗✔️:       }
        // Gemmi❗✔️:     }
        // Gemmi source: third_party/gemmi/include/gemmi/model.hpp.
        // Gemmi❗✔️: template<typename T>
        // Gemmi❗✔️: auto get_id(const T& m) -> decltype(m.name) { return m.name; }
        // Gemmi❗✔️: template<typename Vec, typename S>
        // Gemmi❗✔️: auto find_iter_(Vec& vec, const S& name) {
        // Gemmi❗✔️:   return std::find_if(vec.begin(), vec.end(), [&name](const auto& m) { return get_id(m) == name; });
        // Gemmi❗✔️: }
        // Gemmi❗✔️: template<typename T, typename S>
        // Gemmi❗✔️: T* find_or_null(std::vector<T>& vec, const S& name) {
        // Gemmi❗✔️:   auto it = find_iter_(vec, name);
        // Gemmi❗✔️:   return it != vec.end() ? &*it : nullptr;
        // Gemmi❗✔️: }
        // Gemmi❗✔️: template<typename T, typename S>
        // Gemmi❗✔️: T& find_or_add(std::vector<T>& vec, const S& name) {
        // Gemmi❗✔️:   if (T* ret = find_or_null(vec, name))
        // Gemmi❗✔️:     return *ret;
        // Gemmi❗✔️:   vec.emplace_back(name);
        // Gemmi❗✔️:   return vec.back();
        // Gemmi❗✔️: }
        // Behavior review: retain first-seen sheet order and append each
        // source-order strand to the exact matching sheet id. Preserve the
        // physical `len < 40` and `len > 67` branches, every fixed field and
        // the source endpoint/H-bond assignment order. The owned BIO model
        // requires UTF-8 strings and bounded ASCII chain/residue vocabulary;
        // unrepresentable bytes and source-undefined integer overflow return
        // typed errors before mutating sheet state. Missing H-bond columns
        // retain default empty addresses; these addresses have no altloc.
        // Complexity review: source `find_or_add` and this `position` lookup
        // scan sheets linearly, O(S), followed by a constant number of fixed-
        // width reads and one amortized strand append; no full-row rescan.
        if !gemmi_record_type4(source_line_buffer, 0, *b"SHEE") {
            return Ok(false);
        }
        if line_len < 40 {
            return Ok(true);
        }

        let sheet_id_bytes = read_string(&source_line_buffer[11..14]);
        let sheet_id = String::from_utf8(sheet_id_bytes.to_vec()).map_err(|error| {
            PdbSheetError::TextFieldNotUtf8 {
                field_offset: 11,
                valid_up_to: error.utf8_error().valid_up_to(),
            }
        })?;

        let decode_chain = |field_offset: usize| {
            let bytes = read_string(&source_line_buffer[field_offset..field_offset + 2]);
            PdbChainId::from_ascii(bytes).ok_or_else(|| PdbSheetError::ChainNameNotRepresentable {
                field_offset,
                width: bytes.len(),
                non_ascii: bytes.iter().copied().find(|byte| !byte.is_ascii()),
            })
        };
        let decode_residue = |sequence_offset: usize, name_offset: usize| {
            let sequence_field: &[u8; 5] = source_line_buffer[sequence_offset..sequence_offset + 5]
                .try_into()
                .expect("fixed PDB line buffer contains a five-byte SHEET sequence field");
            let name_field: &[u8; 3] = source_line_buffer[name_offset..name_offset + 3]
                .try_into()
                .expect("fixed PDB line buffer contains a three-byte SHEET residue name");
            read_res_id(sequence_field, name_field).map_err(|source| {
                PdbSheetError::ResidueAddress {
                    sequence_offset,
                    source,
                }
            })
        };
        let decode_atom_name = |field_offset: usize| {
            String::from_utf8(
                read_string(&source_line_buffer[field_offset..field_offset + 4]).to_vec(),
            )
            .map_err(|error| PdbSheetError::TextFieldNotUtf8 {
                field_offset,
                valid_up_to: error.utf8_error().valid_up_to(),
            })
        };

        let start_chain = decode_chain(20)?;
        let start_residue = decode_residue(22, 17)?;
        let end_chain = decode_chain(31)?;
        let end_residue = decode_residue(33, 28)?;
        let sense = read_int(&source_line_buffer[38..40]).ok_or(
            PdbSheetError::IntegerOutsideSourceDefinedRange {
                field_offset: 38,
                width: 2,
            },
        )?;

        let (hbond_atom2, hbond_atom1) = if line_len > 67 {
            let atom2_name = decode_atom_name(41)?;
            let atom2_chain = decode_chain(48)?;
            let atom2_residue = decode_residue(50, 45)?;
            let atom1_name = decode_atom_name(56)?;
            let atom1_chain = decode_chain(63)?;
            let atom1_residue = decode_residue(65, 60)?;
            (
                AtomAddress::new(atom2_chain, atom2_residue, atom2_name, None),
                AtomAddress::new(atom1_chain, atom1_residue, atom1_name, None),
            )
        } else {
            (AtomAddress::default(), AtomAddress::default())
        };

        let start = AtomAddress::new(start_chain, start_residue, "", None);
        let end = AtomAddress::new(end_chain, end_residue, "", None);
        let strand = BioStrand::new(start, end, hbond_atom2, hbond_atom1, sense, String::new());

        let sheet_index = self
            .sheets
            .iter()
            .position(|sheet| sheet.name == sheet_id)
            .unwrap_or_else(|| {
                self.sheets.push(BioSheet::new(&sheet_id));
                self.sheets.len() - 1
            });
        self.sheets[sheet_index].strands.push(strand);
        Ok(true)
    }

    fn hetnam_record(&mut self, source_line_buffer: &[u8; 122], line_len: usize) -> bool {
        // Gemmi source: third_party/gemmi/src/pdb.cpp,
        // populate_structure_from_pdb_stream.
        // Gemmi❗✔️:     } else if (is_record_type4(line, "HETNAM")) {
        // Gemmi❗✔️:       if (len > 71 && line[70] == ' ') {
        // Gemmi❗✔️:         std::string full_code = read_string(line + 71, 8);
        // Gemmi❗✔️:         if (!full_code.empty())
        // Gemmi❗✔️:           st.shortened_ccd_codes.emplace_back(full_code, read_string(line + 11, 3));
        // Gemmi❗✔️:       }
        // A recognized HETNAM row is consumed even when its optional extension
        // gate is closed or the full code trims to empty. Raw bytes are retained
        // here, matching Gemmi's std::string; conversion to the bounded residue
        // name occurs only if restoration actually reaches a residue row.
        // Complexity review: four-byte type check, two fixed field reads and an
        // amortized append; no scan over existing aliases.
        if !gemmi_record_type4(source_line_buffer, 0, *b"HETN") {
            return false;
        }
        if line_len > 71 && source_line_buffer[70] == b' ' {
            let full_code = read_string(&source_line_buffer[71..79]);
            if !full_code.is_empty() {
                self.shortened_ccd_codes.push(PdbCcdAlias {
                    full_code: full_code.to_vec(),
                    short_code: read_string(&source_line_buffer[11..14]).to_vec(),
                });
            }
        }
        true
    }

    fn header_record(
        &mut self,
        source_line_buffer: &[u8; 122],
        line_len: usize,
    ) -> Result<bool, PdbHeaderError> {
        // Gemmi source: third_party/gemmi/src/pdb.cpp,
        // populate_structure_from_pdb_stream.
        // Gemmi❗✔️:     } else if (is_record_type4(line, "HEADER")) {
        // Gemmi❗✔️:       if (len > 50)
        // Gemmi❗✔️:         st.info["_struct_keywords.pdbx_keywords"] = rtrim_str(std::string(line+10, 40));
        // Gemmi❗✔️:       if (len > 59) {
        // Gemmi❗✔️:         std::string date = pdb_date_format_to_iso(std::string(line+50, 9));
        // Gemmi❗✔️:         if (!date.empty())
        // Gemmi❗✔️:           st.info["_pdbx_database_status.recvd_initial_deposition_date"] = date;
        // Gemmi❗✔️:       }
        // Gemmi❗✔️:       if (len > 66) {
        // Gemmi❗✔️:         std::string entry_id = rtrim_str(std::string(line+62, 4));
        // Gemmi❗✔️:         if (!entry_id.empty())
        // Gemmi❗✔️:           st.info["_entry.id"] = entry_id;
        // Gemmi❗✔️:       }
        // Gemmi❗✔️:     } else if (is_record_type4(line, "TITLE")) {
        // Gemmi❗✔️:       if (len > 10)
        // Gemmi❗✔️:         st.info["_struct.title"] += rtrim_str(std::string(line+10, len-10-1));
        // Gemmi❗✔️:     } else if (is_record_type4(line, "KEYWDS")) {
        // Gemmi❗✔️:       if (len > 10)
        // Gemmi❗✔️:         st.info["_struct_keywords.text"] += rtrim_str(std::string(line+10, len-10-1));
        // Gemmi❗✔️:     } else if (is_record_type4(line, "EXPDTA")) {
        // Gemmi❗✔️:       if (len > 10)
        // Gemmi❗✔️:         st.info["_exptl.method"] += trim_str(std::string(line+10, len-10-1));
        // Gemmi❗✔️:     } else if (is_record_type4(line, "AUTHOR") && len > 10) {
        // Gemmi❗✔️:       std::string last;
        // Gemmi❗✔️:       if (!st.meta.authors.empty()) {
        // Gemmi❗✔️:         last = st.meta.authors.back();
        // Gemmi❗✔️:         st.meta.authors.pop_back();
        // Gemmi❗✔️:       }
        // Gemmi❗✔️:       size_t prev_size = st.meta.authors.size();
        // Gemmi❗✔️:       const char* start = skip_blank(line+10);
        // Gemmi❗✔️:       const char* end = rtrim_cstr(start, line+len);
        // Gemmi❗✔️:       split_str_into(std::string(start, end), ',', st.meta.authors);
        // Gemmi❗✔️:       if (!last.empty() && st.meta.authors.size() > prev_size) {
        // Gemmi❗✔️:         // the spaces were trimmed, we may need a space between words
        // Gemmi❗✔️:         if (last.back() != '-' && last.back() != '.')
        // Gemmi❗✔️:           last += ' ';
        // Gemmi❗✔️:         st.meta.authors[prev_size].insert(0, last);
        // Gemmi❗✔️:       }
        // Gemmi❗✔️:     }
        // The source helpers `rtrim_str`, `trim_str`, `skip_blank`,
        // `rtrim_cstr`, and `split_str_into` are anchored in the byte helpers
        // below. The current `BioStructureSourceState::info` and
        // `BioMetadata::authors` are the canonical destinations; the older
        // audit's pre-metadata availability row is historical.
        // Behavior review: preserve Gemmi's record-prefix and length gates,
        // fixed HEADER slices, overwrite-versus-append behavior, exact ASCII
        // trim sets, comma empty fields, continuation merge order, and defer
        // author name normalization until parser finalization. A fixed-column
        // byte slice that is not valid UTF-8 returns a typed Rust-value
        // representation error instead of lossily changing the source bytes.
        // Complexity review: fixed HEADER slices and tags are O(1) apart from
        // their bounded 40-byte conversion; continuation work is O(line_len)
        // plus one amortized vector append per comma-delimited name. BTreeMap
        // updates retain the source std::map logarithmic key lookup.
        if gemmi_record_type4(source_line_buffer, 0, *b"HEAD") {
            let keywords = if line_len > 50 {
                Some(pdb_header_text(
                    gemmi_pdb_rtrim_text(&source_line_buffer[10..50]),
                    10,
                )?)
            } else {
                None
            };
            let date = if line_len > 59 {
                let value = pdb_date_format_to_iso(&source_line_buffer[50..59]);
                (!value.is_empty()).then_some(value)
            } else {
                None
            };
            let entry_id = if line_len > 66 {
                let value = pdb_header_text(gemmi_pdb_rtrim_text(&source_line_buffer[62..66]), 62)?;
                (!value.is_empty()).then_some(value)
            } else {
                None
            };

            if let Some(value) = keywords {
                self.source_state
                    .info
                    .insert("_struct_keywords.pdbx_keywords".to_owned(), value);
            }
            if let Some(value) = date {
                self.source_state.info.insert(
                    "_pdbx_database_status.recvd_initial_deposition_date".to_owned(),
                    value,
                );
            }
            if let Some(value) = entry_id {
                self.source_state.info.insert("_entry.id".to_owned(), value);
            }
            return Ok(true);
        }

        if gemmi_record_type4(source_line_buffer, 0, *b"TITL") {
            if line_len > 10 {
                let value = pdb_header_text(
                    gemmi_pdb_rtrim_text(&source_line_buffer[10..line_len - 1]),
                    10,
                )?;
                self.source_state
                    .info
                    .entry("_struct.title".to_owned())
                    .or_default()
                    .push_str(&value);
            }
            return Ok(true);
        }

        if gemmi_record_type4(source_line_buffer, 0, *b"KEYW") {
            if line_len > 10 {
                let value = pdb_header_text(
                    gemmi_pdb_rtrim_text(&source_line_buffer[10..line_len - 1]),
                    10,
                )?;
                self.source_state
                    .info
                    .entry("_struct_keywords.text".to_owned())
                    .or_default()
                    .push_str(&value);
            }
            return Ok(true);
        }

        if gemmi_record_type4(source_line_buffer, 0, *b"EXPD") {
            if line_len > 10 {
                let value = pdb_header_text(
                    gemmi_pdb_trim_text(&source_line_buffer[10..line_len - 1]),
                    10,
                )?;
                self.source_state
                    .info
                    .entry("_exptl.method".to_owned())
                    .or_default()
                    .push_str(&value);
            }
            return Ok(true);
        }

        if gemmi_record_type4(source_line_buffer, 0, *b"AUTH") && line_len > 10 {
            let names = gemmi_pdb_author_fields(source_line_buffer, line_len)?;
            let mut last = self.metadata.authors.pop().unwrap_or_default();
            let previous_len = self.metadata.authors.len();
            self.metadata.authors.extend(names);
            if !last.is_empty() && self.metadata.authors.len() > previous_len {
                if !matches!(last.as_bytes().last(), Some(b'-' | b'.')) {
                    last.push(' ');
                }
                self.metadata.authors[previous_len].insert_str(0, &last);
            }
            return Ok(true);
        }

        Ok(false)
    }

    fn finalize_author_names(&mut self) {
        // Gemmi source: third_party/gemmi/src/pdb.cpp,
        // populate_structure_from_pdb_stream finalization.
        // Gemmi❗✔️:   for (std::string& name : st.meta.authors)
        // Gemmi❗✔️:     change_author_name_format_to_mmcif(name);
        // `change_author_name_format_to_mmcif` is the existing anchored owner;
        // apply it only after all AUTHOR continuations have been accumulated.
        // Behavior review: finalization order matches Gemmi's post-scan loop.
        // Complexity review: one pass over retained author names, with the
        // source-shaped per-name normalization cost in the shared helper.
        for author in &mut self.metadata.authors {
            change_author_name_format_to_mmcif(author);
        }
    }

    fn read_metadata_from_remarks(&mut self) -> Result<(), PdbRemarkFinalizationError> {
        // Gemmi source: third_party/gemmi/src/pdb.cpp::read_metadata_from_remarks.
        // Gemmi❗✔️:   for (const std::string& remark : st.raw_remarks) {
        // Gemmi❗✔️:     if (remark.size() <= 11)
        // Gemmi❗✔️:       continue;
        // Gemmi❗✔️:     const char* line = remark.c_str();
        // Gemmi❗✔️:     int num = read_int(line + 7, 3);
        // Gemmi❗✔️:     switch (num) {
        // Gemmi❗✔️:       case 2:
        // Gemmi❗✔️:         if (st.resolution == 0.0 && std::strstr(line, "ANGSTROM"))
        // Gemmi❗✔️:           st.resolution = read_double(line + 23, 7);
        // Gemmi❗✔️:         break;
        // Gemmi❗✔️:       case 3:
        // Gemmi❗✔️:         read_remark3_line(line, st.meta, possibly_unfinished_remark3);
        // Gemmi❗✔️:         break;
        // Gemmi❗✔️:       case 200:
        // Gemmi❗✔️:       case 230:
        // Gemmi❗✔️:       case 240:
        // Gemmi❗✔️:         read_remark_200_230_240(line, st.meta, cr_desc);
        // Gemmi❗✔️:         break;
        // Gemmi❗✔️:       case 300:
        // Gemmi❗✔️:         if (!st.meta.remark_300_detail.empty()) {
        // Gemmi❗✔️:           st.meta.remark_300_detail += '\n';
        // Gemmi❗✔️:           st.meta.remark_300_detail += rtrim_str(remark.substr(11));
        // Gemmi❗✔️:         } else if (remark.compare(11, 7, "REMARK:") == 0) {
        // Gemmi❗✔️:           st.meta.remark_300_detail = trim_str(remark.substr(18));
        // Gemmi❗✔️:         }
        // Gemmi❗✔️:         break;
        // Gemmi❗✔️:       case 350: {
        // Gemmi❗✔️:         const char* colon = std::strchr(line+11, ':');
        // Gemmi❗✔️:         if (colon == line+22 && starts_with(line+11, "BIOMOLECULE")) {
        // Gemmi❗✔️:           st.assemblies.emplace_back(read_string(line+23, 20));
        // Gemmi❗✔️:           continue;
        // Gemmi❗✔️:         }
        // Gemmi❗✔️:         if (st.assemblies.empty())
        // Gemmi❗✔️:           continue;
        // Gemmi❗✔️:         Assembly& assembly = st.assemblies.back();
        // Gemmi❗✔️:         auto r350_key = [&](int cpos, const char* text) {
        // Gemmi❗✔️:           return colon == line + cpos && starts_with(line+11, text);
        // Gemmi❗✔️:         };
        // Gemmi❗✔️:         if (starts_with(line+11, "  BIOMT")) {
        // Gemmi❗✔️:           if (read_matrix(matrix, line+13, remark.size()-13) == 3)
        // Gemmi❗✔️:             if (!assembly.generators.empty()) {
        // Gemmi❗✔️:               auto& opers = assembly.generators.back().operators;
        // Gemmi❗✔️:               opers.emplace_back();
        // Gemmi❗✔️:               opers.back().name = read_string(line+20, 3);
        // Gemmi❗✔️:               opers.back().transform = matrix;
        // Gemmi❗✔️:               matrix.set_identity();
        // Gemmi❗✔️:             }
        // Gemmi❗✔️:         } else if (r350_key(44, "AUTHOR DETERMINED")) {
        // Gemmi❗✔️:           assembly.author_determined = true;
        // Gemmi❗✔️:           assembly.oligomeric_details = read_string(line+45, 35);
        // Gemmi❗✔️:         } else if (r350_key(51, "SOFTWARE DETERMINED")) {
        // Gemmi❗✔️:           assembly.software_determined = true;
        // Gemmi❗✔️:           assembly.oligomeric_details = read_string(line+52, 28);
        // Gemmi❗✔️:         } else if (r350_key(24, "SOFTWARE USED")) {
        // Gemmi❗✔️:           assembly.software_name = read_string(line+25, 55);
        // Gemmi❗✔️:         } else if (r350_key(36, "TOTAL BURIED SURFACE AREA")) {
        // Gemmi❗✔️:           assembly.absa = read_double(line+37, 12);
        // Gemmi❗✔️:         } else if (r350_key(38, "SURFACE AREA OF THE COMPLEX")) {
        // Gemmi❗✔️:           assembly.ssa = read_double(line+39, 12);
        // Gemmi❗✔️:         } else if (r350_key(40, "CHANGE IN SOLVENT FREE ENERGY")) {
        // Gemmi❗✔️:           assembly.more = read_double(line+41, 12);
        // Gemmi❗✔️:         } else if (r350_key(40, "APPLY THE FOLLOWING TO CHAINS") ||
        // Gemmi❗✔️:                    r350_key(40, "                   AND CHAINS")) {
        // Gemmi❗✔️:           if (line[11] == 'A')
        // Gemmi❗✔️:             assembly.generators.emplace_back();
        // Gemmi❗✔️:           else if (assembly.generators.empty())
        // Gemmi❗✔️:             continue;
        // Gemmi❗✔️:           split_str_into_multi(read_string(line+41, 39), ", ",
        // Gemmi❗✔️:                                assembly.generators.back().chains);
        // Gemmi❗✔️:         }
        // Gemmi❗✔️:       }
        // Gemmi❗✔️:     }
        // Gemmi❗✔️:     if (st.resolution == 0.0) {
        // Gemmi❗✔️:       for (const RefinementInfo& ref_info : st.meta.refinement)
        // Gemmi❗✔️:         if (!std::isnan(ref_info.resolution_high) && ref_info.resolution_high != 0.) {
        // Gemmi❗✔️:           st.resolution = ref_info.resolution_high;
        // Gemmi❗✔️:           break;
        // Gemmi❗✔️:         }
        // Gemmi❗✔️:     }
        // Gemmi❗✔️:   }
        // Behavior review: source order is preserved and the exact REMARK350
        // continue cases skip the shared resolution fallback. Existing row
        // helpers own their anchored field semantics; the dispatch applies the
        // common fallback only where those helpers do not already apply it.
        // Complexity review: one pass over retained remarks; fixed-width number
        // reads and source helper scans are per-row, with no raw-remark clone.
        let raw_remarks = std::mem::take(&mut self.source_state.raw_remarks);
        let result = (|| {
            for remark in &raw_remarks {
                let bytes = remark.as_bytes();
                if bytes.len() <= 11 {
                    continue;
                }
                let number = read_int(&bytes[7..10])
                    .ok_or(PdbRemarkFinalizationError::IntegerOutsideSourceDefinedRange)?;
                match number {
                    2 | 300 => self
                        .remark_2_300_record(remark)
                        .map_err(PdbRemarkFinalizationError::RemarkMetadata)?,
                    3 => self
                        .remark3_record(remark)
                        .map_err(PdbRemarkFinalizationError::Remark3)?,
                    200 | 230 | 240 => {
                        self.remark_200_230_240_record(remark)
                            .map_err(PdbRemarkFinalizationError::Remark200)?;
                        self.update_remark3_resolution();
                    }
                    350 => {
                        if !self
                            .remark_350_record(remark)
                            .map_err(PdbRemarkFinalizationError::Remark350)?
                        {
                            self.update_remark3_resolution();
                        }
                    }
                    _ => self.update_remark3_resolution(),
                }
            }
            Ok(())
        })();
        self.source_state.raw_remarks = raw_remarks;
        result
    }

    fn restore_full_ccd_codes(
        &mut self,
        hierarchy: &mut PdbHierarchyGrouping,
        connections: &mut [BioConnection],
        cispeps: &mut [BioCisPep],
    ) -> Result<(), PdbResidueKeyError> {
        // Gemmi source: third_party/gemmi/src/polyheur.cpp.
        // Gemmi❗✔️: void restore_full_ccd_codes(Structure& st) {
        // Gemmi❗✔️:   for (const auto& item : st.shortened_ccd_codes)
        // Gemmi❗✔️:     rename_residues(st, item.second, item.first);
        // Gemmi❗✔️:   st.shortened_ccd_codes.clear();
        // Gemmi❗✔️: }
        // Gemmi source helper closure: third_party/gemmi/include/gemmi/modify.hpp.
        // Gemmi❗✔️: inline void rename_residues(Structure& st, const std::string& old_name,
        // Gemmi❗✔️:                                            const std::string& new_name) {
        // Gemmi❗✔️:   auto update = [&](ResidueId& rid) {
        // Gemmi❗✔️:     if (rid.name == old_name)
        // Gemmi❗✔️:       rid.name = new_name;
        // Gemmi❗✔️:   };
        // Gemmi❗✔️:   process_addresses(st, [&](AtomAddress& aa) { update(aa.res_id); });
        // Gemmi❗✔️:   for (ModRes& modres : st.mod_residues)
        // Gemmi❗✔️:     update(modres.res_id);
        // Gemmi❗✔️:   for (Entity& ent : st.entities)
        // Gemmi❗✔️:     for (std::string& mon_ids : ent.full_sequence)
        // Gemmi❗✔️:       for (size_t start = 0;;) {
        // Gemmi❗✔️:         size_t end = mon_ids.find(',', start);
        // Gemmi❗✔️:         if (mon_ids.compare(start, end-start, old_name) == 0) {
        // Gemmi❗✔️:           mon_ids.replace(start, end-start, new_name);
        // Gemmi❗✔️:           if (end != std::string::npos)
        // Gemmi❗✔️:             end = start + new_name.size();
        // Gemmi❗✔️:         }
        // Gemmi❗✔️:         if (end == std::string::npos)
        // Gemmi❗✔️:           break;
        // Gemmi❗✔️:         start = end + 1;
        // Gemmi❗✔️:       }
        // Gemmi❗✔️:   for (Model& model : st.models)
        // Gemmi❗✔️:     for (Chain& chain : model.chains)
        // Gemmi❗✔️:       for (Residue& res : chain.residues)
        // Gemmi❗✔️:         update(res);
        // Gemmi❗✔️: }
        // Gemmi source helper closure: third_party/gemmi/include/gemmi/modify.hpp.
        // Gemmi❗✔️: template<typename Func>
        // Gemmi❗✔️: void process_addresses(Structure& st, Func func) {
        // Gemmi❗✔️:   for (Connection& con : st.connections) {
        // Gemmi❗✔️:     func(con.partner1);
        // Gemmi❗✔️:     func(con.partner2);
        // Gemmi❗✔️:   }
        // Gemmi❗✔️:   for (CisPep& cispep : st.cispeps) {
        // Gemmi❗✔️:     func(cispep.partner_c);
        // Gemmi❗✔️:     func(cispep.partner_n);
        // Gemmi❗✔️:   }
        // Gemmi❗✔️:   for (Helix& helix : st.helices) {
        // Gemmi❗✔️:     func(helix.start);
        // Gemmi❗✔️:     func(helix.end);
        // Gemmi❗✔️:   }
        // Gemmi❗✔️:   for (Sheet& sheet : st.sheets)
        // Gemmi❗✔️:     for (Sheet::Strand& strand : sheet.strands) {
        // Gemmi❗✔️:       func(strand.start);
        // Gemmi❗✔️:       func(strand.end);
        // Gemmi❗✔️:       func(strand.hbond_atom2);
        // Gemmi❗✔️:       func(strand.hbond_atom1);
        // Gemmi❗✔️:     }
        // Gemmi❗✔️: }
        // The source address closure now covers the represented connection,
        // cis-peptide, helix, sheet, MODRES, entity-sequence, and hierarchy
        // values. A source name outside the existing canonical address width
        // remains a structured model-boundary error, not an upstream parity
        // claim. Preflight keeps such a failure atomic.
        // Behavior review: aliases are applied in input order; preflight checks
        // every bounded address after the complete ordered rename chain so an
        // over-width/non-ASCII restored residue fails before changing any
        // staged state. The raw SEQRES byte slots use Gemmi's comma-token
        // replacement order. Alias state is cleared only after successful
        // restoration. Complexity review: for A aliases and R staged addresses,
        // preflight and commit are both O(A*R), matching the source's per-alias
        // full-structure passes; sequence replacement adds only linear scans of
        // the bounded PDB sequence fields.

        for residue in &self.mod_residues {
            mapped_ccd_residue_address(residue.res_id, &self.shortened_ccd_codes)?;
        }
        for connection in connections.iter() {
            mapped_ccd_atom_address(&connection.partner1, &self.shortened_ccd_codes)?;
            mapped_ccd_atom_address(&connection.partner2, &self.shortened_ccd_codes)?;
        }
        for cispep in cispeps.iter() {
            mapped_ccd_atom_address(&cispep.partner_c, &self.shortened_ccd_codes)?;
            mapped_ccd_atom_address(&cispep.partner_n, &self.shortened_ccd_codes)?;
        }
        for helix in &self.helices {
            mapped_ccd_atom_address(&helix.start, &self.shortened_ccd_codes)?;
            mapped_ccd_atom_address(&helix.end, &self.shortened_ccd_codes)?;
        }
        for sheet in &self.sheets {
            for strand in &sheet.strands {
                mapped_ccd_atom_address(&strand.start, &self.shortened_ccd_codes)?;
                mapped_ccd_atom_address(&strand.end, &self.shortened_ccd_codes)?;
                mapped_ccd_atom_address(&strand.hbond_atom2, &self.shortened_ccd_codes)?;
                mapped_ccd_atom_address(&strand.hbond_atom1, &self.shortened_ccd_codes)?;
            }
        }
        for model in &hierarchy.models {
            for chain in &model.chains {
                for residue in &chain.residues {
                    mapped_ccd_residue_address(residue.address, &self.shortened_ccd_codes)?;
                }
            }
        }

        for entity in &mut self.entities.entities {
            for sequence_slot in &mut entity.full_sequence {
                for alias in &self.shortened_ccd_codes {
                    rename_ccd_sequence_tokens(sequence_slot, &alias.short_code, &alias.full_code);
                }
            }
        }
        for residue in &mut self.mod_residues {
            residue.res_id = mapped_ccd_residue_address(residue.res_id, &self.shortened_ccd_codes)?;
        }
        for connection in connections {
            connection.partner1 =
                mapped_ccd_atom_address(&connection.partner1, &self.shortened_ccd_codes)?;
            connection.partner2 =
                mapped_ccd_atom_address(&connection.partner2, &self.shortened_ccd_codes)?;
        }
        for cispep in cispeps {
            cispep.partner_c =
                mapped_ccd_atom_address(&cispep.partner_c, &self.shortened_ccd_codes)?;
            cispep.partner_n =
                mapped_ccd_atom_address(&cispep.partner_n, &self.shortened_ccd_codes)?;
        }
        for helix in &mut self.helices {
            helix.start = mapped_ccd_atom_address(&helix.start, &self.shortened_ccd_codes)?;
            helix.end = mapped_ccd_atom_address(&helix.end, &self.shortened_ccd_codes)?;
        }
        for sheet in &mut self.sheets {
            for strand in &mut sheet.strands {
                strand.start = mapped_ccd_atom_address(&strand.start, &self.shortened_ccd_codes)?;
                strand.end = mapped_ccd_atom_address(&strand.end, &self.shortened_ccd_codes)?;
                strand.hbond_atom2 =
                    mapped_ccd_atom_address(&strand.hbond_atom2, &self.shortened_ccd_codes)?;
                strand.hbond_atom1 =
                    mapped_ccd_atom_address(&strand.hbond_atom1, &self.shortened_ccd_codes)?;
            }
        }
        for model in &mut hierarchy.models {
            for chain in &mut model.chains {
                for residue in &mut chain.residues {
                    residue.address =
                        mapped_ccd_residue_address(residue.address, &self.shortened_ccd_codes)?;
                }
            }
        }
        self.shortened_ccd_codes.clear();
        Ok(())
    }

    fn dbref_record(
        &mut self,
        source_line_buffer: &[u8; 122],
    ) -> Result<PdbDbRefAction, PdbDbRefError> {
        self.entities.dbref_record(source_line_buffer)
    }

    fn conect_record(
        &mut self,
        source_line_buffer: &[u8; 122],
        line_len: usize,
    ) -> Result<bool, PdbConectError> {
        // Gemmi source: third_party/gemmi/src/pdb.cpp,
        // populate_structure_from_pdb_stream.
        // Gemmi❗✔️:     } else if (is_record_type4(line, "CONECT")) {
        // Gemmi❗✔️:       int serial = read_serial(line+6);
        // Gemmi❗✔️:       if (len >= 11 && serial != 0) {
        // Gemmi❗✔️:         std::vector<int>& bonded_atoms = st.conect_map[serial];
        // Gemmi❗✔️:         int limit = std::min(27, (int)len - 1);
        // Gemmi❗✔️:         for (int offset = 11; offset <= limit; offset += 5) {
        // Gemmi❗✔️:           int n = read_serial(line+offset);
        // Gemmi❗✔️:           if (n != 0)
        // Gemmi❗✔️:             bonded_atoms.push_back(n);
        // Gemmi❗✔️:         }
        // Gemmi❗✔️:       }
        // Gemmi❗✔️:     }
        // Gemmi's `std::map<int, std::vector<int>>` maps to the canonical
        // BTreeMap serial metadata and each record appends to its per-serial
        // vector. A missing serial value indicates source-undefined signed
        // integer overflow in the existing fixed-width decoder and remains a
        // typed error, not a fabricated serial.
        // Behavior review: the source-prefix check, primary serial read before
        // the length gate, len>=11/nonzero-source gate, map-entry creation,
        // four possible five-column target offsets, zero filtering, duplicate
        // retention and encounter order all follow the pinned loop. No
        // hierarchy index or chemical bond is synthesized from source serials.
        // Complexity review: one BTreeMap entry lookup (O(log distinct serials))
        // and at most four fixed-width target decodes per record match the
        // source's ordered map and bounded loop; vector appends are amortized
        // O(1), with no full-map or atom scan.
        if !gemmi_record_type4(source_line_buffer, 0, *b"CONE") {
            return Ok(false);
        }

        let source_serial_field: &[u8; 5] = source_line_buffer[6..11]
            .try_into()
            .expect("fixed PDB line buffer contains the CONECT source serial");
        let source_serial = read_serial(source_serial_field)
            .ok_or(PdbConectError::SourceSerialOutsideDefinedRange { field_offset: 6 })?;
        if line_len < 11 || source_serial == 0 {
            return Ok(true);
        }

        let bonded_atoms = self
            .source_state
            .conect_map
            .entry(source_serial)
            .or_default();
        let limit = (line_len - 1).min(27);
        for offset in (11..=limit).step_by(5) {
            let target_serial_field: &[u8; 5] = source_line_buffer[offset..offset + 5]
                .try_into()
                .expect("fixed PDB line buffer contains the CONECT target field");
            let target_serial = read_serial(target_serial_field).ok_or(
                PdbConectError::SourceSerialOutsideDefinedRange {
                    field_offset: offset,
                },
            )?;
            if target_serial != 0 {
                bonded_atoms.push(target_serial);
            }
        }
        Ok(true)
    }

    fn record_line_checked(&mut self, line: &[u8]) -> Result<(), PdbReaderCause> {
        // Gemmi source: third_party/gemmi/src/pdb.cpp,
        // populate_structure_from_pdb_stream.
        // Gemmi❗✔️:     ++line_num;
        // Gemmi❗✔️:     if (options.check_non_ascii && st.non_ascii_line == 0)
        // Gemmi❗✔️:       for (size_t i = 0; i < len; ++i)
        // Gemmi❗✔️:         if (static_cast<unsigned char>(line[i]) >= 0x80) {
        // Gemmi❗✔️:           st.non_ascii_line = line_num;
        // Gemmi❗✔️:           break;
        // Gemmi❗✔️:         }
        // C++ signed-int overflow in line_num is undefined. Keep it outside the
        // claimed source-equivalent domain and return a typed condition to the
        // orchestration error owner instead of panicking or wrapping.
        self.line_number = self
            .line_number
            .checked_add(1)
            .ok_or(PdbReaderCause::LineNumberOutsideSourceIntRange)?;
        if self.check_non_ascii && self.source_state.non_ascii_line == 0 {
            if line.iter().any(|byte| *byte >= 0x80) {
                self.source_state.non_ascii_line = self.line_number;
            }
        }
        Ok(())
    }

    #[cfg(test)]
    fn record_line(&mut self, line: &[u8]) {
        self.record_line_checked(line)
            .expect("fixed regression lines remain within Gemmi's defined int range");
    }

    fn wrong_input_format(
        &self,
        line: &[u8],
        model_active: bool,
        source: &str,
    ) -> Option<PdbInputFormatError> {
        // Gemmi source: third_party/gemmi/src/pdb.cpp,
        // populate_structure_from_pdb_stream.
        // Gemmi❗✔️:     } else if (is_record_type4(line, "data")) {
        // Gemmi❗✔️:       if (line[4] == '_' && !model)
        // Gemmi❗✔️:         fail("Incorrect file format (perhaps it is cif not pdb?): " + source);
        // Gemmi❗✔️:     }
        // Gemmi❗✔️:     } else if (is_record_type4(line, "{\"da")) {
        // Gemmi❗✔️:       if (ialpha3_id(line+4) == ialpha3_id("ta_") && !model)
        // Gemmi❗✔️:       fail("Incorrect file format (perhaps it is mmJSON not pdb?): " + source);
        // Gemmi❗✔️:     }
        if model_active {
            return None;
        }
        if gemmi_record_type4(line, 0, *b"data") && line.get(4) == Some(&b'_') {
            Some(PdbInputFormatError::Cif {
                source: source.to_owned(),
            })
        } else if gemmi_record_type4(line, 0, *b"{\"da")
            && gemmi_ialpha3_id(line, 4) == gemmi_ialpha3_id(b"ta_", 0)
        {
            Some(PdbInputFormatError::Mmjson {
                source: source.to_owned(),
            })
        } else {
            None
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
struct PdbEntity {
    source_name: Vec<u8>,
    entity_kind: EntityKind,
    full_sequence: Vec<Vec<u8>>,
    dbrefs: Vec<PdbEntityDbRef>,
    subchains: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Default)]
struct PdbEntityState {
    entities: Vec<PdbEntity>,
}

#[derive(Debug, Clone, PartialEq)]
struct PdbEntityDbRef {
    db_name: Vec<u8>,
    accession_code: Vec<u8>,
    id_code: Vec<u8>,
    isoform: Vec<u8>,
    seq_begin: PdbSeqId,
    seq_end: PdbSeqId,
    db_begin: PdbSeqId,
    db_end: PdbSeqId,
    label_seq_begin: Option<i32>,
    label_seq_end: Option<i32>,
}

impl Default for PdbEntityDbRef {
    fn default() -> Self {
        // Gemmi source: third_party/gemmi/include/gemmi/metadata.hpp.
        // Gemmi❗✔️: struct DbRef {
        // Gemmi❗✔️:   std::string db_name;
        // Gemmi❗✔️:   std::string accession_code;
        // Gemmi❗✔️:   std::string id_code;
        // Gemmi❗✔️:   std::string isoform;  // pdbx_db_isoform
        // Gemmi❗✔️:   SeqId seq_begin, seq_end;
        // Gemmi❗✔️:   SeqId db_begin, db_end;
        // Gemmi❗✔️:   SeqId::OptionalNum label_seq_begin, label_seq_end;
        // Gemmi❗✔️: };
        // Gemmi❗✔️: struct SeqId {
        // Gemmi❗✔️:   using OptionalNum = OptionalInt<INT_MIN>;
        // Gemmi❗✔️:   OptionalNum num;
        // Gemmi❗✔️:   char icode = ' ';
        // Gemmi❗✔️:   SeqId() = default;
        // Gemmi❗✔️: }
        // The parser-local representation keeps source strings as bytes until
        // BioEntityRow materialization. PdbSeqId uses i32::MIN/None for the
        // source's absent number/space insertion defaults; DBREF's read_int
        // assignment of zero remains a present numeric zero.
        Self {
            db_name: Vec::new(),
            accession_code: Vec::new(),
            id_code: Vec::new(),
            isoform: Vec::new(),
            seq_begin: PdbSeqId::new(i32::MIN, None),
            seq_end: PdbSeqId::new(i32::MIN, None),
            db_begin: PdbSeqId::new(i32::MIN, None),
            db_end: PdbSeqId::new(i32::MIN, None),
            label_seq_begin: None,
            label_seq_end: None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PdbDbRefAction {
    NotDbRef,
    Continue,
    Stop,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PdbDbRefError {
    SourceIntegerOutsideDefinedRange { field_offset: usize },
    MissingCurrentDbRef,
}

#[derive(Debug, Clone, PartialEq)]
enum PdbCrystalError {
    BioStructureData(BioStructureError),
    TextFieldNotUtf8 {
        field_offset: usize,
        valid_up_to: usize,
    },
}

impl From<BioStructureError> for PdbCrystalError {
    fn from(error: BioStructureError) -> Self {
        Self::BioStructureData(error)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PdbMtrixError {
    TextFieldNotUtf8 {
        field_offset: usize,
        valid_up_to: usize,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PdbConectError {
    SourceSerialOutsideDefinedRange { field_offset: usize },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PdbModResError {
    ChainNameNotRepresentable {
        width: usize,
        non_ascii: Option<u8>,
    },
    ResidueAddress(PdbResidueKeyError),
    TextFieldNotUtf8 {
        field_offset: usize,
        valid_up_to: usize,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PdbHelixError {
    ChainNameNotRepresentable {
        field_offset: usize,
        width: usize,
        non_ascii: Option<u8>,
    },
    ResidueAddress(PdbResidueKeyError),
    IntegerOutsideSourceDefinedRange {
        field_offset: usize,
        width: usize,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PdbSheetError {
    TextFieldNotUtf8 {
        field_offset: usize,
        valid_up_to: usize,
    },
    ChainNameNotRepresentable {
        field_offset: usize,
        width: usize,
        non_ascii: Option<u8>,
    },
    ResidueAddress {
        sequence_offset: usize,
        source: PdbResidueKeyError,
    },
    IntegerOutsideSourceDefinedRange {
        field_offset: usize,
        width: usize,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PdbHeaderError {
    TextFieldNotUtf8 {
        field_offset: usize,
        valid_up_to: usize,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PdbRemarkError {
    LineLengthOutsideBuffer { line_len: usize },
    TextFieldNotUtf8 { valid_up_to: usize },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PdbRemark3Error {
    IntegerOutsideSourceDefinedRange,
    TextNotUtf8 { valid_up_to: usize },
    MissingContinuationTarget,
    TlsChainIdNotRepresentable { width: usize },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PdbRemark200Error {
    IntegerOutsideSourceDefinedRange,
    TextNotUtf8 { valid_up_to: usize },
    MissingExperimentalRows,
    UndefinedPhNumericScreening,
    Software(PdbRemark3Error),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PdbRemarkMetadataError {
    IntegerOutsideSourceDefinedRange,
    UndefinedResolutionField { remark_length: usize },
    TextBoundaryNotUtf8 { offset: usize },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PdbRemark350Error {
    IntegerOutsideSourceDefinedRange,
    SourceUndefinedAccess {
        offset: usize,
        width: usize,
        record_length: usize,
    },
    TextNotUtf8 {
        field_offset: usize,
        valid_up_to: usize,
    },
    StateInvariantLost,
}

fn remark3_int(value: &[u8]) -> Result<i32, PdbRemark3Error> {
    // Gemmi❗✔️:       ref_info.reflection_count = std::atoi(value);
    // Gemmi❗✔️:       ref_info.bin_count = std::atoi(value);
    // Gemmi❗✔️:       ref_info.rfree_set_count = atoi(value);
    // Gemmi source uses `std::atoi`/`atoi` here. Within the source-defined int
    // range, the existing PDB `read_int` primitive has the same C-locale
    // whitespace, optional sign and decimal-prefix behavior. Overflow is
    // undefined in the source routines, so do not invent a wrapped or zero
    // value for it in this canonical helper.
    read_int(value).ok_or(PdbRemark3Error::IntegerOutsideSourceDefinedRange)
}

fn remark3_text(value: &[u8]) -> Result<String, PdbRemark3Error> {
    // Gemmi❗✔️:       ref_info.cross_validation_method = std::string(value, end);
    // Gemmi❗✔️:       ref_info.rfree_selection_method = std::string(value, end);
    std::str::from_utf8(value)
        .map(str::to_owned)
        .map_err(|error| PdbRemark3Error::TextNotUtf8 {
            valid_up_to: error.valid_up_to(),
        })
}

pub(super) fn gemmi_no_sign_atoi(value: &[u8]) -> Option<(i32, usize)> {
    // Gemmi source: third_party/gemmi/include/gemmi/atox.hpp.
    // Gemmi❗✔️: inline int no_sign_atoi(const char* p, const char** endptr=nullptr) {
    // Gemmi❗✔️:   int n = 0;
    // Gemmi❗✔️:   while (is_space(*p))
    // Gemmi❗✔️:     ++p;
    // Gemmi❗✔️:   for (; is_digit(*p); ++p)
    // Gemmi❗✔️:     n = n * 10 + (*p - '0');
    // Gemmi❗✔️:   if (endptr)
    // Gemmi❗✔️:     *endptr = p;
    // Gemmi❗✔️:   return n;
    // Gemmi❗✔️: }
    // Behavior review: C-locale whitespace is skipped, no sign is consumed,
    // a leading digit prefix is accumulated, and no-conversion yields zero.
    // Checked arithmetic returns a typed source-undefined-overflow boundary.
    // Complexity review: two forward scans, O(n) time and O(1) storage, as in
    // the pinned helper; no allocation.
    let mut index = 0;
    while index < value.len() && gemmi_is_space(value[index]) {
        index += 1;
    }
    let mut number = 0_i32;
    while let Some(&digit) = value.get(index).filter(|&&byte| gemmi_is_digit(byte)) {
        number = number
            .checked_mul(10)?
            .checked_add(i32::from(digit - b'0'))?;
        index += 1;
    }
    Some((number, index))
}

fn gemmi_add_restraint_count_weight(
    refinement: &mut BioRefinementInfo,
    name: &str,
    value: &[u8],
) -> Result<(), PdbRemark3Error> {
    // Gemmi source: third_party/gemmi/src/pdb.cpp, add_restraint_count_weight.
    // Gemmi❗✔️: void add_restraint_count_weight(RefinementInfo& ref_info, const char* key, const char* value) {
    // Gemmi❗✔️:   if (*value == 'N') // NULL instead of number
    // Gemmi❗✔️:     return;
    // Gemmi❗✔️:   ref_info.restr_stats.emplace_back(key);
    // Gemmi❗✔️:   RefinementInfo::Restr& restr = ref_info.restr_stats.back();
    // Gemmi❗✔️:   const char* endptr;
    // Gemmi❗✔️:   restr.count = no_sign_atoi(value, &endptr);
    // Gemmi❗✔️:   if (const char* sep = std::strchr(endptr, ';'))
    // Gemmi❗✔️:     restr.weight = fast_atof(sep + 1, &endptr);
    // Gemmi❗✔️:   if (const char* sep = std::strchr(endptr, ';'))
    // Gemmi❗✔️:     restr.function = read_string(sep+1, 50);
    // Gemmi❗✔️: }
    // Behavior review: only an initial uppercase N suppresses insertion;
    // every other value appends a fresh row, even when repeated. The integer
    // cursor is the source no-sign digit-prefix end; each semicolon search
    // begins at the precise prior conversion end; absent fields retain -1,
    // NaN and empty-string defaults. Source-int overflow and unrepresentable
    // text remain explicit Rust boundaries rather than fabricated values.
    // Complexity review: two linear C-string separator searches, one integer
    // prefix scan and one bounded float-prefix conversion match the source's
    // linear scans; only the canonical row and function text are owned.
    if value.first() == Some(&b'N') {
        return Ok(());
    }

    let (count, count_end) =
        gemmi_no_sign_atoi(value).ok_or(PdbRemark3Error::IntegerOutsideSourceDefinedRange)?;
    let mut restraint = BioRefinementRestraint::new(name);
    restraint.count = count;

    let find_separator = |start: usize| {
        value
            .get(start..)?
            .iter()
            .take_while(|byte| **byte != 0)
            .position(|byte| *byte == b';')
            .map(|offset| start + offset)
    };

    if let Some(weight_separator) = find_separator(count_end) {
        let weight_start = weight_separator + 1;
        let (weight, consumed) =
            crate::bio_numeric::gemmi_fast_atof_with_end(&value[weight_start..]);
        restraint.weight = weight;
        if let Some(function_separator) = find_separator(weight_start + consumed) {
            let function_start = function_separator + 1;
            let function_field = value.get(function_start..).unwrap_or_default();
            let function_field = &function_field[..function_field.len().min(50)];
            restraint.function = remark3_text(read_string(function_field))?;
        }
    }

    refinement.restr_stats.push(restraint);
    Ok(())
}

fn gemmi_is_tls_item(key: &[u8]) -> bool {
    // Gemmi source: third_party/gemmi/src/pdb.cpp, is_tls_item.
    // Gemmi❗✔️: bool is_tls_item(const std::string& key) {
    // Gemmi❗✔️:   return key.size() == 3 &&
    // Gemmi❗✔️:     (key[0] == 'T' || key[0] == 'L' || key[0] == 'S') &&
    // Gemmi❗✔️:     (key[1] == '1' || key[1] == '2' || key[1] == '3') &&
    // Gemmi❗✔️:     (key[2] == '1' || key[2] == '2' || key[2] == '3');
    // Gemmi❗✔️: }
    // Behavior review: retain the exact three-byte uppercase family and
    // one-based 1..=3 row/column grammar. Complexity review: both paths inspect
    // at most three bytes with no allocation, matching the source's bounded
    // comparisons.
    key.len() == 3
        && matches!(key[0], b'T' | b'L' | b'S')
        && matches!(key[1], b'1'..=b'3')
        && matches!(key[2], b'1'..=b'3')
}

fn gemmi_split_str_multi(input: &[u8]) -> Vec<&[u8]> {
    // Gemmi source: third_party/gemmi/include/gemmi/util.hpp,
    // split_str_into_multi and split_str_multi.
    // Gemmi❗🔝: inline void split_str_into_multi(const std::string& str, const char* seps,
    // Gemmi❗🔝:                                  std::vector<std::string>& result) {
    // Gemmi❗🔝:   std::size_t start = str.find_first_not_of(seps);
    // Gemmi❗🔝:   while (start != std::string::npos) {
    // Gemmi❗🔝:     std::size_t end = str.find_first_of(seps, start);
    // Gemmi❗🔝:     result.emplace_back(str, start, end - start);
    // Gemmi❗🔝:     start = str.find_first_not_of(seps, end);
    // Gemmi❗🔝:   }
    // Gemmi❗🔝: }
    // Gemmi❗🔝: inline std::vector<std::string> split_str_multi(const std::string& str,
    // Gemmi❗🔝:                                                 const char* seps=" \t") {
    // Gemmi❗🔝:   std::vector<std::string> result;
    // Gemmi❗🔝:   split_str_into_multi(str, seps, result);
    // Gemmi❗🔝:   return result;
    // Gemmi❗🔝: }
    // Behavior review: separator bytes are exactly ASCII space and tab and
    // empty runs produce no tokens. Complexity review: one O(n) scan plus a
    // token vector; borrowed slices avoid the source's per-token owned-string
    // copies without changing token content or order.
    let mut tokens = Vec::new();
    let mut start = None;
    for (index, byte) in input.iter().copied().enumerate() {
        if matches!(byte, b' ' | b'\t') {
            if let Some(start) = start.take() {
                tokens.push(&input[start..index]);
            }
        } else if start.is_none() {
            start = Some(index);
        }
    }
    if let Some(start) = start {
        tokens.push(&input[start..]);
    }
    tokens
}

fn gemmi_tls_seq_id(value: &[u8]) -> Result<Option<PdbSeqId>, PdbRemark3Error> {
    // Gemmi source: third_party/gemmi/include/gemmi/seqid.hpp.
    // Gemmi❗✔️:   explicit SeqId(const std::string& str) {
    // Gemmi❗✔️:     char* endptr;
    // Gemmi❗✔️:     num = std::strtol(str.c_str(), &endptr, 10);
    // Gemmi❗✔️:     if (endptr == str.c_str() || (*endptr != '\0' && endptr[1] != '\0'))
    // Gemmi❗✔️:       throw std::invalid_argument("Not a seqid: " + str);
    // Gemmi❗✔️:     icode = (*endptr | 0x20);
    // Gemmi❗✔️:   }
    // Gemmi❗✔️: Behavior review: after the source fixed-field trim, accept a
    // decimal prefix plus at most one insertion byte; a no-digit or longer
    // suffix is the source invalid_argument path, caught by the caller to omit
    // that selection. The canonical PdbSeqId is i32; source-int overflow is a
    // typed boundary, not a guessed or narrowed value. Complexity review:
    // one bounded fixed-field scan and the existing read_int primitive.
    let mut numeric_end = 0;
    while numeric_end < value.len() && gemmi_is_space(value[numeric_end]) {
        numeric_end += 1;
    }
    if matches!(value.get(numeric_end), Some(b'+' | b'-')) {
        numeric_end += 1;
    }
    let digit_start = numeric_end;
    while numeric_end < value.len() && gemmi_is_digit(value[numeric_end]) {
        numeric_end += 1;
    }
    if numeric_end == digit_start || (numeric_end < value.len() && numeric_end + 1 != value.len()) {
        return Ok(None);
    }
    let sequence_number =
        read_int(&value[..numeric_end]).ok_or(PdbRemark3Error::IntegerOutsideSourceDefinedRange)?;
    let insertion = value.get(numeric_end).copied().unwrap_or(0) | 0x20;
    Ok(Some(PdbSeqId::new(
        sequence_number,
        (insertion != b' ').then_some(insertion),
    )))
}

fn remark3_fixed_field(line: &[u8], offset: usize, width: usize) -> &[u8] {
    // Gemmi source: third_party/gemmi/src/pdb.cpp, read_string.
    // Gemmi❗✔️:   // left trim
    // Gemmi❗✔️:   while (field_length != 0 && is_space(*p)) {
    // Gemmi❗✔️:     ++p;
    // Gemmi❗✔️:     --field_length;
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   // EOL/EOF ends the string
    // Gemmi❗✔️:   for (int i = 0; i < field_length; ++i)
    // Gemmi❗✔️:     if (p[i] == '\n' || p[i] == '\r' || p[i] == '\0') {
    // Gemmi❗✔️:       field_length = i;
    // Gemmi❗✔️:       break;
    // Gemmi❗✔️:     }
    // Gemmi❗✔️:   // right trim
    // Gemmi❗✔️:   while (field_length != 0 && is_space(p[field_length-1]))
    // Gemmi❗✔️:     --field_length;
    // Gemmi❗✔️:   return std::string(p, field_length);
    // The PDB reader's source line buffer is zero-filled past the copied row.
    // A Rust slice ends at that same first zero byte, so a short fixed-width
    // suffix is equivalent to the available bytes before the zero terminator.
    // The caller's source end-pointer test ensures a field start is present.
    // Complexity review: this takes one bounded slice and delegates the
    // source-shaped trim to read_string, O(width), without allocating.
    let suffix = line.get(offset..).unwrap_or_default();
    read_string(&suffix[..suffix.len().min(width)])
}

fn gemmi_tls_symmetric_component_index(row: usize, column: usize) -> usize {
    // Gemmi source: third_party/gemmi/include/gemmi/math.hpp, SMat33::unchecked_ref.
    // Gemmi❗✔️: T* ptrs[9] = {&u11, &u12, &u13, &u12, &u22, &u23, &u13, &u23, &u33};
    // Gemmi❗✔️: return *ptrs[3 * i + j];
    // Behavior review: map the source's row-major pointer table to the
    // canonical six-element u11,u22,u33,u12,u13,u23 storage; both off-diagonal
    // orientations select the same slot and later writes overwrite earlier
    // ones. Complexity review: constant-time bounded index mapping, no alloc.
    [[0, 3, 4], [3, 1, 5], [4, 5, 2]][row][column]
}

impl std::fmt::Display for PdbRemarkError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::LineLengthOutsideBuffer { line_len } => write!(
                formatter,
                "PDB remark line length {line_len} exceeds its fixed source buffer"
            ),
            Self::TextFieldNotUtf8 { valid_up_to } => write!(
                formatter,
                "PDB remark is not UTF-8 after {valid_up_to} valid bytes"
            ),
        }
    }
}

impl std::error::Error for PdbRemarkError {}

impl std::fmt::Display for PdbHeaderError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::TextFieldNotUtf8 {
                field_offset,
                valid_up_to,
            } => write!(
                formatter,
                "PDB header text at byte offset {field_offset} is not UTF-8 (valid prefix {valid_up_to} bytes)"
            ),
        }
    }
}

impl std::error::Error for PdbHeaderError {}

impl PdbEntityState {
    fn seqres_record(&mut self, line: &[u8]) -> bool {
        // Gemmi source: third_party/gemmi/src/pdb.cpp,
        // populate_structure_from_pdb_stream.
        // Gemmi✔️✔️:     } else if (is_record_type4(line, "SEQRES")) {
        // Gemmi✔️✔️:       std::string chain_name = read_string(line+10, 2);
        // Gemmi✔️✔️:       Entity& ent = impl::find_or_add(st.entities, chain_name);
        // Gemmi✔️✔️:       ent.entity_type = EntityType::Polymer;
        // Gemmi✔️✔️:       for (int i = 19; i < 68 && i < (int)len; i += 4) {
        // Gemmi✔️✔️:         std::string res_name = read_string(line+i, 3);
        // Gemmi✔️✔️:         if (!res_name.empty())
        // Gemmi✔️✔️:           ent.full_sequence.emplace_back(res_name);
        // Gemmi✔️✔️:       }
        // Gemmi source helper closure: third_party/gemmi/include/gemmi/pdb.hpp and util.hpp.
        // Gemmi✔️✔️: inline bool is_record_type4(const char* s, const char* record) {
        // Gemmi✔️✔️:   return ialpha4_id(s) == ialpha4_id(record);
        // Gemmi✔️✔️: }
        // Gemmi✔️✔️: constexpr int ialpha4_id(const char* s) {
        // Gemmi✔️✔️:   return (s[0] << 24 | s[1] << 16 | s[2] << 8 | s[3]) & ~0x20202020;
        // Gemmi✔️✔️: }
        // Gemmi source helper closure: third_party/gemmi/include/gemmi/model.hpp.
        // Gemmi✔️✔️: template<typename T>
        // Gemmi✔️✔️: auto get_id(const T& m) -> decltype(m.name) { return m.name; }
        // Gemmi✔️✔️: template<typename Vec, typename S>
        // Gemmi✔️✔️: auto find_iter_(Vec& vec, const S& name) {
        // Gemmi✔️✔️:   return std::find_if(vec.begin(), vec.end(), [&name](const auto& m) { return get_id(m) == name; });
        // Gemmi✔️✔️: }
        // Gemmi✔️✔️: template<typename T, typename S>
        // Gemmi✔️✔️: T* find_or_null(std::vector<T>& vec, const S& name) {
        // Gemmi✔️✔️:   auto it = find_iter_(vec, name);
        // Gemmi✔️✔️:   return it != vec.end() ? &*it : nullptr;
        // Gemmi✔️✔️: }
        // Gemmi✔️✔️: template<typename T, typename S>
        // Gemmi✔️✔️: T& find_or_add(std::vector<T>& vec, const S& name) {
        // Gemmi✔️✔️:   if (T* ret = find_or_null(vec, name))
        // Gemmi✔️✔️:     return *ret;
        // Gemmi✔️✔️:   vec.emplace_back(name);
        // Gemmi✔️✔️:   return vec.back();
        // Gemmi✔️✔️: }
        // Gemmi source helper closure: third_party/gemmi/include/gemmi/metadata.hpp.
        // Gemmi✔️✔️: std::string name;
        // Gemmi✔️✔️: EntityType entity_type = EntityType::Unknown;
        // Gemmi✔️✔️: std::vector<std::string> full_sequence;
        // Gemmi✔️✔️: Entity() = default;
        // Gemmi✔️✔️: explicit Entity(const std::string& name_) noexcept : name(name_) {}
        // Behavior review: record matching uses Gemmi's four-byte case-folded
        // prefix; entity identity is the exact trimmed source byte string.
        // Existing entities retain order and prior sequence, while every
        // matching record sets Polymer before appending only nonempty slots.
        // Short fields are NUL-padded by the source line buffer; slicing only
        // the available prefix has the same `read_string` result because the
        // first missing byte is the source NUL terminator.
        // Complexity review: entity selection is one linear scan, matching
        // Gemmi's find_or_add/find_or_null/find_if closure; each record reads
        // at most thirteen fixed three-byte slots and allocates only retained
        // names, with no repeated sequence copying.
        if !gemmi_record_type4(line, 0, *b"SEQR") {
            return false;
        }

        let name_end = line.len().min(12);
        let name_start = 10.min(name_end);
        let source_name = read_string(&line[name_start..name_end]).to_vec();
        let entity_index = if let Some(index) = self
            .entities
            .iter()
            .position(|entity| entity.source_name == source_name)
        {
            index
        } else {
            self.entities.push(PdbEntity {
                source_name,
                entity_kind: EntityKind::Unknown,
                full_sequence: Vec::new(),
                dbrefs: Vec::new(),
                subchains: Vec::new(),
            });
            self.entities.len() - 1
        };

        let entity = &mut self.entities[entity_index];
        entity.entity_kind = EntityKind::Polymer;
        for offset in (19..68).step_by(4) {
            if offset >= line.len() {
                break;
            }
            let field_end = line.len().min(offset + 3);
            let residue_name = read_string(&line[offset..field_end]);
            if !residue_name.is_empty() {
                entity.full_sequence.push(residue_name.to_vec());
            }
        }
        true
    }

    fn dbref_record(&mut self, line: &[u8; 122]) -> Result<PdbDbRefAction, PdbDbRefError> {
        // Gemmi source: third_party/gemmi/src/pdb.cpp,
        // populate_structure_from_pdb_stream.
        // Gemmi❗✔️:     } else if (is_record_type4(line, "DBREF")) { // DBREF or DBREF1 or DBREF2
        // Gemmi❗✔️:       std::string chain_name = read_string(line+11, 2);
        // Gemmi❗✔️:       Entity& ent = impl::find_or_add(st.entities, chain_name);
        // Gemmi❗✔️:       ent.entity_type = EntityType::Polymer;
        // Gemmi❗✔️:       if (line[5] == ' ' || line[5] == '1')
        // Gemmi❗✔️:         ent.dbrefs.emplace_back();
        // Gemmi❗✔️:       else if (ent.dbrefs.empty()) // DBREF2 without DBREF1?
        // Gemmi❗✔️:         break;
        // Gemmi❗✔️:       Entity::DbRef& dbref = ent.dbrefs.back();
        // Gemmi❗✔️:       if (line[5] == ' ' || line[5] == '1') {
        // Gemmi❗✔️:         dbref.seq_begin = read_seq_id(line+14);
        // Gemmi❗✔️:         dbref.seq_end = read_seq_id(line+20);
        // Gemmi❗✔️:         dbref.db_name = read_string(line+26, 6);
        // Gemmi❗✔️:         if (line[5] == ' ') {
        // Gemmi❗✔️:           dbref.accession_code = read_string(line+33, 8);
        // Gemmi❗✔️:           dbref.id_code = read_string(line+42, 12);
        // Gemmi❗✔️:           dbref.db_begin.num = read_int(line+55, 5);
        // Gemmi❗✔️:           dbref.db_begin.icode = line[60];
        // Gemmi❗✔️:           dbref.db_end.num = read_int(line+62, 5);
        // Gemmi❗✔️:           dbref.db_end.icode = line[67];
        // Gemmi❗✔️:         } else {  // line[5] == '1'
        // Gemmi❗✔️:           dbref.id_code = read_string(line+47, 20);
        // Gemmi❗✔️:         }
        // Gemmi❗✔️:       } else if (line[5] == '2') {
        // Gemmi❗✔️:         dbref.accession_code = read_string(line+18, 22);
        // Gemmi❗✔️:         dbref.db_begin.num = read_int(line+45, 10);
        // Gemmi❗✔️:         dbref.db_end.num = read_int(line+57, 10);
        // Gemmi❗✔️:       }
        // Gemmi❗✔️:     }
        // Gemmi source helper: third_party/gemmi/include/gemmi/pdb.hpp.
        // Gemmi❗✔️: inline bool is_record_type4(const char* s, const char* record) {
        // Gemmi❗✔️:   return ialpha4_id(s) == ialpha4_id(record);
        // Gemmi❗✔️: }
        // Gemmi❗✔️: constexpr int ialpha4_id(const char* s) {
        // Gemmi❗✔️:   return (s[0] << 24 | s[1] << 16 | s[2] << 8 | s[3]) & ~0x20202020;
        // Gemmi❗✔️: }
        // Gemmi source helper: third_party/gemmi/src/pdb.cpp::read_string.
        // Gemmi❗✔️: std::string read_string(const char* p, int field_length) {
        // Gemmi❗✔️:   // left trim
        // Gemmi❗✔️:   while (field_length != 0 && is_space(*p)) {
        // Gemmi❗✔️:     ++p;
        // Gemmi❗✔️:     --field_length;
        // Gemmi❗✔️:   }
        // Gemmi❗✔️:   // EOL/EOF ends the string
        // Gemmi❗✔️:   for (int i = 0; i < field_length; ++i)
        // Gemmi❗✔️:     if (p[i] == '\n' || p[i] == '\r' || p[i] == '\0') {
        // Gemmi❗✔️:       field_length = i;
        // Gemmi❗✔️:       break;
        // Gemmi❗✔️:     }
        // Gemmi❗✔️:   // right trim
        // Gemmi❗✔️:   while (field_length != 0 && is_space(p[field_length-1]))
        // Gemmi❗✔️:     --field_length;
        // Gemmi❗✔️:   return std::string(p, field_length);
        // Gemmi❗✔️: }
        // Gemmi source helper: third_party/gemmi/include/gemmi/model.hpp.
        // Gemmi❗✔️: template<typename T>
        // Gemmi❗✔️: auto get_id(const T& m) -> decltype(m.name) { return m.name; }
        // Gemmi❗✔️: template<typename Vec, typename S>
        // Gemmi❗✔️: auto find_iter_(Vec& vec, const S& name) {
        // Gemmi❗✔️:   return std::find_if(vec.begin(), vec.end(), [&name](const auto& m) { return get_id(m) == name; });
        // Gemmi❗✔️: }
        // Gemmi❗✔️: template<typename T, typename S>
        // Gemmi❗✔️: T* find_or_null(std::vector<T>& vec, const S& name) {
        // Gemmi❗✔️:   auto it = find_iter_(vec, name);
        // Gemmi❗✔️:   return it != vec.end() ? &*it : nullptr;
        // Gemmi❗✔️: }
        // Gemmi❗✔️: template<typename T, typename S>
        // Gemmi❗✔️: T& find_or_add(std::vector<T>& vec, const S& name) {
        // Gemmi❗✔️:   if (T* ret = find_or_null(vec, name))
        // Gemmi❗✔️:     return *ret;
        // Gemmi❗✔️:   vec.emplace_back(name);
        // Gemmi❗✔️:   return vec.back();
        // Gemmi❗✔️: }
        // Gemmi source helper: third_party/gemmi/src/pdb.cpp::read_seq_id.
        // Gemmi❗✔️: SeqId read_seq_id(const char* str) {
        // Gemmi❗✔️:   SeqId seqid;
        // Gemmi❗✔️:   if (str[4] != '\r' && str[4] != '\n')
        // Gemmi❗✔️:     seqid.icode = str[4];
        // Gemmi❗✔️:   // We support hybrid-36 extension, although it is never used in practice
        // Gemmi❗✔️:   // as 9999 residues per chain are enough.
        // Gemmi❗✔️:   if (str[0] < 'A') {
        // Gemmi❗✔️:     for (int i = 4; i != 0; --i, ++str)
        // Gemmi❗✔️:       if (!is_space(*str)) {
        // Gemmi❗✔️:         seqid.num = read_int(str, i);
        // Gemmi❗✔️:         break;
        // Gemmi❗✔️:       }
        // Gemmi❗✔️:   } else {
        // Gemmi❗✔️:     seqid.num = read_base36<4>(str) - 466560 + 10000;
        // Gemmi❗✔️:   }
        // Gemmi❗✔️:   return seqid;
        // Gemmi❗✔️: }
        // Gemmi source helper: third_party/gemmi/src/pdb.cpp::read_base36.
        // Gemmi❗✔️: template<int N> int read_base36(const char* p) {
        // Gemmi❗✔️:   char zstr[N+1] = {0};
        // Gemmi❗✔️:   std::memcpy(zstr, p, N);
        // Gemmi❗✔️:   return std::strtol(zstr, nullptr, 36);
        // Gemmi❗✔️: }
        // Gemmi source helper: third_party/gemmi/src/pdb.cpp / atox.hpp.
        // Gemmi❗✔️: int read_int(const char* p, int field_length) {
        // Gemmi❗✔️:   return string_to_int(p, false, field_length);
        // Gemmi❗✔️: }
        // Gemmi❗✔️: // no checking for overflow
        // Gemmi❗✔️: int string_to_int(const char* p, bool checked, size_t length=0) {
        // Gemmi❗✔️:   int mult = -1;
        // Gemmi❗✔️:   int n = 0;
        // Gemmi❗✔️:   size_t i = 0;
        // Gemmi❗✔️:   while ((length == 0 || i < length) && is_space(p[i]))
        // Gemmi❗✔️:     ++i;
        // Gemmi❗✔️:   if (p[i] == '-') {
        // Gemmi❗✔️:     mult = 1;
        // Gemmi❗✔️:     ++i;
        // Gemmi❗✔️:   } else if (p[i] == '+') {
        // Gemmi❗✔️:     ++i;
        // Gemmi❗✔️:   }
        // Gemmi❗✔️:   bool has_digits = false;
        // Gemmi❗✔️:   // use negative numbers because INT_MIN < -INT_MAX
        // Gemmi❗✔️:   for (; (length == 0 || i < length) && is_digit(p[i]); ++i) {
        // Gemmi❗✔️:     n = n * 10 - (p[i] - '0');
        // Gemmi❗✔️:     has_digits = true;
        // Gemmi❗✔️:   }
        // Gemmi❗✔️:   if (checked) {
        // Gemmi❗✔️:     while ((length == 0 || i < length) && is_space(p[i]))
        // Gemmi❗✔️:       ++i;
        // Gemmi❗✔️:     if (!has_digits || p[i] != '\0')
        // Gemmi❗✔️:       throw std::invalid_argument("not an integer: " +
        // Gemmi❗✔️:                           std::string(p, length ? length : i+1));
        // Gemmi❗✔️:   }
        // Gemmi❗✔️:   return mult * n;
        // Gemmi❗✔️: }
        // Gemmi source helper: third_party/gemmi/include/gemmi/atox.hpp.
        // Gemmi❗✔️: // equivalent of std::isspace for C locale (no handling of EOF)
        // Gemmi❗✔️: inline bool is_space(char c) {
        // Gemmi❗✔️:   static const std::uint8_t table[256] = { // 1 for 9-13 and 32
        // Gemmi❗✔️:     0,0,0,0,0,0,0,0, 0,1,1,1,1,1,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
        // Gemmi❗✔️:     1,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
        // Gemmi❗✔️:     0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
        // Gemmi❗✔️:     0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
        // Gemmi❗✔️:     0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
        // Gemmi❗✔️:     0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
        // Gemmi❗✔️:     0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
        // Gemmi❗✔️:     0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0
        // Gemmi❗✔️:   };
        // Gemmi❗✔️:   return table[(std::uint8_t)c] != 0;
        // Gemmi❗✔️: }
        // Gemmi❗✔️: // equivalent of std::isdigit for C locale (no handling of EOF)
        // Gemmi❗✔️: inline bool is_digit(char c) {
        // Gemmi❗✔️:   return c >= '0' && c <= '9';
        // Gemmi❗✔️: }
        // Gemmi source helper: third_party/gemmi/include/gemmi/metadata.hpp.
        // Gemmi❗✔️: struct Entity {
        // Gemmi❗✔️:   struct DbRef {
        // Gemmi❗✔️:     std::string db_name;
        // Gemmi❗✔️:     std::string accession_code;
        // Gemmi❗✔️:     std::string id_code;
        // Gemmi❗✔️:     std::string isoform;  // pdbx_db_isoform
        // Gemmi❗✔️:     SeqId seq_begin, seq_end;
        // Gemmi❗✔️:     SeqId db_begin, db_end;
        // Gemmi❗✔️:     SeqId::OptionalNum label_seq_begin, label_seq_end;
        // Gemmi❗✔️:   };
        // Gemmi❗✔️:   std::string name;
        // Gemmi❗✔️:   std::vector<std::string> subchains;
        // Gemmi❗✔️:   EntityType entity_type = EntityType::Unknown;
        // Gemmi❗✔️:   PolymerType polymer_type = PolymerType::Unknown;
        // Gemmi❗✔️:   // In case of microheterogeneity, PDB SEQRES has only the first residue name.
        // Gemmi❗✔️:   bool reflects_microhetero = false;
        // Gemmi❗✔️:   std::vector<DbRef> dbrefs;
        // Gemmi❗✔️:   /// List of SIFTS Uniprot ACs referenced by SiftsUnpResidue::acc_index
        // Gemmi❗✔️:   std::vector<std::string> sifts_unp_acc;
        // Gemmi❗✔️:   /// SEQRES or entity_poly_seq with microheterogeneity as comma-separated names
        // Gemmi❗✔️:   std::vector<std::string> full_sequence;
        // Gemmi❗✔️:
        // Gemmi❗✔️:   Entity() = default;
        // Gemmi❗✔️:   explicit Entity(const std::string& name_) noexcept : name(name_) {}
        // Gemmi❗✔️:   static std::string first_mon(const std::string& mon_list) {
        // Gemmi❗✔️:     return mon_list.substr(0, mon_list.find(','));
        // Gemmi❗✔️:   }
        // Gemmi❗✔️: };
        // The source's overflow arithmetic is undefined; the current helper
        // returns a typed boundary error for it. The source-defined no-digit
        // result remains Some(0). A Stop outcome must terminate the enclosing
        // PDB line loop after preserving the source's entity creation/promotion.
        // Behavior review: source branch order, byte columns, new-row versus
        // continuation behavior, prior values, and early termination are
        // represented explicitly; test evidence remains due in Steps276/278.
        // Complexity review: one linear source-name lookup matches Gemmi's
        // find_or_add closure; each accepted record performs bounded field
        // scans and O(1) updates, with allocations only for stored byte fields.
        if !gemmi_record_type4(line, 0, *b"DBRE") {
            return Ok(PdbDbRefAction::NotDbRef);
        }

        let source_name = read_string(&line[11..13]).to_vec();
        let entity_index = if let Some(index) = self
            .entities
            .iter()
            .position(|entity| entity.source_name == source_name)
        {
            index
        } else {
            self.entities.push(PdbEntity {
                source_name,
                entity_kind: EntityKind::Unknown,
                full_sequence: Vec::new(),
                dbrefs: Vec::new(),
                subchains: Vec::new(),
            });
            self.entities.len() - 1
        };

        let entity = &mut self.entities[entity_index];
        entity.entity_kind = EntityKind::Polymer;
        let variant = line[5];
        if variant == b' ' || variant == b'1' {
            entity.dbrefs.push(PdbEntityDbRef::default());
        } else if entity.dbrefs.is_empty() {
            return Ok(PdbDbRefAction::Stop);
        }

        let Some(dbref) = entity.dbrefs.last_mut() else {
            return Err(PdbDbRefError::MissingCurrentDbRef);
        };
        if variant == b' ' || variant == b'1' {
            let seq_begin: &[u8; 5] = line[14..19]
                .try_into()
                .expect("fixed PDB source line buffer contains DBREF sequence start");
            dbref.seq_begin = read_seq_id(seq_begin)
                .ok_or(PdbDbRefError::SourceIntegerOutsideDefinedRange { field_offset: 14 })?;
            let seq_end: &[u8; 5] = line[20..25]
                .try_into()
                .expect("fixed PDB source line buffer contains DBREF sequence end");
            dbref.seq_end = read_seq_id(seq_end)
                .ok_or(PdbDbRefError::SourceIntegerOutsideDefinedRange { field_offset: 20 })?;
            dbref.db_name = read_string(&line[26..32]).to_vec();
            if variant == b' ' {
                dbref.accession_code = read_string(&line[33..41]).to_vec();
                dbref.id_code = read_string(&line[42..54]).to_vec();
                let db_begin_num = read_int(&line[55..60])
                    .ok_or(PdbDbRefError::SourceIntegerOutsideDefinedRange { field_offset: 55 })?;
                dbref.db_begin =
                    PdbSeqId::new(db_begin_num, (line[60] != b' ').then_some(line[60]));
                let db_end_num = read_int(&line[62..67])
                    .ok_or(PdbDbRefError::SourceIntegerOutsideDefinedRange { field_offset: 62 })?;
                dbref.db_end = PdbSeqId::new(db_end_num, (line[67] != b' ').then_some(line[67]));
            } else {
                dbref.id_code = read_string(&line[47..67]).to_vec();
            }
        } else if variant == b'2' {
            dbref.accession_code = read_string(&line[18..40]).to_vec();
            let db_begin_icode = dbref.db_begin.ins_code();
            let db_begin_num = read_int(&line[45..55])
                .ok_or(PdbDbRefError::SourceIntegerOutsideDefinedRange { field_offset: 45 })?;
            dbref.db_begin = PdbSeqId::new(db_begin_num, db_begin_icode);
            let db_end_icode = dbref.db_end.ins_code();
            let db_end_num = read_int(&line[57..67])
                .ok_or(PdbDbRefError::SourceIntegerOutsideDefinedRange { field_offset: 57 })?;
            dbref.db_end = PdbSeqId::new(db_end_num, db_end_icode);
        }
        Ok(PdbDbRefAction::Continue)
    }
}

#[derive(Debug, Clone, PartialEq)]
struct PdbAtomFields {
    coordinate_text: [u8; 24],
    serial: PdbAtomSerial,
    name: AtomName,
    altloc: Option<AltLocLabel>,
    position: [f64; 3],
    occupancy: f64,
    b_iso: f64,
    anisou: [f64; 6],
    element: Element,
    isotope_mass_number: Option<u16>,
    formal_charge: i8,
}

#[derive(Debug, Clone, Copy)]
struct PdbResidueMapKey(ResidueAddress);

impl PartialEq for PdbResidueMapKey {
    fn eq(&self, other: &Self) -> bool {
        // Gemmi source: third_party/gemmi/include/gemmi/seqid.hpp
        // Gemmi❗✔️:   bool operator==(const SeqId& o) const {
        // Gemmi❗✔️:     return num == o.num && ((icode ^ o.icode) & ~0x20) == 0;
        // Gemmi❗✔️:   }
        // Gemmi❗✔️:   bool matches(const ResidueId& o) const {
        // Gemmi❗✔️:     return seqid == o.seqid && segment == o.segment && name == o.name;
        // Gemmi❗✔️:   }
        // Behavior review: delegate all key equality to the canonical source-shaped address.
        // Complexity review: sequence and insertion are scalar; bounded segment/name checks are O(1).
        self.0.matches(&other.0)
    }
}

impl Eq for PdbResidueMapKey {}

impl Hash for PdbResidueMapKey {
    fn hash<H: Hasher>(&self, state: &mut H) {
        // Gemmi source: third_party/gemmi/include/gemmi/seqid.hpp
        // Gemmi❗✔️:   size_t seqid_hash = (*r.seqid.num << 7) + (r.seqid.icode | 0x20);
        // Gemmi❗✔️:   return seqid_hash ^ hash<string>()(r.segment) ^ hash<string>()(r.name);
        // Gemmi's signed left shift is undefined for negative/sentinel and some
        // large sequence numbers. Hash bits are not observed by the reader;
        // preserve the exact source equality classes with a safe Rust hash.
        // Behavior review: the insertion byte is canonicalized by the one bit
        // ignored by SeqId equality; sequence, segment, and name remain exact.
        // Complexity review: hashing touches only scalar fields and source-bounded
        // segment/name text; HashMap lookup retains expected O(1) complexity.
        self.0.sequence_number().hash(state);
        self.0.insertion_code().map(|code| code & !0x20).hash(state);
        self.0.segment().hash(state);
        self.0.name().hash(state);
    }
}

#[derive(Debug, Clone, PartialEq)]
struct PdbResidueGroup {
    row_id: BioResidueId,
    address: ResidueAddress,
    het_flag: u8,
    entity_kind: EntityKind,
    subchain_id: String,
    atoms: Vec<PdbAtomFields>,
}

#[derive(Debug, Clone, PartialEq)]
struct PdbChainGroup {
    row_id: BioChainId,
    source_id: PdbChainId,
    residues: Vec<PdbResidueGroup>,
    residue_indices: HashMap<PdbResidueMapKey, usize>,
    current_residue: Option<usize>,
}

#[derive(Debug, Clone, PartialEq)]
struct PdbModelGroup {
    row_id: BioModelId,
    source_model_number: Option<i32>,
    chains: Vec<PdbChainGroup>,
    current_chain: Option<usize>,
    after_ter: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct PdbSsbondResolutionError {
    partner: u8,
    atom_name: String,
}

impl std::fmt::Display for PdbSsbondResolutionError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(formatter, "No such atom: {}", self.atom_name)
    }
}

impl std::error::Error for PdbSsbondResolutionError {}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PdbLinkNameCounter {
    Covale,
    MetalC,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PdbLinkRecordError {
    ChainNameNotAscii {
        partner: u8,
        byte: u8,
    },
    AtomNameNotAscii {
        partner: u8,
        byte: u8,
    },
    LinkIdNotUtf8 {
        valid_up_to: usize,
    },
    ResidueAddress {
        partner: u8,
        source: PdbResidueKeyError,
    },
    NameCounterOutsideSourceDefinedRange {
        counter: PdbLinkNameCounter,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PdbCisPepRecordError {
    RequiredFieldsPastRecord {
        record_len: usize,
        required_len: usize,
    },
    NoModelsAfterSourceDefault,
    SingleModelNumberUnavailable,
    ChainNameNotAscii {
        partner: u8,
        byte: u8,
    },
    ResidueAddress {
        partner: u8,
        source: PdbResidueKeyError,
    },
    ModelNumberOutsideSourceDefinedRange,
}

#[derive(Debug, Default)]
struct PdbLinkNameCounters {
    covale: i32,
    metalc: i32,
}

fn decode_pdb_link_address(
    record: &[u8],
    partner: u8,
    offset: usize,
) -> Result<(AtomAddress, [u8; 4]), PdbLinkRecordError> {
    // Gemmi source: third_party/gemmi/src/pdb.cpp::process_conn LINK branch.
    // Gemmi❗✔️:       for (int i : {0, 1}) {
    // Gemmi❗✔️:         const char* t = record.c_str() + 30 * i;
    // Gemmi❗✔️:         AtomAddress& ad = (i == 0 ? c.partner1 : c.partner2);
    // Gemmi❗✔️:         ad.chain_name = read_string(t + 20, 2);
    // Gemmi❗✔️:         ad.res_id = read_res_id(t + 22, t + 17);
    // Gemmi❗✔️:         ad.atom_name = read_string(t + 12, 4);
    // Gemmi❗✔️:         ad.altloc = read_altloc(t[16]);
    // Gemmi❗✔️:       }
    // Behavior review: `offset` is exactly 0 or 30 from the source's ordered
    // two-partner loop. Fixed slices retain the source's raw four-column atom
    // spelling for element inference while AtomAddress receives read_string's
    // trimmed logical name. Residue segment remains the empty value from
    // read_res_id, chain and residue encounter order are unchanged, and a
    // blank altloc maps to source NUL. Non-ASCII chain/atom bytes and a
    // nonrepresentable residue are returned as typed model-boundary errors.
    // Complexity review: constant-width slices and conversions per partner;
    // only the logical atom name allocates, as the source strings do.
    let raw_atom_name: [u8; 4] = record[offset + 12..offset + 16]
        .try_into()
        .expect("the LINK minimum-length gate covers both atom-name fields");
    AtomName::from_ascii(&raw_atom_name).ok_or_else(|| PdbLinkRecordError::AtomNameNotAscii {
        partner,
        byte: raw_atom_name
            .iter()
            .copied()
            .find(|byte| !byte.is_ascii())
            .expect("AtomName::from_ascii rejects only non-ASCII bytes"),
    })?;

    let chain_field = &record[offset + 20..offset + 22];
    let chain_name = read_string(chain_field);
    let chain_name = PdbChainId::from_ascii(chain_name).ok_or_else(|| {
        PdbLinkRecordError::ChainNameNotAscii {
            partner,
            byte: chain_name
                .iter()
                .copied()
                .find(|byte| !byte.is_ascii())
                .expect("PdbChainId::from_ascii rejects only non-ASCII bytes"),
        }
    })?;
    let sequence_field: &[u8; 5] = record[offset + 22..offset + 27]
        .try_into()
        .expect("the LINK minimum-length gate covers both sequence fields");
    let residue_name_field: &[u8; 3] = record[offset + 17..offset + 20]
        .try_into()
        .expect("the LINK minimum-length gate covers both residue-name fields");
    let residue = read_res_id(sequence_field, residue_name_field)
        .map_err(|source| PdbLinkRecordError::ResidueAddress { partner, source })?;
    let logical_atom_name = read_string(&raw_atom_name);
    let logical_atom_name = std::str::from_utf8(logical_atom_name)
        .expect("the raw AtomName ASCII invariant also covers its trimmed spelling")
        .to_owned();
    let altloc = read_altloc(record[offset + 16]).map_or(0, |label| label.value());

    Ok((
        AtomAddress::new(chain_name, residue, logical_atom_name, Some(altloc)),
        raw_atom_name,
    ))
}

fn find_pdb_link_atom<'a>(
    model: &'a PdbModelGroup,
    address: &AtomAddress,
) -> Option<&'a PdbAtomFields> {
    // Gemmi source helper: third_party/gemmi/include/gemmi/model.hpp::Model::find_cra.
    // Gemmi❗✔️:   CRA find_cra(const AtomAddress& address, bool ignore_segment=false) {
    // Gemmi❗✔️:     for (Chain& chain : chains)
    // Gemmi❗✔️:       if (chain.name == address.chain_name) {
    // Gemmi❗✔️:         for (Residue& res : chain.residues)
    // Gemmi❗✔️:           if (address.res_id.matches_noseg(res) &&
    // Gemmi❗✔️:               (ignore_segment || address.res_id.segment == res.segment)) {
    // Gemmi❗✔️:             Atom *at = nullptr;
    // Gemmi❗✔️:             if (!address.atom_name.empty())
    // Gemmi❗✔️:               at = res.find_atom(address.atom_name, address.altloc);
    // Gemmi❗✔️:             return {&chain, &res, at};
    // Gemmi❗✔️:           }
    // Gemmi❗✔️:       }
    // Gemmi❗✔️:     return {nullptr, nullptr, nullptr};
    // Gemmi❗✔️:   }
    // Gemmi source helper: third_party/gemmi/include/gemmi/model.hpp::Residue::find_atom.
    // Gemmi❗✔️:   Atom* find_atom(const std::string& atom_name, char altloc, El el=El::X,
    // Gemmi❗✔️:                   bool strict_altloc=true) {
    // Gemmi❗✔️:     if (!strict_altloc && altloc == '\0')
    // Gemmi❗✔️:       altloc = '*';
    // Gemmi❗✔️:     for (Atom& a : atoms)
    // Gemmi❗✔️:       if (a.name == atom_name && a.altloc_matches(altloc) &&
    // Gemmi❗✔️:           (el == El::X || a.element == el))
    // Gemmi❗✔️:         return &a;
    // Gemmi❗✔️:     return nullptr;
    // Gemmi❗✔️:   }
    // Behavior review: select the first source-order matching chain/residue,
    // with exact segment equality, then the first logical atom-name match
    // whose stored blank/requested altloc passes Gemmi's default strict
    // matching rule. A matching residue with no atom terminates lookup, just
    // as find_cra returns a residue with a null atom rather than continuing.
    // Complexity review: one ordered hierarchy scan and one atom scan in the
    // selected residue; no temporary indexes or hierarchy copies.
    for chain in &model.chains {
        if chain.source_id != address.chain_name() {
            continue;
        }
        for residue in &chain.residues {
            let requested_residue = address.residue();
            if !requested_residue.matches_without_segment(&residue.address)
                || requested_residue.segment() != residue.address.segment()
            {
                continue;
            }
            if address.logical_atom_name().is_empty() {
                return None;
            }
            let request = if address.altloc() == b'*' {
                AltLocRequest::Any
            } else {
                AltLocRequest::Exact(
                    (address.altloc() != 0).then(|| AltLocLabel::new(address.altloc())),
                )
            };
            return residue.atoms.iter().find(|atom| {
                read_string(atom.name.as_bytes()).as_ref() == address.logical_atom_name().as_bytes()
                    && cosmolkit_bio::altloc_matches(atom.altloc, request)
            });
        }
    }
    None
}

fn pdb_link_partner_element(
    model: &PdbModelGroup,
    address: &AtomAddress,
    raw_atom_name: &[u8; 4],
) -> Element {
    // Gemmi source: third_party/gemmi/src/pdb.cpp::process_conn LINK branch.
    // Gemmi❗✔️:       auto get_elem = [&](const char* name, const AtomAddress& ad) {
    // Gemmi❗✔️:         if (element_from_padded_name_is_ambiguous(name)) {
    // Gemmi❗✔️:           const_CRA cra = st.first_model().find_cra(ad);
    // Gemmi❗✔️:           if (cra.atom)
    // Gemmi❗✔️:             return cra.atom->element.elem;
    // Gemmi❗✔️:         }
    // Gemmi❗✔️:         return infer_element_from_padded_name(name);
    // Gemmi❗✔️:       };
    // Behavior review: only source-ambiguous four-column names consult the
    // first-model address lookup. A found atom's explicit parsed element wins;
    // otherwise the existing canonical padded-name inference is used. The D
    // isotope representation does not affect the following metal predicate.
    // Complexity review: fixed-width ambiguity checks are constant time; an
    // ambiguous name incurs the single ordered model/address scan above.
    let bytes = raw_atom_name;
    let ambiguous = bytes[0] != b' '
        && bytes[3] != b' '
        && !gemmi_is_digit(bytes[0])
        && !gemmi_is_digit(bytes[1]);
    if ambiguous && let Some(atom) = find_pdb_link_atom(model, address) {
        return atom.element;
    }
    infer_element_from_padded_name(raw_atom_name).0
}

fn gemmi_default_element_is_metal(element: Element) -> bool {
    // Gemmi source: third_party/gemmi/include/gemmi/elem.hpp::is_metal_value.
    // Gemmi❗✔️: inline bool& is_metal_value(El el) {
    // Gemmi❗✔️:   static bool table[] = {
    // Gemmi❗✔️:     // X     H     He
    // Gemmi❗✔️:     false, false, false,
    // Gemmi❗✔️:     // Li  Be     B      C      N      O      F     Ne
    // Gemmi❗✔️:     true, true, false, false, false, false, false, false,
    // Gemmi❗✔️:     // Na  Mg    Al     Si     P      S      Cl     Ar
    // Gemmi❗✔️:     true, true, true, false, false, false, false, false,
    // Gemmi❗✔️:     // K   Ca    Sc    Ti    V     Cr    Mn    Fe    Co    Ni    Cu    Zn
    // Gemmi❗✔️:     true, true, true, true, true, true, true, true, true, true, true, true,
    // Gemmi❗✔️:     // Ga  Ge    As     Se     Br     Kr
    // Gemmi❗✔️:     true, true, false, false, false, false,
    // Gemmi❗✔️:     // Rb  Sr    Y     Zr    Nb    Mo    Tc    Ru    Rh    Pd    Ag    Cd
    // Gemmi❗✔️:     true, true, true, true, true, true, true, true, true, true, true, true,
    // Gemmi❗✔️:     // In  Sn    Sb    Te      I     Xe
    // Gemmi❗✔️:     true, true, true, false, false, false,
    // Gemmi❗✔️:     // Cs  Ba    La    Ce    Pr    Nd    Pm    Sm    Eu    Gd    Tb    Dy
    // Gemmi❗✔️:     true, true, true, true, true, true, true, true, true, true, true, true,
    // Gemmi❗✔️:     // Ho  Er    Tm    Yb    Lu    Hf    Ta    W     Re    Os    Ir    Pt
    // Gemmi❗✔️:     true, true, true, true, true, true, true, true, true, true, true, true,
    // Gemmi❗✔️:     // Au  Hg    Tl    Pb    Bi    Po    At     Rn
    // Gemmi❗✔️:     true, true, true, true, true, false, false, false,
    // Gemmi❗✔️:     // Fr  Ra    Ac    Th    Pa    U     Np    Pu    Am    Cm    Bk    Cf
    // Gemmi❗✔️:     true, true, true, true, true, true, true, true, true, true, true, true,
    // Gemmi❗✔️:     // Es  Fm    Md    No    Lr    Rf    Db    Sg    Bh    Hs    Mt    Ds
    // Gemmi❗✔️:     true, true, true, true, true, true, true, true, true, true, true, true,
    // Gemmi❗✔️:     // Rg  Cn    Nh    Fl    Mc    Lv    Ts     Og
    // Gemmi❗✔️:     true, true, true, true, true, true, false, false,
    // Gemmi❗✔️:     // D    END
    // Gemmi❗✔️:     false, false
    // Gemmi❗✔️:   };
    // Gemmi❗✔️:   return table[static_cast<int>(el)];
    // Gemmi❗✔️: }
    // Gemmi source: third_party/gemmi/include/gemmi/elem.hpp.
    // Gemmi❗✔️: inline bool is_metal(El el) { return is_metal_value(el); }
    // Behavior review: these atomic-number ranges encode the pinned table's
    // initialized default values. Gemmi also exposes set_is_metal(), which
    // mutates its process-global table; COSMolKit does not model that separate
    // global override, so this helper does not claim parity after such an
    // external mutation. Under the parser's pinned default element state, the
    // ranges include the source's nonintuitive Ge/Sb and superheavy entries.
    // Complexity review: a fixed series of integer-range comparisons is
    // O(1), replacing the source's O(1) array lookup with no allocation.
    matches!(
        element.atomic_number(),
        3..=4 | 11..=13 | 19..=32 | 37..=51 | 55..=83 | 87..=116
    )
}

fn compare_link_symops(record: &[u8], reported_sym: &mut [i16; 4]) -> BioAsu {
    // Gemmi source: third_party/gemmi/src/pdb.cpp::compare_link_symops.
    // Gemmi❗✔️: Asu compare_link_symops(const std::string& record, short* reported_sym) {
    // Gemmi❗✔️:   if (record.size() < 72)
    // Gemmi❗✔️:     return Asu::Any;  // it could be interpreted as Same
    // Gemmi❗✔️:   std::string s1 = read_string(&record[59], 6);
    // Gemmi❗✔️:   std::string s2 = read_string(&record[66], 6);
    // Gemmi❗✔️:   if (s1 == s2)
    // Gemmi❗✔️:     return Asu::Same;
    // Gemmi❗✔️:   size_t len1 = s1.length();
    // Gemmi❗✔️:   size_t len2 = s2.length();
    // Gemmi❗✔️:   if (len1 >= 4 && len1 < 6 && len2 >= 4 && len2 < 6) {
    // Gemmi❗✔️:     // for 5 digits, we assume here that two digits are for sym_idx
    // Gemmi❗✔️:     if (s1[0] == '1' && len1 == 4)
    // Gemmi❗✔️:       reported_sym[0] = (short) read_int(s2.c_str(), len2 - 3);
    // Gemmi❗✔️:     else
    // Gemmi❗✔️:       reported_sym[0] = 99;
    // Gemmi❗✔️:     for (size_t i = 1; i <= 3; ++i)
    // Gemmi❗✔️:       reported_sym[i] = s2[len2 - 4 + i] - s1[len1 - 4 + i];
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   return Asu::Different;
    // Gemmi❗✔️: }
    // Behavior review: short records leave the caller's zero-initialized
    // vector unchanged and return Any; equal trimmed fields return Same before
    // decoding; only unequal 4/5-byte fields write deltas, with source's
    // special `1xxx`/four-byte rule and fixed-width read_int prefix conversion.
    // Complexity review: two six-byte trims and at most four scalar writes;
    // constant time and auxiliary space.
    if record.len() < 72 {
        return BioAsu::Any;
    }
    let symmetry1 = read_string(&record[59..65]);
    let symmetry2 = read_string(&record[66..72]);
    if symmetry1 == symmetry2 {
        return BioAsu::Same;
    }
    let len1 = symmetry1.len();
    let len2 = symmetry2.len();
    if (4..6).contains(&len1) && (4..6).contains(&len2) {
        reported_sym[0] = if symmetry1[0] == b'1' && len1 == 4 {
            let sym_index_width = len2 - 3;
            read_int(&symmetry2[..sym_index_width])
                .expect("one or two source decimal bytes cannot overflow i32") as i16
        } else {
            99
        };
        for index in 1..=3 {
            reported_sym[index] =
                i16::from(symmetry2[len2 - 4 + index]) - i16::from(symmetry1[len1 - 4 + index]);
        }
    }
    BioAsu::Different
}

fn decode_pdb_link_record(
    record: &[u8],
    first_model: &PdbModelGroup,
    counters: &mut PdbLinkNameCounters,
) -> Result<Option<BioConnection>, PdbLinkRecordError> {
    // Gemmi source: third_party/gemmi/src/pdb.cpp::process_conn LINK branch.
    // Gemmi❗✔️:     } else if (record[0] == 'L' || record[0] == 'l') { // LINK
    // Gemmi❗✔️:       if (record.length() < 57)
    // Gemmi❗✔️:         continue;
    // Gemmi❗✔️:       Connection c;
    // Gemmi❗✔️:       auto get_elem = [&](const char* name, const AtomAddress& ad) {
    // Gemmi❗✔️:         if (element_from_padded_name_is_ambiguous(name)) {
    // Gemmi❗✔️:           const_CRA cra = st.first_model().find_cra(ad);
    // Gemmi❗✔️:           if (cra.atom)
    // Gemmi❗✔️:             return cra.atom->element.elem;
    // Gemmi❗✔️:         }
    // Gemmi❗✔️:         return infer_element_from_padded_name(name);
    // Gemmi❗✔️:       };
    // Gemmi❗✔️:       if (is_metal(get_elem(&record[12], c.partner1)) ||
    // Gemmi❗✔️:           is_metal(get_elem(&record[42], c.partner2))) {
    // Gemmi❗✔️:         c.name = "metalc" + std::to_string(++metalc_count);
    // Gemmi❗✔️:         c.type = Connection::MetalC;
    // Gemmi❗✔️:       } else {
    // Gemmi❗✔️:         c.name = "covale" + std::to_string(++covale_count);
    // Gemmi❗✔️:         c.type = Connection::Covale;
    // Gemmi❗✔️:       }
    // Gemmi❗✔️:       c.asu = compare_link_symops(record, c.reported_sym);
    // Gemmi❗✔️:       if (record.length() > 73) {
    // Gemmi❗✔️:         if (record[4] == 'R')
    // Gemmi❗✔️:           c.link_id = read_string(&record[72], 8);
    // Gemmi❗✔️:         else
    // Gemmi❗✔️:           c.reported_distance = read_double(&record[73], 5);
    // Gemmi❗✔️:       }
    // Gemmi❗✔️:       st.connections.emplace_back(c);
    // Behavior review: callers dispatch LINK/LINKR records here in source
    // order and supply counters initialized by process_conn. Records shorter
    // than 57 bytes are skipped before field reads; accepted rows construct
    // both addresses, source-classify by the first-model element lookup or
    // padded-name inference, increment exactly one type-specific name counter,
    // preserve compare_link_symops' ASU/symmetry state, and only inspect the
    // optional link-id/distance region when length exceeds 73. A short optional
    // suffix terminates at the Rust slice boundary, matching the source's
    // NUL-terminated record string for defined fields. The mutable Gemmi metal
    // override remains outside the modeled default-state contract documented
    // by `gemmi_default_element_is_metal`.
    // Complexity review: fixed-width parsing is O(1); at most two ambiguous
    // partner lookups scan the first model in source order. Connection and
    // names allocate the same retained values as Gemmi; no complete model copy.
    if record.len() < 57 {
        return Ok(None);
    }

    let (partner1, atom_name1) = decode_pdb_link_address(record, 1, 0)?;
    let (partner2, atom_name2) = decode_pdb_link_address(record, 2, 30)?;
    let first_is_metal = gemmi_default_element_is_metal(pdb_link_partner_element(
        first_model,
        &partner1,
        &atom_name1,
    ));
    let is_metal = first_is_metal
        || gemmi_default_element_is_metal(pdb_link_partner_element(
            first_model,
            &partner2,
            &atom_name2,
        ));

    let (name, kind) = if is_metal {
        counters.metalc = counters.metalc.checked_add(1).ok_or(
            PdbLinkRecordError::NameCounterOutsideSourceDefinedRange {
                counter: PdbLinkNameCounter::MetalC,
            },
        )?;
        (
            format!("metalc{}", counters.metalc),
            BioConnectionKind::MetalC,
        )
    } else {
        counters.covale = counters.covale.checked_add(1).ok_or(
            PdbLinkRecordError::NameCounterOutsideSourceDefinedRange {
                counter: PdbLinkNameCounter::Covale,
            },
        )?;
        (
            format!("covale{}", counters.covale),
            BioConnectionKind::Covale,
        )
    };

    let mut reported_sym = [0; 4];
    let asu = compare_link_symops(record, &mut reported_sym);
    let mut connection = BioConnection {
        name,
        kind,
        asu,
        partner1,
        partner2,
        reported_sym,
        ..BioConnection::default()
    };
    if record.len() > 73 {
        if record[4] == b'R' {
            let link_id_end = record.len().min(80);
            let link_id = read_string(&record[72..link_id_end]);
            connection.link_id = std::str::from_utf8(link_id)
                .map_err(|error| PdbLinkRecordError::LinkIdNotUtf8 {
                    valid_up_to: error.valid_up_to(),
                })?
                .to_owned();
        } else {
            let distance_end = record.len().min(78);
            connection.reported_distance = read_double(&record[73..distance_end]);
        }
    }
    Ok(Some(connection))
}

fn decode_pdb_cispep_record(
    record: &[u8],
    models: &[PdbModelGroup],
) -> Result<Option<BioCisPep>, PdbCisPepRecordError> {
    // Gemmi source: third_party/gemmi/src/pdb.cpp::process_conn.
    // Gemmi❗✔️:     } else if (record[0] == 'C' || record[0] == 'c') { // CISPEP
    // Gemmi❗✔️:       if (record.length() < 22)
    // Gemmi❗✔️:         continue;
    // Gemmi❗✔️:       const char* r = record.c_str();
    // Gemmi❗✔️:       CisPep cispep;
    // Gemmi❗✔️:       cispep.partner_c.chain_name = read_string(r + 14, 2);
    // Gemmi❗✔️:       cispep.partner_c.res_id = read_res_id(r + 17, r + 11);
    // Gemmi❗✔️:       cispep.partner_n.chain_name = read_string(r + 28, 2);
    // Gemmi❗✔️:       cispep.partner_n.res_id = read_res_id(r + 31, r + 25);
    // Gemmi❗✔️:       // In files with a single model in the PDB CISPEP modNum is 0,
    // Gemmi❗✔️:       // but _struct_mon_prot_cis.pdbx_PDB_model_num is 1.
    // Gemmi❗✔️:       cispep.model_num = st.models.size() == 1 ? st.models[0].num : read_int(r + 43, 3);
    // Gemmi❗✔️:       cispep.reported_angle = read_double(r + 53, 6);
    // Gemmi❗✔️:       st.cispeps.push_back(cispep);
    // Gemmi❗✔️:     }
    // Behavior review: the source dispatch calls this branch only for a
    // captured CISPEP row. Preserve the explicit <22 skip, the two independent
    // residue-only AtomAddresses (empty atom names/NUL altloc), and the
    // single-model rule that uses that model's source number instead of the
    // record modNum field. PDB never writes BioCisPep.only_altloc, so its
    // canonical NUL default remains. Gemmi's branch reads through byte 58
    // without checking lengths 22..58; those inputs cross the std::string
    // storage boundary and are not source-defined. Rust reports that exact
    // safety boundary rather than reading beyond the record or inventing
    // padding. Chain/residue representation failures and an impossible
    // missing model number remain structured errors; read_double preserves
    // Gemmi's initialized-zero/ignored-status semantics for a present field.
    // Complexity review: all source fields have fixed width; two canonical
    // address values and one result are built in O(1) time/storage, without
    // copying or scanning the structure's models beyond its length/first row.
    if record.len() < 22 {
        return Ok(None);
    }
    if record.len() < 59 {
        return Err(PdbCisPepRecordError::RequiredFieldsPastRecord {
            record_len: record.len(),
            required_len: 59,
        });
    }
    if models.is_empty() {
        return Err(PdbCisPepRecordError::NoModelsAfterSourceDefault);
    }

    let decode_partner = |partner: u8,
                          chain_offset: usize,
                          sequence_offset: usize,
                          residue_name_offset: usize|
     -> Result<AtomAddress, PdbCisPepRecordError> {
        let chain_bytes = read_string(&record[chain_offset..chain_offset + 2]);
        let chain_name = PdbChainId::from_ascii(chain_bytes).ok_or_else(|| {
            PdbCisPepRecordError::ChainNameNotAscii {
                partner,
                byte: chain_bytes
                    .iter()
                    .copied()
                    .find(|byte| !byte.is_ascii())
                    .expect("PdbChainId::from_ascii rejects only non-ASCII bytes"),
            }
        })?;
        let sequence: &[u8; 5] = record[sequence_offset..sequence_offset + 5]
            .try_into()
            .expect("the fixed CISPEP record-length gate covers the residue sequence");
        let residue_name: &[u8; 3] = record[residue_name_offset..residue_name_offset + 3]
            .try_into()
            .expect("the fixed CISPEP record-length gate covers the residue name");
        let residue = read_res_id(sequence, residue_name)
            .map_err(|source| PdbCisPepRecordError::ResidueAddress { partner, source })?;
        Ok(AtomAddress::new(chain_name, residue, "", None))
    };

    let partner_c = decode_partner(1, 14, 17, 11)?;
    let partner_n = decode_partner(2, 28, 31, 25)?;
    let model_num = if models.len() == 1 {
        models[0]
            .source_model_number
            .ok_or(PdbCisPepRecordError::SingleModelNumberUnavailable)?
    } else {
        read_int(&record[43..46])
            .ok_or(PdbCisPepRecordError::ModelNumberOutsideSourceDefinedRange)?
    };
    let reported_angle = read_double(&record[53..59]);

    let mut cispep = BioCisPep::default();
    cispep.partner_c = partner_c;
    cispep.partner_n = partner_n;
    cispep.model_num = model_num;
    cispep.reported_angle = reported_angle;
    Ok(Some(cispep))
}

fn ssbond_atom_group<'a>(
    address: &'a AtomAddress,
    residue: &'a PdbResidueGroup,
) -> impl Iterator<Item = &'a PdbAtomFields> + 'a {
    // Gemmi source: third_party/gemmi/include/gemmi/model.hpp::Residue::get.
    // Gemmi❗✔️: AtomGroup get(const std::string& atom_name) {
    // Gemmi❗✔️:   for (Atom& atom : atoms)
    // Gemmi❗✔️:     if (atom.name == atom_name)
    // Gemmi❗✔️:       return AtomGroup(&atom, atoms.data() + atoms.size());
    // Gemmi❗✔️:   fail("No such atom: " + atom_name);
    // Gemmi❗✔️: }
    // Gemmi source: third_party/gemmi/include/gemmi/iterator.hpp::ItemGroup.
    // Gemmi❗✔️: ItemGroup(Item* start, const Item* end)
    // Gemmi❗✔️:     : size_(int(end - start)), extent_(int(end - start)), start_(start) {
    // Gemmi❗✔️:   for (const Item* i = start + 1; i != end; ++i)
    // Gemmi❗✔️:     if (i->group_key() != start->group_key())
    // Gemmi❗✔️:       --size_;
    // Gemmi❗✔️: }
    // Gemmi❗✔️: iterator& operator++() {
    // Gemmi❗✔️:   const Item* prev = ptr++;
    // Gemmi❗✔️:   while (ptr != end && ptr->group_key() != prev->group_key())
    // Gemmi❗✔️:     ++ptr;
    // Gemmi❗✔️:   return *this;
    // Gemmi❗✔️: }
    // Behavior review: locate the first requested logical name, then yield all
    // later atoms with that same name in encounter order, including nonadjacent
    // alternatives. An absent group is distinguished by the caller and becomes
    // the source-shaped `No such atom` error; no altloc filtering occurs here.
    // Complexity review: one lazy linear scan per iteration, with no copied
    // group or allocation; repeated partner-two scans match the source loop.
    residue.atoms.iter().filter(move |atom| {
        read_string(atom.name.as_bytes()) == address.logical_atom_name().as_bytes()
    })
}

fn complete_ssbond_atom<'a>(
    address: AtomAddress,
    model: &'a PdbModelGroup,
) -> (AtomAddress, Option<&'a PdbResidueGroup>) {
    // Gemmi source: third_party/gemmi/src/pdb.cpp::complete_ssbond_atom.
    // Gemmi❗✔️: const Residue* complete_ssbond_atom(AtomAddress& ad, const Model& mdl) {
    // Gemmi❗✔️:   ad.atom_name = "SG";
    // Gemmi❗✔️:   const_CRA cra = mdl.find_cra(ad);
    // Gemmi❗✔️:   if (cra.residue && (!cra.atom || cra.atom->element != El::S))
    // Gemmi❗✔️:     if (const Atom* a = cra.residue->find_by_element(El::S)) {
    // Gemmi❗✔️:       ad.atom_name = a->name;
    // Gemmi❗✔️:       ad.altloc = a->altloc;
    // Gemmi❗✔️:     }
    // Gemmi❗✔️:   return cra.residue;
    // Gemmi❗✔️: }
    // Gemmi source helper: third_party/gemmi/include/gemmi/model.hpp::Model::find_cra.
    // Gemmi❗✔️: CRA find_cra(const AtomAddress& address, bool ignore_segment=false) {
    // Gemmi❗✔️:   for (Chain& chain : chains)
    // Gemmi❗✔️:     if (chain.name == address.chain_name) {
    // Gemmi❗✔️:       for (Residue& res : chain.residues)
    // Gemmi❗✔️:         if (address.res_id.matches_noseg(res) &&
    // Gemmi❗✔️:             (ignore_segment || address.res_id.segment == res.segment)) {
    // Gemmi❗✔️:           Atom *at = nullptr;
    // Gemmi❗✔️:           if (!address.atom_name.empty())
    // Gemmi❗✔️:             at = res.find_atom(address.atom_name, address.altloc);
    // Gemmi❗✔️:           return {&chain, &res, at};
    // Gemmi❗✔️:         }
    // Gemmi❗✔️:     }
    // Gemmi❗✔️:   return {nullptr, nullptr, nullptr};
    // Gemmi❗✔️: }
    // Gemmi source helper: third_party/gemmi/include/gemmi/model.hpp::Residue::find_atom.
    // Gemmi❗✔️: Atom* find_atom(const std::string& atom_name, char altloc, El el=El::X,
    // Gemmi❗✔️:                 bool strict_altloc=true) {
    // Gemmi❗✔️:   if (!strict_altloc && altloc == '\0')
    // Gemmi❗✔️:     altloc = '*';
    // Gemmi❗✔️:   for (Atom& a : atoms)
    // Gemmi❗✔️:     if (a.name == atom_name && a.altloc_matches(altloc) && (el == El::X || a.element == el))
    // Gemmi❗✔️:       return &a;
    // Gemmi❗✔️:   return nullptr;
    // Gemmi❗✔️: }
    // Gemmi source helper: third_party/gemmi/include/gemmi/model.hpp::Residue::find_by_element.
    // Gemmi❗✔️: const Atom* find_by_element(El el) const {
    // Gemmi❗✔️:   for (const Atom& a : atoms)
    // Gemmi❗✔️:     if (a.element == el)
    // Gemmi❗✔️:       return &a;
    // Gemmi❗✔️:   return nullptr;
    // Gemmi❗✔️: }
    // Gemmi source helper: third_party/gemmi/include/gemmi/model.hpp::Atom::altloc_matches.
    // Gemmi❗✔️: bool altloc_matches(char request) const {
    // Gemmi❗✔️:   return request == '*' || altloc == '\0' || altloc == request;
    // Gemmi❗✔️: }
    // Behavior review: reset the logical name to SG, select the first source-
    // order chain/residue matching find_cra (including its exact segment gate),
    // then use find_atom's altloc rule. A missing or non-sulfur SG falls back
    // to the first sulfur regardless of altloc and copies that atom's name and
    // altloc. A missing residue remains absent; no connection is discarded.
    // Complexity review: one ordered chain/residue scan and at most two linear
    // atom scans, matching Model::find_cra/find_atom/find_by_element without
    // materializing an auxiliary hierarchy or atom index.
    let mut resolved = AtomAddress::new(
        address.chain_name(),
        address.residue(),
        "SG",
        Some(address.altloc()),
    );
    let target_residue = resolved.residue();
    let target_chain = resolved.chain_name();
    let residue = model
        .chains
        .iter()
        .filter(|chain| chain.source_id == target_chain)
        .flat_map(|chain| chain.residues.iter())
        .find(|candidate| {
            candidate.address.matches_without_segment(&target_residue)
                && candidate.address.segment() == target_residue.segment()
        });
    let Some(residue) = residue else {
        return (resolved, None);
    };

    let request = if resolved.altloc() == b'*' {
        AltLocRequest::Any
    } else {
        AltLocRequest::Exact((resolved.altloc() != 0).then(|| AltLocLabel::new(resolved.altloc())))
    };
    let atom = residue.atoms.iter().find(|atom| {
        read_string(atom.name.as_bytes()) == resolved.logical_atom_name().as_bytes()
            && cosmolkit_bio::altloc_matches(atom.altloc, request)
    });
    if atom.is_none_or(|atom| atom.element != Element::S)
        && let Some(sulfur) = residue.atoms.iter().find(|atom| atom.element == Element::S)
    {
        resolved = AtomAddress::new(
            resolved.chain_name(),
            resolved.residue(),
            std::str::from_utf8(read_string(sulfur.name.as_bytes()))
                .expect("PDB atom names are validated ASCII"),
            Some(sulfur.altloc.map_or(0, AltLocLabel::value)),
        );
    }
    (resolved, Some(residue))
}

fn complete_ssbond(
    connection: &mut BioConnection,
    model: &PdbModelGroup,
    crystal: &BioCrystalInfo,
) -> Result<(), PdbSsbondResolutionError> {
    // Gemmi source: third_party/gemmi/src/pdb.cpp::complete_ssbond.
    // Gemmi❗✔️: void complete_ssbond(Connection& con, const Model& mdl, const UnitCell& cell) {
    // Gemmi❗✔️:   const Residue* res1 = complete_ssbond_atom(con.partner1, mdl);
    // Gemmi❗✔️:   const Residue* res2 = complete_ssbond_atom(con.partner2, mdl);
    // Gemmi❗✔️:   if (res1 && res2 && (con.partner1.altloc != '\0' || con.partner2.altloc != '\0')) {
    // Gemmi❗✔️:     // pick a pair of atoms in the shortest distance
    // Gemmi❗✔️:     double min_dist_sq = INFINITY;
    // Gemmi❗✔️:     for (const Atom& a1 : const_cast<Residue*>(res1)->get(con.partner1.atom_name))
    // Gemmi❗✔️:       for (const Atom& a2 : const_cast<Residue*>(res2)->get(con.partner2.atom_name))
    // Gemmi❗✔️:         if (a2.same_conformer(a1)) {
    // Gemmi❗✔️:           double dist_sq = cell.find_nearest_image(a1.pos, a2.pos, con.asu).dist_sq;
    // Gemmi❗✔️:           if (dist_sq < min_dist_sq) {
    // Gemmi❗✔️:             con.partner1.altloc = a1.altloc;
    // Gemmi❗✔️:             con.partner2.altloc = a2.altloc;
    // Gemmi❗✔️:             min_dist_sq = dist_sq;
    // Gemmi❗✔️:           }
    // Gemmi❗✔️:         }
    // Gemmi❗✔️:   }
    // Gemmi❗✔️: }
    // Behavior review: both endpoints are independently resolved before the
    // altloc gate. Missing residues do not fail. If both exist and either
    // resolved altloc is nonblank, each source-name group must exist (partner
    // one first); absent groups preserve Gemmi's `No such atom` failure as a
    // typed error. All matching-name atoms participate in source order;
    // same-conformer filtering, nearest-image distance, strict minimum and
    // first-winner ties are preserved. Address updates are staged locally so a
    // returned error does not expose a partially mutated connection.
    // Complexity review: residue/address resolution is linear in the ordered
    // model hierarchy. Pair selection retains the source's repeated scan of
    // partner two for each partner-one atom, O(A1 + K1*A2), without allocating
    // atom-group copies; the lattice query has the already-audited fixed image
    // traversal cost.
    let (mut partner1, residue1) = complete_ssbond_atom(connection.partner1.clone(), model);
    let (mut partner2, residue2) = complete_ssbond_atom(connection.partner2.clone(), model);
    let mut selected_altlocs = None;

    if let (Some(residue1), Some(residue2)) = (residue1, residue2)
        && (partner1.altloc() != 0 || partner2.altloc() != 0)
    {
        if ssbond_atom_group(&partner1, residue1).next().is_none() {
            return Err(PdbSsbondResolutionError {
                partner: 1,
                atom_name: partner1.logical_atom_name().to_owned(),
            });
        }
        if ssbond_atom_group(&partner2, residue2).next().is_none() {
            return Err(PdbSsbondResolutionError {
                partner: 2,
                atom_name: partner2.logical_atom_name().to_owned(),
            });
        }

        let mut min_dist_sq = f64::INFINITY;
        for atom1 in ssbond_atom_group(&partner1, residue1) {
            for atom2 in ssbond_atom_group(&partner2, residue2) {
                if is_same_conformer(atom2.altloc, atom1.altloc) {
                    let dist_sq =
                        find_nearest_image(crystal, atom1.position, atom2.position, connection.asu)
                            .dist_sq();
                    if dist_sq < min_dist_sq {
                        selected_altlocs = Some((atom1.altloc, atom2.altloc));
                        min_dist_sq = dist_sq;
                    }
                }
            }
        }
    }

    if let Some((altloc1, altloc2)) = selected_altlocs {
        partner1 = AtomAddress::new(
            partner1.chain_name(),
            partner1.residue(),
            partner1.logical_atom_name(),
            Some(altloc1.map_or(0, AltLocLabel::value)),
        );
        partner2 = AtomAddress::new(
            partner2.chain_name(),
            partner2.residue(),
            partner2.logical_atom_name(),
            Some(altloc2.map_or(0, AltLocLabel::value)),
        );
    }
    connection.partner1 = partner1;
    connection.partner2 = partner2;
    Ok(())
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct PdbGroupLocation {
    model_id: BioModelId,
    chain_id: BioChainId,
    residue_id: BioResidueId,
    atom_index: usize,
}

#[derive(Debug, Clone, PartialEq, Default)]
struct PdbHierarchyGrouping {
    models: Vec<PdbModelGroup>,
    chain_row_count: usize,
    residue_row_count: usize,
}

impl PdbHierarchyGrouping {
    fn add_model(
        &mut self,
        source_model_number: Option<i32>,
    ) -> Result<BioModelId, BioStructureError> {
        // Gemmi source: third_party/gemmi/src/pdb.cpp, populate_structure_from_pdb_stream.
        // Gemmi❗✔️:           st.models.emplace_back(num);
        // Gemmi❗✔️:           model = &st.models.back();
        // Behavior review: source model numbers are stored in first-seen model-row order.
        // Complexity review: row append and id creation are amortized O(1).
        let index = self.models.len();
        let value = u32::try_from(index)
            .map_err(|_| BioStructureError::RowIndexTooLarge { value: index })?;
        let row_id = BioModelId::new(value);
        self.models.push(PdbModelGroup {
            row_id,
            source_model_number,
            chains: Vec::new(),
            current_chain: None,
            after_ter: false,
        });
        Ok(row_id)
    }

    fn add_atom_record(
        &mut self,
        model_id: BioModelId,
        chain_name: PdbChainId,
        residue: ResidueAddress,
        record_type: u8,
        fields: PdbAtomFields,
    ) -> Result<PdbGroupLocation, BioStructureError> {
        // Gemmi source: third_party/gemmi/src/pdb.cpp,
        // populate_structure_from_pdb_stream.
        // Gemmi❗✔️:       if (!chain || chain_name != chain->name) {
        // Gemmi❗✔️:         const Chain* prev_part = model->find_chain(chain_name);
        // Gemmi❗✔️:         after_ter = prev_part &&
        // Gemmi❗✔️:                     prev_part->residues[0].entity_type == EntityType::Polymer;
        // Gemmi❗✔️:         model->chains.emplace_back(chain_name);
        // Gemmi❗✔️:         chain = &model->chains.back();
        // Gemmi❗✔️:         resmap.clear();
        // Gemmi❗✔️:         resi = nullptr;
        // Gemmi❗✔️:       }
        // Gemmi source helper closure: third_party/gemmi/include/gemmi/model.hpp.
        // Gemmi❗✔️:   Chain* find_chain(const std::string& chain_name) {
        // Gemmi❗✔️:     return impl::find_or_null(chains, chain_name);
        // Gemmi❗✔️:   }
        // Gemmi❗✔️: T* find_or_null(std::vector<T>& vec, const S& name) {
        // Gemmi❗✔️:   auto it = find_iter_(vec, name);
        // Gemmi❗✔️:   return it != vec.end() ? &*it : nullptr;
        // Gemmi❗✔️: }
        // Gemmi❗✔️: auto find_iter_(Vec& vec, const S& name) {
        // Gemmi❗✔️:   return std::find_if(vec.begin(), vec.end(), [&name](const auto& m) { return get_id(m) == name; });
        // Gemmi❗✔️: }
        // Gemmi❗✔️: auto get_id(const T& m) -> decltype(m.name) { return m.name; }
        // Gemmi❗✔️:       if (!resi || !resi->matches(rid)) {
        // Gemmi❗✔️:         auto it = resmap.find(rid);
        // Gemmi❗✔️:         if (it == resmap.end()) {
        // Gemmi❗✔️:           resmap.emplace(rid, (int) chain->residues.size());
        // Gemmi❗✔️:           chain->residues.emplace_back(rid);
        // Gemmi❗✔️:           resi = &chain->residues.back();
        // Gemmi❗✔️:           resi->het_flag = line[0] & ~0x20;
        // Gemmi❗✔️:         } else {
        // Gemmi❗✔️:           resi = &chain->residues[it->second];
        // Gemmi❗✔️:         }
        // Gemmi❗✔️:       }
        // Behavior review: new chain parts preserve source order and recompute
        // after_ter from the first earlier same-named part, just as Gemmi's
        // first-match find_chain does. An empty prior chain would make Gemmi's
        // residues[0] undefined; parser-created chain parts always receive a
        // residue with their first atom, and this safe check does not claim an
        // output for a malformed internal empty-chain state. New residues
        // inherit Water/NonPolymer only after a prior Polymer part. The current
        // residue fast path and first-seen HET flag remain unchanged.
        // Complexity review: source current-chain checks are O(1); a part
        // transition performs one first-match linear scan, matching Gemmi's
        // find_chain/find_or_null/find_iter_ closure. The residue HashMap keeps
        // expected O(1) lookup and existing rows append once.
        let model_index = model_id.index();
        let Some(model) = self.models.get_mut(model_index) else {
            return Err(BioStructureError::RowReferenceOutOfBounds {
                table: "models",
                index: model_id.value(),
                table_len: self.models.len(),
            });
        };

        let chain_index = match model.current_chain {
            Some(index) if model.chains[index].source_id == chain_name => index,
            _ => {
                let after_ter = model
                    .chains
                    .iter()
                    .find(|chain| chain.source_id == chain_name)
                    .and_then(|chain| chain.residues.first())
                    .is_some_and(|residue| residue.entity_kind == EntityKind::Polymer);
                model.after_ter = after_ter;
                let index = self.chain_row_count;
                let value = u32::try_from(index)
                    .map_err(|_| BioStructureError::RowIndexTooLarge { value: index })?;
                let row_id = BioChainId::new(value);
                self.chain_row_count = index
                    .checked_add(1)
                    .ok_or(BioStructureError::RowIndexTooLarge { value: index })?;
                let chain_index = model.chains.len();
                model.chains.push(PdbChainGroup {
                    row_id,
                    source_id: chain_name,
                    residues: Vec::new(),
                    residue_indices: HashMap::new(),
                    current_residue: None,
                });
                model.current_chain = Some(chain_index);
                chain_index
            }
        };

        let chain = &mut model.chains[chain_index];
        let residue_index = match chain.current_residue {
            Some(index) if chain.residues[index].address.matches(&residue) => index,
            _ => {
                let key = PdbResidueMapKey(residue);
                match chain.residue_indices.get(&key).copied() {
                    Some(index) => index,
                    None => {
                        let index = self.residue_row_count;
                        let value = u32::try_from(index)
                            .map_err(|_| BioStructureError::RowIndexTooLarge { value: index })?;
                        let row_id = BioResidueId::new(value);
                        self.residue_row_count = index
                            .checked_add(1)
                            .ok_or(BioStructureError::RowIndexTooLarge { value: index })?;
                        let local_index = chain.residues.len();
                        let entity_kind = if model.after_ter {
                            if gemmi_pdb_residue_is_water(residue.name()) {
                                EntityKind::Water
                            } else {
                                EntityKind::NonPolymer
                            }
                        } else {
                            EntityKind::Unknown
                        };
                        chain.residues.push(PdbResidueGroup {
                            row_id,
                            address: residue,
                            het_flag: record_type & !0x20,
                            entity_kind,
                            subchain_id: String::new(),
                            atoms: Vec::new(),
                        });
                        chain.residue_indices.insert(key, local_index);
                        local_index
                    }
                }
            }
        };
        chain.current_residue = Some(residue_index);
        let residue_group = &mut chain.residues[residue_index];
        let atom_index = residue_group.atoms.len();
        residue_group.atoms.push(fields);
        Ok(PdbGroupLocation {
            model_id: model.row_id,
            chain_id: chain.row_id,
            residue_id: residue_group.row_id,
            atom_index,
        })
    }
}

fn pdb_chain_entity_and_subchain_status(chain: &PdbChainGroup) -> (bool, bool) {
    // Gemmi source: third_party/gemmi/src/polyheur.cpp::has_entity_types_and_subchains.
    // Gemmi❗✔️: static std::pair<bool,bool> has_entity_types_and_subchains(const Chain& chain) {
    // Gemmi❗✔️:   bool has_entity_types = true;
    // Gemmi❗✔️:   bool has_subchains = true;
    // Gemmi❗✔️:   for (const Residue& res : chain.residues) {
    // Gemmi❗✔️:     if (res.subchain.empty())
    // Gemmi❗✔️:       has_subchains = false;
    // Gemmi❗✔️:     if (res.entity_type == EntityType::Unknown)
    // Gemmi❗✔️:       has_entity_types = false;
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   return {has_entity_types, has_subchains};
    // Gemmi❗✔️: }
    // Behavior review: both booleans are accumulated across every residue;
    // an empty chain therefore retains the source's (true, true) identity and
    // is not assigned a generated subchain by the fixed PDB call.
    // Complexity review: one forward pass with two boolean updates per row,
    // O(residues) time and O(1) auxiliary space, matching the source.
    let mut has_entity_types = true;
    let mut has_subchains = true;
    for residue in &chain.residues {
        if residue.subchain_id.is_empty() {
            has_subchains = false;
        }
        if residue.entity_kind == EntityKind::Unknown {
            has_entity_types = false;
        }
    }
    (has_entity_types, has_subchains)
}

fn assign_pdb_subchain_names(
    chain: &mut PdbChainGroup,
    nonpolymer_counter: &mut i32,
) -> Result<(), PdbModelTransitionError> {
    // Gemmi source: third_party/gemmi/src/polyheur.cpp::assign_subchain_names.
    // Gemmi❗✔️: void assign_subchain_names(Chain& chain, int& nonpolymer_counter) {
    // Gemmi❗✔️:   for (Residue& res : chain.residues) {
    // Gemmi❗✔️:     res.subchain = chain.name;
    // Gemmi❗✔️:     // We'd use '-' as a separator (A-p or B-4 is more clear), but although
    // Gemmi❗✔️:     // such names are valid in mmCIF, OneDep refuses to accept them.
    // Gemmi❗✔️:     res.subchain += "x";
    // Gemmi❗✔️:     switch (res.entity_type) {
    // Gemmi❗✔️:       case EntityType::Polymer:
    // Gemmi❗✔️:         res.subchain += 'p';
    // Gemmi❗✔️:         break;
    // Gemmi❗✔️:       case EntityType::NonPolymer:
    // Gemmi❗✔️:         ++nonpolymer_counter;
    // Gemmi❗✔️:         // to keep the name short use base36 for 2+ digit numbers:
    // Gemmi❗✔️:         // 1, 2, ..., 9, 00, 01, ..., 09, 0A, 0B, ..., 0Z, 10, ...
    // Gemmi❗✔️:         if (nonpolymer_counter < 10) {
    // Gemmi❗✔️:           res.subchain += char('0' + nonpolymer_counter);
    // Gemmi❗✔️:         } else {
    // Gemmi❗✔️:           const char base36[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    // Gemmi❗✔️:           int n = nonpolymer_counter - 10;
    // Gemmi❗✔️:           if (n < 36)
    // Gemmi❗✔️:             res.subchain += '0';
    // Gemmi❗✔️:           size_t pos = res.subchain.size();
    // Gemmi❗✔️:           while (n != 0) {
    // Gemmi❗✔️:             res.subchain.insert(res.subchain.begin() + pos, base36[n % 36]);
    // Gemmi❗✔️:             n /= 36;
    // Gemmi❗✔️:           }
    // Gemmi❗✔️:         }
    // Gemmi❗✔️:         break;
    // Gemmi❗✔️:       case EntityType::Water:
    // Gemmi❗✔️:         res.subchain += 'w';
    // Gemmi❗✔️:         break;
    // Gemmi❗✔️:       // In the wwPDB branched are kept each in separate auth/label chain.
    // Gemmi❗✔️:       // So we have one subchain in chain.
    // Gemmi❗✔️:       case EntityType::Branched:
    // Gemmi❗✔️:         res.subchain += 'b';
    // Gemmi❗✔️:         break;
    // Gemmi❗✔️:       case EntityType::Unknown:
    // Gemmi❗✔️:         break;
    // Gemmi❗✔️:     }
    // Gemmi❗✔️:   }
    // Gemmi❗✔️: }
    // Behavior review: the PDB caller invokes this only when all residue
    // kinds are known. It replaces every subchain in that chain, uses one
    // model-local counter shared by repeated chain names, and follows the
    // source's unusual counter-10/base36 insertion order exactly. Signed-int
    // counter overflow is undefined upstream; this port returns a typed error
    // rather than wrapping or saturating.
    // Complexity review: each residue is visited once. Source String insert
    // at a fixed position can shift an O(log36 n) digit suffix; Rust uses the
    // same insertion position and work, with one generated String per row.
    let chain_name = chain.source_id;
    for residue in &mut chain.residues {
        let mut subchain = chain_name.as_str().to_owned();
        subchain.push('x');
        match residue.entity_kind {
            EntityKind::Polymer => subchain.push('p'),
            EntityKind::NonPolymer => {
                *nonpolymer_counter = (*nonpolymer_counter).checked_add(1).ok_or(
                    PdbModelTransitionError::NonpolymerSubchainCounterOutsideSourceIntRange {
                        chain: chain_name,
                    },
                )?;
                if *nonpolymer_counter < 10 {
                    subchain.push(char::from(b'0' + *nonpolymer_counter as u8));
                } else {
                    const BASE36: &[u8; 36] = b"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
                    let mut n = *nonpolymer_counter - 10;
                    if n < 36 {
                        subchain.push('0');
                    }
                    let position = subchain.len();
                    while n != 0 {
                        subchain.insert(position, char::from(BASE36[(n % 36) as usize]));
                        n /= 36;
                    }
                }
            }
            EntityKind::Water => subchain.push('w'),
            EntityKind::Branched => subchain.push('b'),
            EntityKind::Unknown => {}
        }
        residue.subchain_id = subchain;
    }
    Ok(())
}

fn assign_pdb_subchains(models: &mut [PdbModelGroup]) -> Result<(), PdbModelTransitionError> {
    // Gemmi source: third_party/gemmi/src/polyheur.cpp::assign_subchains.
    // Gemmi❗🔝: void assign_subchains(Structure& st, bool force, bool fail_if_unknown) {
    // Gemmi❗🔝:   for (Model& model : st.models) {
    // Gemmi❗🔝:     std::map<std::string, int> counters;
    // Gemmi❗🔝:     for (Chain& chain : model.chains) {
    // Gemmi❗🔝:       auto has = has_entity_types_and_subchains(chain);
    // Gemmi❗🔝:       if (force || !has.second) {
    // Gemmi❗🔝:         if (has.first)  // all chain's residues have known entity_type
    // Gemmi❗🔝:           assign_subchain_names(chain, counters[chain.name]);
    // Gemmi❗🔝:         else if (fail_if_unknown)
    // Gemmi❗🔝:           fail("assign_subchains(): missing entity_type in chain " + chain.name);
    // Gemmi❗🔝:       }
    // Gemmi❗🔝:     }
    // Gemmi❗🔝:   }
    // Gemmi❗🔝: }
    // Behavior review: Gemmi's PDB caller fixes force=false and
    // fail_if_unknown=false. Only chains lacking any subchain are considered;
    // fully typed chains are rebuilt, while a chain containing Unknown is
    // left unchanged. Counters reset per model and persist across same-name
    // chain parts. HashMap iteration is never observed; only exact-key counter
    // lookup matters.
    // Complexity review: the source std::map performs O(log c) key access and
    // stores allocated chain-name keys. This fixed Rust path uses a hash of the
    // bounded Copy PdbChainId for expected O(1) access without string-key
    // allocation; chain visitation and residue work preserve source order.
    for model in models {
        let mut counters = HashMap::<PdbChainId, i32>::new();
        for chain in &mut model.chains {
            let (has_entity_types, has_subchains) = pdb_chain_entity_and_subchain_status(chain);
            if !has_subchains && has_entity_types {
                let counter = counters.entry(chain.source_id).or_default();
                assign_pdb_subchain_names(chain, counter)?;
            }
        }
    }
    Ok(())
}

fn first_polymer_span_subchain(chain: &PdbChainGroup) -> Option<&str> {
    // Gemmi source: third_party/gemmi/include/gemmi/model.hpp::Chain::get_polymer.
    // Gemmi❗✔️:   ResidueSpan get_polymer() {
    // Gemmi❗✔️:     auto begin = residues.begin();
    // Gemmi❗✔️:     while (begin != residues.end() && begin->entity_type != EntityType::Polymer)
    // Gemmi❗✔️:       ++begin;
    // Gemmi❗✔️:     auto end = begin;
    // Gemmi❗✔️:     while (end != residues.end() && end->entity_type == EntityType::Polymer
    // Gemmi❗✔️:                                  && end->subchain == begin->subchain)
    // Gemmi❗✔️:       ++end;
    // Gemmi❗✔️:     return ResidueSpan(residues, &*begin, end - begin);
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   const std::string& subchain_id() const {
    // Gemmi❗✔️:     if (this->empty())
    // Gemmi❗✔️:       throw std::out_of_range("subchain_id(): empty span");
    // Gemmi❗✔️:     if (this->size() > 1 && this->front().subchain != this->back().subchain)
    // Gemmi❗✔️:       fail("subchain id varies in a residue span: ", this->front().subchain,
    // Gemmi❗✔️:            " vs ", this->back().subchain);
    // Gemmi❗✔️:     return this->begin()->subchain;
    // Gemmi❗✔️:   }
    // Behavior review: get_polymer selects only the first contiguous Polymer
    // span and stops at either a non-Polymer row or a changed subchain. The
    // subsequent truth test excludes an empty span; for a nonempty selected
    // span, all rows share its first subchain, so subchain_id's error branch is
    // unreachable. Rust returns None when no Polymer row exists rather than
    // evaluating Gemmi's `&*end` expression for that source-undefined case.
    // Complexity review: two bounded forward scans visit at most the chain's
    // residue rows once in total; no allocation or clone is required.
    let start = chain
        .residues
        .iter()
        .position(|residue| residue.entity_kind == EntityKind::Polymer)?;
    let subchain_id = chain.residues[start].subchain_id.as_str();
    let mut end = start;
    while let Some(residue) = chain.residues.get(end) {
        if residue.entity_kind != EntityKind::Polymer || residue.subchain_id != subchain_id {
            break;
        }
        end += 1;
    }
    (end > start).then_some(subchain_id)
}

fn associate_pdb_entity_subchains(models: &[PdbModelGroup], entities: &mut PdbEntityState) {
    // Gemmi source: third_party/gemmi/src/pdb.cpp::populate_structure_from_pdb_stream.
    // Gemmi❗✔️:   for (Chain& ch : st.models[0].chains)
    // Gemmi❗✔️:     if (Entity* entity = st.get_entity(ch.name))
    // Gemmi❗✔️:       if (auto polymer = ch.get_polymer())
    // Gemmi❗✔️:         entity->subchains.emplace_back(polymer.subchain_id());
    // Gemmi source helper closure: third_party/gemmi/include/gemmi/model.hpp.
    // Gemmi❗✔️:   Entity* get_entity(const std::string& ent_id) {
    // Gemmi❗✔️:     return impl::find_or_null(entities, ent_id);
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   T* find_or_null(std::vector<T>& vec, const S& name) {
    // Gemmi❗✔️:     auto it = find_iter_(vec, name);
    // Gemmi❗✔️:     return it != vec.end() ? &*it : nullptr;
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   auto find_iter_(Vec& vec, const S& name) {
    // Gemmi❗✔️:     return std::find_if(vec.begin(), vec.end(), [&name](const auto& m) { return get_id(m) == name; });
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   auto get_id(const T& m) -> decltype(m.name) { return m.name; }
    // Behavior review: association is limited to first-model chain order,
    // exact source-name entity lookup, and one first Polymer span per chain.
    // Every match appends (never sorts or deduplicates) the span ID, including
    // a source-empty ID for a nonempty Polymer span whose assignment was
    // skipped due to an Unknown residue elsewhere in that chain.
    // Complexity review: entity lookup is the same linear entity-vector scan
    // as Gemmi; polymer search is linear in each visited first-model chain.
    let Some(first_model) = models.first() else {
        return;
    };
    for chain in &first_model.chains {
        let Some(entity) = entities
            .entities
            .iter_mut()
            .find(|entity| entity.source_name.as_slice() == chain.source_id.as_bytes())
        else {
            continue;
        };
        if let Some(subchain_id) = first_polymer_span_subchain(chain) {
            entity.subchains.push(subchain_id.to_owned());
        }
    }
}

#[derive(Debug)]
enum PdbModelTransitionError {
    ModelWithoutEnd,
    DuplicateModelNumber { number: i32 },
    AtomBetweenModels,
    ModelNumberOutsideSourceIntRange,
    NonpolymerSubchainCounterOutsideSourceIntRange { chain: PdbChainId },
    AnisouWithoutAtom,
    DuplicateAnisou,
    AnisouIntegerOutsideSourceIntRange { offset: usize },
    ModelIdOutOfBounds { id: BioModelId, model_count: usize },
    Structure(BioStructureError),
}

impl std::fmt::Display for PdbModelTransitionError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ModelWithoutEnd => formatter.write_str("MODEL without ENDMDL?"),
            Self::DuplicateModelNumber { number } => {
                write!(formatter, "duplicate MODEL number: {number}")
            }
            Self::AtomBetweenModels => formatter.write_str("ATOM/HETATM between models"),
            Self::ModelNumberOutsideSourceIntRange => {
                formatter.write_str("PDB MODEL number is outside the source integer range")
            }
            Self::NonpolymerSubchainCounterOutsideSourceIntRange { chain } => write!(
                formatter,
                "PDB nonpolymer subchain counter for chain {:?} is outside the source integer range",
                chain.as_str()
            ),
            Self::AnisouWithoutAtom => {
                formatter.write_str("ANISOU record not directly after ATOM/HETATM.")
            }
            Self::DuplicateAnisou => {
                formatter.write_str("Duplicated ANISOU record or not directly after ATOM/HETATM.")
            }
            Self::AnisouIntegerOutsideSourceIntRange { offset } => write!(
                formatter,
                "ANISOU integer field at byte offset {offset} is outside the source integer range"
            ),
            Self::ModelIdOutOfBounds { id, model_count } => write!(
                formatter,
                "PDB active model {} is outside {} model rows",
                id.value(),
                model_count
            ),
            Self::Structure(error) => error.fmt(formatter),
        }
    }
}

impl std::error::Error for PdbModelTransitionError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Structure(error) => Some(error),
            _ => None,
        }
    }
}

impl From<BioStructureError> for PdbModelTransitionError {
    fn from(error: BioStructureError) -> Self {
        Self::Structure(error)
    }
}

#[derive(Debug, Clone, PartialEq, Default)]
struct PdbModelTransition {
    grouping: PdbHierarchyGrouping,
    active_model: Option<BioModelId>,
    stopped_at_end_record: bool,
}

impl PdbModelTransition {
    fn model_record(
        &mut self,
        source_line_buffer: &[u8; 122],
    ) -> Result<BioModelId, PdbModelTransitionError> {
        // Gemmi source: third_party/gemmi/src/pdb.cpp,
        // populate_structure_from_pdb_stream.
        // Gemmi❗✔️:     } else if (is_record_type4(line, "MODEL")) {
        // Gemmi❗✔️:       if (model && chain)
        // Gemmi❗✔️:         wrong("MODEL without ENDMDL?");
        // Gemmi❗✔️:       int num = read_int(line+6, 8);
        // Gemmi❗✔️:       model = &st.find_or_add_model(num);
        // Gemmi❗✔️:       if (!model->chains.empty())
        // Gemmi❗✔️:         wrong("duplicate MODEL number: " + std::to_string(num));
        // Gemmi❗✔️:       chain = nullptr;
        // Behavior review: lookup is by exact source model number, preserving
        // first model-row order; an existing empty model is reusable, while a
        // populated duplicate is rejected. The active model changes only after
        // source checks succeed. Fixed-width read_int overflow is not source-
        // defined and is returned structurally rather than assigned a value.
        // Complexity review: vector position lookup scans the model rows once,
        // matching Gemmi's linear find_or_null/find_or_add; append is amortized
        // O(1), with no per-record model-number allocation.
        if let Some(active_id) = self.active_model {
            let active = self.grouping.models.get(active_id.index()).ok_or(
                PdbModelTransitionError::ModelIdOutOfBounds {
                    id: active_id,
                    model_count: self.grouping.models.len(),
                },
            )?;
            if active.current_chain.is_some() {
                return Err(PdbModelTransitionError::ModelWithoutEnd);
            }
        }

        let number = read_int(&source_line_buffer[6..14])
            .ok_or(PdbModelTransitionError::ModelNumberOutsideSourceIntRange)?;
        let model_id = if let Some(model) = self
            .grouping
            .models
            .iter()
            .find(|model| model.source_model_number == Some(number))
        {
            if !model.chains.is_empty() {
                return Err(PdbModelTransitionError::DuplicateModelNumber { number });
            }
            model.row_id
        } else {
            self.grouping.add_model(Some(number))?
        };
        let Some(model) = self.grouping.models.get_mut(model_id.index()) else {
            return Err(PdbModelTransitionError::ModelIdOutOfBounds {
                id: model_id,
                model_count: self.grouping.models.len(),
            });
        };
        model.current_chain = None;
        self.active_model = Some(model_id);
        Ok(model_id)
    }

    fn end_model_record(&mut self) {
        // Gemmi source: third_party/gemmi/src/pdb.cpp,
        // populate_structure_from_pdb_stream.
        // Gemmi❗✔️:     } else if (is_record_type4(line, "ENDMDL")) {
        // Gemmi❗✔️:       model = nullptr;
        // Gemmi❗✔️:       chain = nullptr;
        // Behavior review: ENDMDL is accepted even without an active model;
        // the parser's next ATOM row follows the implicit-model source branch.
        // Complexity review: only two optional state slots are cleared.
        if let Some(model_id) = self.active_model
            && let Some(model) = self.grouping.models.get_mut(model_id.index())
        {
            model.current_chain = None;
        }
        self.active_model = None;
    }

    fn end_record(&mut self) {
        // Gemmi source: third_party/gemmi/src/pdb.cpp,
        // populate_structure_from_pdb_stream.
        // Gemmi❗✔️:     } else if (is_record_type3(line, "END")) {
        // Gemmi❗✔️:       break;
        // Behavior review: the record terminates dispatch without changing
        // already accumulated rows or the active model's source state.
        // Complexity review: a single flag write is O(1).
        self.stopped_at_end_record = true;
    }

    fn anisou_record(
        &mut self,
        source_line_buffer: &[u8; 122],
    ) -> Result<(), PdbModelTransitionError> {
        // Gemmi source: third_party/gemmi/src/pdb.cpp,
        // populate_structure_from_pdb_stream.
        // Gemmi❗✔️:     } else if (is_record_type4(line, "ANISOU")) {
        // Gemmi❗✔️:       if (!model || !chain || !resi || resi->atoms.empty())
        // Gemmi❗✔️:         wrong("ANISOU record not directly after ATOM/HETATM.");
        // Gemmi❗✔️:       // We assume that ANISOU refers to the last atom.
        // Gemmi❗✔️:       // Can it not be the case?
        // Gemmi❗✔️:       Atom &atom = resi->atoms.back();
        // Gemmi❗✔️:       if (atom.aniso.u11 != 0.)
        // Gemmi❗✔️:         wrong("Duplicated ANISOU record or not directly after ATOM/HETATM.");
        // Gemmi❗✔️:       atom.aniso.u11 = read_int(line+28, 7) * 1e-4f;
        // Gemmi❗✔️:       atom.aniso.u22 = read_int(line+35, 7) * 1e-4f;
        // Gemmi❗✔️:       atom.aniso.u33 = read_int(line+42, 7) * 1e-4f;
        // Gemmi❗✔️:       atom.aniso.u12 = read_int(line+49, 7) * 1e-4f;
        // Gemmi❗✔️:       atom.aniso.u13 = read_int(line+56, 7) * 1e-4f;
        // Gemmi❗✔️:       atom.aniso.u23 = read_int(line+63, 7) * 1e-4f;
        // Behavior review: association follows the active model, its current
        // chain and residue, and that residue's final atom. The duplicate
        // sentinel is specifically prior u11; zero u11 allows a later ANISOU
        // row to overwrite. `read_int` reproduces the source's width-limited
        // prefix/no-conversion conversion. Every value is rounded as f32 and
        // then widened into the existing f64 IO carrier because Gemmi stores
        // the six components as `float`. A seven-byte field cannot overflow
        // signed i32; the structured range variant is unreachable on this row
        // format and avoids inventing an overflow value if assumptions change.
        // Complexity review: hierarchy navigation and last-element access are
        // O(1); six scans each inspect at most seven bytes, with no allocation.
        let model_id = self
            .active_model
            .ok_or(PdbModelTransitionError::AnisouWithoutAtom)?;
        let model_index = model_id.index();
        let model = self.grouping.models.get(model_index).ok_or(
            PdbModelTransitionError::ModelIdOutOfBounds {
                id: model_id,
                model_count: self.grouping.models.len(),
            },
        )?;
        let chain_index = model
            .current_chain
            .ok_or(PdbModelTransitionError::AnisouWithoutAtom)?;
        let chain = model
            .chains
            .get(chain_index)
            .ok_or(PdbModelTransitionError::AnisouWithoutAtom)?;
        let residue_index = chain
            .current_residue
            .ok_or(PdbModelTransitionError::AnisouWithoutAtom)?;
        let residue = chain
            .residues
            .get(residue_index)
            .ok_or(PdbModelTransitionError::AnisouWithoutAtom)?;
        let atom = residue
            .atoms
            .last()
            .ok_or(PdbModelTransitionError::AnisouWithoutAtom)?;
        if atom.anisou[0] != 0.0 {
            return Err(PdbModelTransitionError::DuplicateAnisou);
        }

        let mut anisou = [0.0; 6];
        for (component, offset) in anisou.iter_mut().zip([28_usize, 35, 42, 49, 56, 63]) {
            let value = read_int(&source_line_buffer[offset..offset + 7])
                .ok_or(PdbModelTransitionError::AnisouIntegerOutsideSourceIntRange { offset })?;
            *component = f64::from(value as f32 * 1e-4_f32);
        }

        let atom = self
            .grouping
            .models
            .get_mut(model_index)
            .and_then(|model| model.chains.get_mut(chain_index))
            .and_then(|chain| chain.residues.get_mut(residue_index))
            .and_then(|residue| residue.atoms.last_mut())
            .ok_or(PdbModelTransitionError::AnisouWithoutAtom)?;
        atom.anisou = anisou;
        Ok(())
    }

    fn ter_record(
        &mut self,
        line: &[u8],
        options: PdbInputOptions,
        source_state: &mut BioStructureSourceState,
    ) -> Result<(), PdbModelTransitionError> {
        // Gemmi source: third_party/gemmi/src/pdb.cpp,
        // populate_structure_from_pdb_stream.
        // Gemmi❗✔️:     } else if (is_record_type3(line, "TER") && !options.ignore_ter) {
        // Gemmi❗✔️:       if (!chain || st.ter_status == 'e')
        // Gemmi❗✔️:         continue;
        // Gemmi❗✔️:       st.ter_status = 'y';
        // Gemmi❗✔️:       if (options.split_chain_on_ter) {
        // Gemmi❗✔️:         chain = nullptr;
        // Gemmi❗✔️:         // split_chain_on_ter is used for AMBER files that can have TER records
        // Gemmi❗✔️:         // in various places. So in such case TER doesn't imply entity_type.
        // Gemmi❗✔️:         continue;
        // Gemmi❗✔️:       }
        // Gemmi❗✔️:       // If we have 2+ TER records in one chain, they are used in non-standard
        // Gemmi❗✔️:       // way and should be better ignored (in all the chains).
        // Gemmi❗✔️:       if (after_ter) {
        // Gemmi❗✔️:         st.ter_status = 'e';  // all entity_types will be later set to Unknown
        // Gemmi❗✔️:         continue;
        // Gemmi❗✔️:       }
        // Gemmi❗✔️:       for (Residue& res : chain->residues) {
        // Gemmi❗✔️:         res.entity_type = EntityType::Polymer;
        // Gemmi❗✔️:         // Sanity check: water should not be marked as a polymer.
        // Gemmi❗✔️:         if GEMMI_UNLIKELY(res.is_water())
        // Gemmi❗✔️:           st.ter_status = 'e';  // all entity_types will be later set to Unknown
        // Gemmi❗✔️:       }
        // Gemmi❗✔️:       after_ter = true;
        // Gemmi❗✔️:     }
        // Behavior review: record matching, ignore, absent-chain and prior-error
        // paths are no-ops. A handled record makes the source status durable as
        // `y`; split mode clears only the current chain pointer. In normal mode
        // an already-after-TER chain sets `e`; otherwise every current residue
        // is marked Polymer, while a water residue sets `e` and finalization
        // later clears all entity kinds. The state changes follow source order.
        // Complexity review: record matching and status branches are O(1); a
        // normal TER visits each residue of the current chain once, as source.
        if !gemmi_record_type3(line, 0, *b"TER\0") || options.ignore_ter {
            return Ok(());
        }
        if source_state.ter_status == b'e' {
            return Ok(());
        }

        let Some(model_id) = self.active_model else {
            return Ok(());
        };
        let model_index = model_id.index();
        let Some(model) = self.grouping.models.get_mut(model_index) else {
            return Err(PdbModelTransitionError::ModelIdOutOfBounds {
                id: model_id,
                model_count: self.grouping.models.len(),
            });
        };
        let Some(chain_index) = model.current_chain else {
            return Ok(());
        };
        if model.chains.get(chain_index).is_none() {
            let index = u32::try_from(chain_index)
                .map_err(|_| BioStructureError::RowIndexTooLarge { value: chain_index })?;
            return Err(BioStructureError::RowReferenceOutOfBounds {
                table: "chains",
                index,
                table_len: model.chains.len(),
            }
            .into());
        }

        source_state.ter_status = b'y';
        if options.split_chain_on_ter {
            model.current_chain = None;
            return Ok(());
        }
        if model.after_ter {
            source_state.ter_status = b'e';
            return Ok(());
        }

        for residue in &mut model.chains[chain_index].residues {
            residue.entity_kind = EntityKind::Polymer;
            if gemmi_pdb_residue_is_water(residue.address.name()) {
                source_state.ter_status = b'e';
            }
        }
        model.after_ter = true;
        Ok(())
    }

    fn finalize_ter_entity_types(&mut self, source_state: &BioStructureSourceState) {
        // Gemmi source: third_party/gemmi/src/pdb.cpp and polyheur.cpp.
        // Gemmi❗✔️:   if (st.ter_status == 'e')
        // Gemmi❗✔️:     remove_entity_types(st);
        // Gemmi❗✔️: void remove_entity_types(Structure& st) {
        // Gemmi❗✔️:   for (Model& model : st.models)
        // Gemmi❗✔️:     for (Chain& chain : model.chains)
        // Gemmi❗✔️:       for (Residue& res : chain.residues)
        // Gemmi❗✔️:         res.entity_type = EntityType::Unknown;
        // Gemmi❗✔️: }
        // Behavior review: only durable error status `e` rolls back every
        // residue's entity kind; `0` and `y` preserve their accumulated state.
        // Complexity review: one nested visit over model/chain/residue rows,
        // matching the source cleanup traversal without allocating.
        if source_state.ter_status == b'e' {
            for model in &mut self.grouping.models {
                for chain in &mut model.chains {
                    for residue in &mut chain.residues {
                        residue.entity_kind = EntityKind::Unknown;
                    }
                }
            }
        }
    }

    fn finalize_subchains_and_entities(
        &mut self,
        source_state: &BioStructureSourceState,
        entities: &mut PdbEntityState,
    ) -> Result<(), PdbModelTransitionError> {
        // Gemmi source: third_party/gemmi/src/pdb.cpp::populate_structure_from_pdb_stream.
        // Gemmi❗✔️:   if (st.ter_status == 'e')
        // Gemmi❗✔️:     remove_entity_types(st);
        // Gemmi❗✔️:
        // Gemmi❗✔️:   // Here we assign Residue::subchain, but only for chains with all
        // Gemmi❗✔️:   // Residue::entity_type assigned, i.e. for chains with TER.
        // Gemmi❗✔️:   assign_subchains(st, /*force=*/false, /*fail_if_unknown=*/false);
        // Gemmi❗✔️:
        // Gemmi❗✔️:   for (Chain& ch : st.models[0].chains)
        // Gemmi❗✔️:     if (Entity* entity = st.get_entity(ch.name))
        // Gemmi❗✔️:       if (auto polymer = ch.get_polymer())
        // Gemmi❗✔️:         entity->subchains.emplace_back(polymer.subchain_id());
        // Behavior review: this wrapper preserves the source ordering by
        // clearing TER-inferred kinds first, assigning names next, and finally
        // appending first-model polymer subchains to exact-name entities. The
        // stream finalizer calls it after default-model insertion. Existing
        // entity subchain entries remain and source order/duplicates append.
        // Complexity review: rollback, assignment and association each scan
        // only their required rows; no whole-structure copy or sorting occurs.
        self.finalize_ter_entity_types(source_state);
        assign_pdb_subchains(&mut self.grouping.models)?;
        associate_pdb_entity_subchains(&self.grouping.models, entities);
        Ok(())
    }

    fn stopped_at_end_record(&self) -> bool {
        self.stopped_at_end_record
    }

    fn implicit_model_for_atom(&mut self) -> Result<BioModelId, PdbModelTransitionError> {
        // Gemmi source: third_party/gemmi/src/pdb.cpp,
        // populate_structure_from_pdb_stream.
        // Gemmi❗✔️:         if (!model) {
        // Gemmi❗✔️:           int num = (int) st.models.size() + 1;
        // Gemmi❗✔️:           if (st.find_model(num))
        // Gemmi❗✔️:             wrong("ATOM/HETATM between models");
        // Gemmi❗✔️:           st.models.emplace_back(num);
        // Gemmi❗✔️:           model = &st.models.back();
        // Gemmi❗✔️:         }
        // Behavior review: an active explicit/implicit model is reused. With
        // no active model, the candidate is current model count plus one; any
        // existing row with that number causes the source error regardless of
        // whether that row contains chains. The fixed-width model count cannot
        // be represented by an invented wraparound.
        // Complexity review: one linear model-number search matches Gemmi's
        // find_model; adding the row is amortized O(1).
        if let Some(model_id) = self.active_model {
            if self.grouping.models.get(model_id.index()).is_none() {
                return Err(PdbModelTransitionError::ModelIdOutOfBounds {
                    id: model_id,
                    model_count: self.grouping.models.len(),
                });
            }
            return Ok(model_id);
        }

        let next_number = self
            .grouping
            .models
            .len()
            .checked_add(1)
            .and_then(|number| i32::try_from(number).ok())
            .ok_or(PdbModelTransitionError::ModelNumberOutsideSourceIntRange)?;
        if self
            .grouping
            .models
            .iter()
            .any(|model| model.source_model_number == Some(next_number))
        {
            return Err(PdbModelTransitionError::AtomBetweenModels);
        }
        let model_id = self.grouping.add_model(Some(next_number))?;
        self.active_model = Some(model_id);
        Ok(model_id)
    }

    fn add_atom_record(
        &mut self,
        chain_name: PdbChainId,
        residue: ResidueAddress,
        record_type: u8,
        fields: PdbAtomFields,
    ) -> Result<PdbGroupLocation, PdbModelTransitionError> {
        let model_id = self.implicit_model_for_atom()?;
        self.grouping
            .add_atom_record(model_id, chain_name, residue, record_type, fields)
            .map_err(Into::into)
    }

    fn finish_eof(&mut self) -> Result<(), PdbModelTransitionError> {
        // Gemmi source: third_party/gemmi/src/pdb.cpp,
        // populate_structure_from_pdb_stream.
        // Gemmi❗✔️:   if (st.models.empty())
        // Gemmi❗✔️:     st.models.emplace_back(1);
        // Behavior review: only a wholly model-free input gets the source's
        // empty model 1; explicit or implicit prior models remain untouched.
        // Complexity review: checking vector emptiness and one append are O(1).
        if self.grouping.models.is_empty() {
            self.grouping.add_model(Some(1))?;
        }
        Ok(())
    }

    fn wrong_input_format(
        &self,
        reader: &PdbReaderState,
        source_line_buffer: &[u8; 122],
        source: &str,
    ) -> Option<PdbInputFormatError> {
        reader.wrong_input_format(source_line_buffer, self.active_model.is_some(), source)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PdbRemarkFinalizationError {
    IntegerOutsideSourceDefinedRange,
    Remark3(PdbRemark3Error),
    Remark200(PdbRemark200Error),
    RemarkMetadata(PdbRemarkMetadataError),
    Remark350(PdbRemark350Error),
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum PdbSsbondRecordError {
    EmptyRecord,
    ChainNameNotAscii {
        partner: u8,
        byte: u8,
    },
    ResidueAddress {
        partner: u8,
        source: PdbResidueKeyError,
    },
    CounterOutsideSourceIntRange,
}

#[derive(Debug, Clone, PartialEq)]
enum PdbConnectionFinalizationError {
    SsbondRecord(PdbSsbondRecordError),
    LinkRecord(PdbLinkRecordError),
    CisPepRecord(PdbCisPepRecordError),
    SsbondResolution(PdbSsbondResolutionError),
}

#[derive(Debug)]
enum PdbStructureFinalizationError {
    Model(PdbModelTransitionError),
    Connection(PdbConnectionFinalizationError),
    Remarks(PdbRemarkFinalizationError),
    Ccd(PdbResidueKeyError),
    NoModelAfterDefault,
    SourceIdNotRepresentable {
        field: &'static str,
    },
    Utf8 {
        field: &'static str,
        valid_up_to: usize,
    },
    Structure(BioStructureError),
    StagingRowOrder {
        table: &'static str,
        expected: usize,
        actual: usize,
    },
}

impl From<BioStructureError> for PdbStructureFinalizationError {
    fn from(error: BioStructureError) -> Self {
        Self::Structure(error)
    }
}

impl From<PdbModelTransitionError> for PdbStructureFinalizationError {
    fn from(error: PdbModelTransitionError) -> Self {
        Self::Model(error)
    }
}

fn decode_pdb_ssbond_record(
    record: &[u8],
    disulf_count: &mut i32,
    first_model: &PdbModelGroup,
    crystal: &BioCrystalInfo,
) -> Result<Option<BioConnection>, PdbConnectionFinalizationError> {
    // Gemmi source: third_party/gemmi/src/pdb.cpp::process_conn SSBOND branch.
    // Gemmi❗✔️:     if (record[0] == 'S' || record[0] == 's') { // SSBOND
    // Gemmi❗✔️:       if (record.length() < 32)
    // Gemmi❗✔️:         continue;
    // Gemmi❗✔️:       Connection c;
    // Gemmi❗✔️:       c.name = "disulf" + std::to_string(++disulf_count);
    // Gemmi❗✔️:       c.type = Connection::Disulf;
    // Gemmi❗✔️:       const char* r = record.c_str();
    // Gemmi❗✔️:       c.partner1.chain_name = read_string(r + 14, 2);
    // Gemmi❗✔️:       c.partner1.res_id = read_res_id(r + 17, r + 11);
    // Gemmi❗✔️:       c.partner2.chain_name = read_string(r + 28, 2);
    // Gemmi❗✔️:       char res_id2[5] = {' ', ' ', ' ', ' ', ' '};
    // Gemmi❗✔️:       std::memcpy(res_id2, r + 31, std::min((size_t)5, record.length() - 31));
    // Gemmi❗✔️:       c.partner2.res_id = read_res_id(res_id2, r + 25);
    // Gemmi❗✔️:       c.asu = compare_link_symops(record, c.reported_sym);
    // Gemmi❗✔️:       if (record.length() > 73)
    // Gemmi❗✔️:         c.reported_distance = read_double(r + 73, 5);
    // Gemmi❗✔️:       complete_ssbond(c, st.first_model(), st.cell);
    // Gemmi❗✔️:       st.connections.emplace_back(c);
    // Gemmi❗✔️:     }
    // Behavior review: source order/counter, short-row skip, NUL-terminated
    // fixed fields, the padded second sequence buffer, ASU/symmetry parsing,
    // optional five-byte distance and first-model SG completion are retained.
    // Signed counter overflow and values outside the existing typed address
    // representation return structured boundaries rather than invented state.
    // Complexity review: fixed-width decoding is O(1); SSBOND atom completion
    // performs only the already source-shaped ordered hierarchy scans.
    if record.first().is_none() {
        return Err(PdbConnectionFinalizationError::SsbondRecord(
            PdbSsbondRecordError::EmptyRecord,
        ));
    }
    if record.len() < 32 {
        return Ok(None);
    }

    *disulf_count = disulf_count.checked_add(1).ok_or_else(|| {
        PdbConnectionFinalizationError::SsbondRecord(
            PdbSsbondRecordError::CounterOutsideSourceIntRange,
        )
    })?;
    let decode_chain = |partner: u8, offset: usize| {
        let bytes = read_string(&record[offset..offset + 2]);
        PdbChainId::from_ascii(bytes).ok_or_else(|| {
            PdbConnectionFinalizationError::SsbondRecord(PdbSsbondRecordError::ChainNameNotAscii {
                partner,
                byte: bytes
                    .iter()
                    .copied()
                    .find(|byte| !byte.is_ascii())
                    .expect("PdbChainId::from_ascii rejects only non-ASCII bytes"),
            })
        })
    };
    let decode_residue = |partner: u8,
                          sequence: &[u8; 5],
                          name_offset: usize|
     -> Result<ResidueAddress, PdbConnectionFinalizationError> {
        let name: &[u8; 3] = record[name_offset..name_offset + 3]
            .try_into()
            .expect("the SSBOND length gate covers fixed residue names");
        read_res_id(sequence, name).map_err(|source| {
            PdbConnectionFinalizationError::SsbondRecord(PdbSsbondRecordError::ResidueAddress {
                partner,
                source,
            })
        })
    };

    let chain1 = decode_chain(1, 14)?;
    let sequence1: &[u8; 5] = record[17..22]
        .try_into()
        .expect("the SSBOND length gate covers the first sequence id");
    let residue1 = decode_residue(1, sequence1, 11)?;
    let chain2 = decode_chain(2, 28)?;
    let mut sequence2 = [b' '; 5];
    let sequence2_len = (record.len() - 31).min(sequence2.len());
    sequence2[..sequence2_len].copy_from_slice(&record[31..31 + sequence2_len]);
    let residue2 = decode_residue(2, &sequence2, 25)?;

    let partner1 = AtomAddress::new(chain1, residue1, "", None);
    let partner2 = AtomAddress::new(chain2, residue2, "", None);
    let mut reported_sym = [0; 4];
    let mut connection = BioConnection {
        name: format!("disulf{disulf_count}"),
        kind: BioConnectionKind::Disulf,
        asu: compare_link_symops(record, &mut reported_sym),
        partner1,
        partner2,
        reported_sym,
        ..BioConnection::default()
    };
    if record.len() > 73 {
        let mut distance_field = [0; 5];
        let available = record.len().saturating_sub(73).min(distance_field.len());
        distance_field[..available].copy_from_slice(&record[73..73 + available]);
        connection.reported_distance = read_double(&distance_field);
    }
    complete_ssbond(&mut connection, first_model, crystal)
        .map_err(|error| PdbConnectionFinalizationError::SsbondResolution(error))?;
    Ok(Some(connection))
}

fn process_pdb_connection_records(
    records: &[Vec<u8>],
    first_model: &PdbModelGroup,
    models: &[PdbModelGroup],
    crystal: &BioCrystalInfo,
) -> Result<(Vec<BioConnection>, Vec<BioCisPep>), PdbConnectionFinalizationError> {
    // Gemmi source: third_party/gemmi/src/pdb.cpp::process_conn.
    // Gemmi❗✔️: void process_conn(Structure& st, const std::vector<std::string>& conn_records) {
    // Gemmi❗✔️:   int disulf_count = 0;
    // Gemmi❗✔️:   int covale_count = 0;
    // Gemmi❗✔️:   int metalc_count = 0;
    // Gemmi❗✔️:   for (const std::string& record : conn_records) {
    // Gemmi❗✔️:     if (record[0] == 'S' || record[0] == 's') { // SSBOND
    // Gemmi❗✔️:       if (record.length() < 32)
    // Gemmi❗✔️:         continue;
    // Gemmi❗✔️:       Connection c;
    // Gemmi❗✔️:       c.name = "disulf" + std::to_string(++disulf_count);
    // Gemmi❗✔️:       c.type = Connection::Disulf;
    // Gemmi❗✔️:       const char* r = record.c_str();
    // Gemmi❗✔️:       c.partner1.chain_name = read_string(r + 14, 2);
    // Gemmi❗✔️:       c.partner1.res_id = read_res_id(r + 17, r + 11);
    // Gemmi❗✔️:       c.partner2.chain_name = read_string(r + 28, 2);
    // Gemmi❗✔️:       char res_id2[5] = {' ', ' ', ' ', ' ', ' '};
    // Gemmi❗✔️:       std::memcpy(res_id2, r + 31, std::min((size_t)5, record.length() - 31));
    // Gemmi❗✔️:       c.partner2.res_id = read_res_id(res_id2, r + 25);
    // Gemmi❗✔️:       c.asu = compare_link_symops(record, c.reported_sym);
    // Gemmi❗✔️:       if (record.length() > 73)
    // Gemmi❗✔️:         c.reported_distance = read_double(r + 73, 5);
    // Gemmi❗✔️:       complete_ssbond(c, st.first_model(), st.cell);
    // Gemmi❗✔️:       st.connections.emplace_back(c);
    // Gemmi❗✔️:     } else if (record[0] == 'L' || record[0] == 'l') { // LINK
    // Gemmi❗✔️:       if (record.length() < 57)
    // Gemmi❗✔️:         continue;
    // Gemmi❗✔️:       Connection c;
    // Gemmi❗✔️:       for (int i : {0, 1}) {
    // Gemmi❗✔️:         const char* t = record.c_str() + 30 * i;
    // Gemmi❗✔️:         AtomAddress& ad = (i == 0 ? c.partner1 : c.partner2);
    // Gemmi❗✔️:         ad.chain_name = read_string(t + 20, 2);
    // Gemmi❗✔️:         ad.res_id = read_res_id(t + 22, t + 17);
    // Gemmi❗✔️:         ad.atom_name = read_string(t + 12, 4);
    // Gemmi❗✔️:         ad.altloc = read_altloc(t[16]);
    // Gemmi❗✔️:       }
    // Gemmi❗✔️:       auto get_elem = [&](const char* name, const AtomAddress& ad) {
    // Gemmi❗✔️:         if (element_from_padded_name_is_ambiguous(name)) {
    // Gemmi❗✔️:           const_CRA cra = st.first_model().find_cra(ad);
    // Gemmi❗✔️:           if (cra.atom)
    // Gemmi❗✔️:             return cra.atom->element.elem;
    // Gemmi❗✔️:         }
    // Gemmi❗✔️:         return infer_element_from_padded_name(name);
    // Gemmi❗✔️:       };
    // Gemmi❗✔️:       // emulating names used in wwPDB mmCIFs (covaleN and metalcN)
    // Gemmi❗✔️:       if (is_metal(get_elem(&record[12], c.partner1)) ||
    // Gemmi❗✔️:           is_metal(get_elem(&record[42], c.partner2))) {
    // Gemmi❗✔️:         c.name = "metalc" + std::to_string(++metalc_count);
    // Gemmi❗✔️:         c.type = Connection::MetalC;
    // Gemmi❗✔️:       } else {
    // Gemmi❗✔️:         c.name = "covale" + std::to_string(++covale_count);
    // Gemmi❗✔️:         c.type = Connection::Covale;
    // Gemmi❗✔️:       }
    // Gemmi❗✔️:       c.asu = compare_link_symops(record, c.reported_sym);
    // Gemmi❗✔️:       if (record.length() > 73) {
    // Gemmi❗✔️:         if (record[4] == 'R')
    // Gemmi❗✔️:           c.link_id = read_string(&record[72], 8);
    // Gemmi❗✔️:         else
    // Gemmi❗✔️:           c.reported_distance = read_double(&record[73], 5);
    // Gemmi❗✔️:       }
    // Gemmi❗✔️:       st.connections.emplace_back(c);
    // Gemmi❗✔️:     } else if (record[0] == 'C' || record[0] == 'c') { // CISPEP
    // Gemmi❗✔️:       if (record.length() < 22)
    // Gemmi❗✔️:         continue;
    // Gemmi❗✔️:       const char* r = record.c_str();
    // Gemmi❗✔️:       CisPep cispep;
    // Gemmi❗✔️:       cispep.partner_c.chain_name = read_string(r + 14, 2);
    // Gemmi❗✔️:       cispep.partner_c.res_id = read_res_id(r + 17, r + 11);
    // Gemmi❗✔️:       cispep.partner_n.chain_name = read_string(r + 28, 2);
    // Gemmi❗✔️:       cispep.partner_n.res_id = read_res_id(r + 31, r + 25);
    // Gemmi❗✔️:       // In files with a single model in the PDB CISPEP modNum is 0,
    // Gemmi❗✔️:       // but _struct_mon_prot_cis.pdbx_PDB_model_num is 1.
    // Gemmi❗✔️:       cispep.model_num = st.models.size() == 1 ? st.models[0].num : read_int(r + 43, 3);
    // Gemmi❗✔️:       cispep.reported_angle = read_double(r + 53, 6);
    // Gemmi❗✔️:       st.cispeps.push_back(cispep);
    // Gemmi❗✔️:     }
    // Gemmi❗✔️:   }
    // Gemmi❗✔️: }
    // Behavior review: each source record is dispatched by its first byte in
    // input order; skipped short rows do not increment counters, and the
    // three output families preserve their independent source vector order.
    // Complexity review: one pass over conn_records with fixed-field parsing;
    // Link ambiguity lookups and SSBOND address completion use their anchored
    // source-shaped scans and no complete hierarchy clone is made.
    let mut disulf_count = 0_i32;
    let mut link_counters = PdbLinkNameCounters::default();
    let mut connections = Vec::new();
    let mut cispeps = Vec::new();
    for record in records {
        let Some(kind) = record.first().copied() else {
            return Err(PdbConnectionFinalizationError::SsbondRecord(
                PdbSsbondRecordError::EmptyRecord,
            ));
        };
        match kind {
            b'S' | b's' => {
                if let Some(connection) =
                    decode_pdb_ssbond_record(record, &mut disulf_count, first_model, crystal)?
                {
                    connections.push(connection);
                }
            }
            b'L' | b'l' => {
                if let Some(connection) =
                    decode_pdb_link_record(record, first_model, &mut link_counters)
                        .map_err(PdbConnectionFinalizationError::LinkRecord)?
                {
                    connections.push(connection);
                }
            }
            b'C' | b'c' => {
                if let Some(cispep) = decode_pdb_cispep_record(record, models)
                    .map_err(PdbConnectionFinalizationError::CisPepRecord)?
                {
                    cispeps.push(cispep);
                }
            }
            _ => {}
        }
    }
    Ok((connections, cispeps))
}

fn pdb_finalization_text(
    field: &'static str,
    bytes: &[u8],
) -> Result<String, PdbStructureFinalizationError> {
    std::str::from_utf8(bytes)
        .map(str::to_owned)
        .map_err(|error| PdbStructureFinalizationError::Utf8 {
            field,
            valid_up_to: error.valid_up_to(),
        })
}

fn pdb_segment_source_id(segment: &str) -> Option<[u8; 4]> {
    if segment.is_empty() {
        return None;
    }
    let mut bytes = [b' '; 4];
    bytes[..segment.len()].copy_from_slice(segment.as_bytes());
    Some(bytes)
}

fn materialize_pdb_bio_structure(
    reader: PdbReaderState,
    transition: PdbModelTransition,
    connections: Vec<BioConnection>,
    cispeps: Vec<BioCisPep>,
) -> Result<BioStructureData, PdbStructureFinalizationError> {
    let mut models = Vec::with_capacity(transition.grouping.models.len());
    let mut chains = Vec::with_capacity(transition.grouping.chain_row_count);
    let mut residues = Vec::with_capacity(transition.grouping.residue_row_count);
    let atom_count = transition
        .grouping
        .models
        .iter()
        .flat_map(|model| &model.chains)
        .flat_map(|chain| &chain.residues)
        .map(|residue| residue.atoms.len())
        .sum();
    let mut atoms = Vec::with_capacity(atom_count);
    let mut positions = Vec::with_capacity(atom_count);

    for model in &transition.grouping.models {
        if model.row_id.index() != models.len() {
            return Err(PdbStructureFinalizationError::StagingRowOrder {
                table: "models",
                expected: models.len(),
                actual: model.row_id.index(),
            });
        }
        let chain_start = chains.len();
        for chain in &model.chains {
            if chain.row_id.index() != chains.len() {
                return Err(PdbStructureFinalizationError::StagingRowOrder {
                    table: "chains",
                    expected: chains.len(),
                    actual: chain.row_id.index(),
                });
            }
            let residue_start = residues.len();
            let entity_id = reader
                .entities
                .entities
                .iter()
                .position(|entity| entity.source_name == chain.source_id.as_bytes())
                .map(|index| {
                    u32::try_from(index)
                        .map(BioEntityId::new)
                        .map_err(|_| BioStructureError::RowIndexTooLarge { value: index })
                })
                .transpose()?;

            for residue in &chain.residues {
                if residue.row_id.index() != residues.len() {
                    return Err(PdbStructureFinalizationError::StagingRowOrder {
                        table: "residues",
                        expected: residues.len(),
                        actual: residue.row_id.index(),
                    });
                }
                let atom_start = atoms.len();
                for atom in &residue.atoms {
                    atoms.push(
                        BioAtomRow::new(
                            residue.row_id,
                            atom.name,
                            atom.element,
                            atom.isotope_mass_number,
                            atom.altloc,
                            atom.formal_charge,
                            BioCalcFlag::NotSet,
                            atom.occupancy,
                            atom.b_iso,
                            atom.anisou,
                            -1,
                            0.0,
                            AtomSourceIds::new(Some(atom.serial)),
                        )
                        .with_pdb_coordinate_text(atom.coordinate_text),
                    );
                    positions.push(atom.position);
                }
                let atom_span = BioRowSpan::from_usize(atom_start, atoms.len() - atom_start)?;
                let sequence_id = residue
                    .address
                    .sequence_number()
                    .map(|number| PdbSeqId::new(number, residue.address.insertion_code()));
                let segment_id = pdb_segment_source_id(residue.address.segment());
                let residue_source = ResidueSourceIds::new(
                    sequence_id,
                    None,
                    segment_id,
                    Some(residue.subchain_id.clone()),
                    None,
                )
                .ok_or(
                    PdbStructureFinalizationError::SourceIdNotRepresentable {
                        field: "PDB residue segment",
                    },
                )?;
                let residue_name =
                    pdb_finalization_text("PDB residue name", residue.address.name().as_bytes())?;
                let residue_info_kind: ResidueInfoKind = find_residue_info(&residue_name).kind;
                residues.push(BioResidueRow::new(
                    chain.row_id,
                    atom_span,
                    residue.address.name(),
                    residue_info_kind,
                    residue.entity_kind,
                    None,
                    Some(residue.het_flag),
                    residue_source,
                    BioSiftsUnpResidue::default(),
                ));
            }
            let residue_span =
                BioRowSpan::from_usize(residue_start, residues.len() - residue_start)?;
            chains.push(BioChainRow::new(
                model.row_id,
                entity_id,
                residue_span,
                ChainKind::Unknown,
                ChainSourceIds::new(Some(chain.source_id), None),
            ));
        }
        let chain_span = BioRowSpan::from_usize(chain_start, chains.len() - chain_start)?;
        models.push(BioModelRow::new(chain_span, model.source_model_number));
    }

    let mut entities = Vec::with_capacity(reader.entities.entities.len());
    for entity in reader.entities.entities {
        let source_entity_id =
            pdb_finalization_text("PDB entity source name", &entity.source_name)?;
        let full_sequence = entity
            .full_sequence
            .iter()
            .map(|row| pdb_finalization_text("PDB entity full sequence", row))
            .collect::<Result<Vec<_>, _>>()?;
        let dbrefs = entity
            .dbrefs
            .into_iter()
            .map(|dbref| {
                Ok(BioEntityDbRef {
                    db_name: pdb_finalization_text("PDB DBREF database", &dbref.db_name)?,
                    accession_code: pdb_finalization_text(
                        "PDB DBREF accession",
                        &dbref.accession_code,
                    )?,
                    id_code: pdb_finalization_text("PDB DBREF id", &dbref.id_code)?,
                    isoform: pdb_finalization_text("PDB DBREF isoform", &dbref.isoform)?,
                    seq_begin: dbref.seq_begin,
                    seq_end: dbref.seq_end,
                    db_begin: dbref.db_begin,
                    db_end: dbref.db_end,
                    label_seq_begin: dbref.label_seq_begin,
                    label_seq_end: dbref.label_seq_end,
                })
            })
            .collect::<Result<Vec<_>, PdbStructureFinalizationError>>()?;
        entities.push(BioEntityRow::new(
            entity.entity_kind,
            cosmolkit_bio::PolymerKind::Unknown,
            false,
            full_sequence,
            dbrefs,
            Vec::new(),
            entity.subchains,
            EntitySourceIds::new(source_entity_id),
        ));
    }

    let parts = BioStructureParts {
        input_format: BioCoordinateFormat::Pdb,
        models,
        chains,
        residues,
        atoms,
        entities,
        connections,
        cispeps,
        mod_residues: reader.mod_residues,
        helices: reader.helices,
        sheets: reader.sheets,
        metadata: reader.metadata,
        source_state: reader.source_state,
        coordinates: BioCoordinateBlock::new(positions),
        crystal: Some(reader.crystal_info),
        ncs_operators: reader.ncs_operators,
        assemblies: reader.assemblies,
    };
    BioStructureData::from_parts(parts).map_err(PdbStructureFinalizationError::Structure)
}

fn finalize_pdb_structure(
    mut reader: PdbReaderState,
    mut transition: PdbModelTransition,
    connection_records: &[Vec<u8>],
    skip_remarks: bool,
) -> Result<BioStructureData, PdbStructureFinalizationError> {
    // Gemmi source: third_party/gemmi/src/pdb.cpp::populate_structure_from_pdb_stream.
    // Gemmi❗✔️:   if (st.models.empty())
    // Gemmi❗✔️:     st.models.emplace_back(1);
    // Gemmi❗✔️:
    // Gemmi❗✔️:   if (st.ter_status == 'e')
    // Gemmi❗✔️:     remove_entity_types(st);
    // Gemmi❗✔️:
    // Gemmi❗✔️:   assign_subchains(st, /*force=*/false, /*fail_if_unknown=*/false);
    // Gemmi❗✔️:   for (Chain& ch : st.models[0].chains)
    // Gemmi❗✔️:     if (Entity* entity = st.get_entity(ch.name))
    // Gemmi❗✔️:       if (auto polymer = ch.get_polymer())
    // Gemmi❗✔️:         entity->subchains.emplace_back(polymer.subchain_id());
    // Gemmi❗✔️:
    // Gemmi❗✔️:   st.setup_cell_images();
    // Gemmi❗✔️:   process_conn(st, conn_records);
    // Gemmi❗✔️:
    // Gemmi❗✔️:   for (std::string& name : st.meta.authors)
    // Gemmi❗✔️:     change_author_name_format_to_mmcif(name);
    // Gemmi❗✔️:
    // Gemmi❗✔️:   if (!options.skip_remarks)
    // Gemmi❗✔️:     read_metadata_from_remarks(st);
    // Gemmi❗✔️:
    // Gemmi❗✔️:   restore_full_ccd_codes(st);
    // Behavior review: every source finalizer runs once in the pinned order.
    // Default model creation precedes TER rollback/subchains, cell images are
    // ready before SSBOND nearest-image completion, and CCD expansion occurs
    // only after all address-bearing records exist. The result is returned only
    // after the canonical BioStructureData validator accepts every row/span/map.
    // Complexity review: source helper traversals are retained in order; row
    // materialization is linear in the hierarchy and moves owned metadata once.
    transition.finish_eof()?;
    transition.finalize_subchains_and_entities(&reader.source_state, &mut reader.entities)?;
    setup_cell_images(&mut reader.crystal_info, &reader.ncs_operators);
    let first_model = transition
        .grouping
        .models
        .first()
        .ok_or(PdbStructureFinalizationError::NoModelAfterDefault)?;
    let (mut connections, mut cispeps) = process_pdb_connection_records(
        connection_records,
        first_model,
        &transition.grouping.models,
        &reader.crystal_info,
    )
    .map_err(PdbStructureFinalizationError::Connection)?;
    reader.finalize_author_names();
    if !skip_remarks {
        reader
            .read_metadata_from_remarks()
            .map_err(PdbStructureFinalizationError::Remarks)?;
    }
    reader
        .restore_full_ccd_codes(&mut transition.grouping, &mut connections, &mut cispeps)
        .map_err(PdbStructureFinalizationError::Ccd)?;
    materialize_pdb_bio_structure(reader, transition, connections, cispeps)
}

#[derive(Debug)]
enum PdbAtomFieldError {
    TooShort { len: usize },
    SerialOutsideSourceDefinedRange,
    ChainNameNotRepresentable { bytes: [u8; 2] },
    AtomNameNotRepresentable { bytes: [u8; 4] },
    ResidueKey(PdbResidueKeyError),
    Charge(PdbChargeError),
}

fn decode_pdb_atom_fields(
    line: &[u8],
    source_line_buffer: &[u8; 122],
) -> Result<PdbAtomFields, PdbAtomFieldError> {
    // Gemmi source: third_party/gemmi/src/pdb.cpp,
    // populate_structure_from_pdb_stream.
    // Gemmi❗✔️:       if (len < 55)
    // Gemmi❗✔️:         wrong("The line is too short to be correct:\n" + std::string(line));
    // The caller maps this typed boundary to Gemmi's line-numbered diagnostic;
    // no atom field is read before the same minimum-length gate.
    let len = line.len();
    if len < 55 {
        return Err(PdbAtomFieldError::TooShort { len });
    }

    // Gemmi❗✔️:       Atom atom;
    // Gemmi source: third_party/gemmi/include/gemmi/model.hpp, Atom.
    // Gemmi❗✔️:   SMat33<float> aniso = {0, 0, 0, 0, 0, 0};
    // Gemmi❗✔️:       atom.serial = read_serial(line+6);
    // Gemmi❗✔️:       atom.name = read_string(line+12, 4);
    // Gemmi❗✔️:       atom.altloc = read_altloc(line[16]);
    // Gemmi❗✔️:       atom.pos.x = read_double(line+30, 8);
    // Gemmi❗✔️:       atom.pos.y = read_double(line+38, 8);
    // Gemmi❗✔️:       atom.pos.z = read_double(line+46, 8);
    // Gemmi❗✔️:       if (len > 58)
    // Gemmi❗✔️:         atom.occ = (float) read_double(line+54, 6);
    // Gemmi❗✔️:       if (len > 64)
    // Gemmi❗✔️:         atom.b_iso = (float) read_double(line+60, 6);
    // Gemmi❗✔️:       if (len > 76 && (std::isalpha(line[76]) || std::isalpha(line[77])))
    // Gemmi❗✔️:         atom.element = Element(line + 76);
    // Gemmi❗✔️:       else
    // Gemmi❗✔️:         atom.element = infer_element_from_padded_name(line+12);
    // Gemmi❗✔️:       atom.charge = (len > 78 ? read_charge(line[78], line[79]) : 0);
    // The canonical AtomName deliberately retains the four raw columns under
    // PUB-bio_io; Gemmi's internal Atom.name is the separate read_string-trimmed
    // spelling. Element inference consumes the exact raw columns and carries
    // the approved H+Some(2) representation for Gemmi El::D. Fixed-width
    // read_double receives the reusable NUL-terminated source buffer, and
    // occupancy/B values are rounded through the source float storage type.
    let name_bytes: [u8; 4] = source_line_buffer[12..16]
        .try_into()
        .expect("fixed PDB line buffer contains the atom-name field");
    let name = raw_atom_name(&name_bytes)
        .ok_or(PdbAtomFieldError::AtomNameNotRepresentable { bytes: name_bytes })?;
    let serial_field: &[u8; 5] = source_line_buffer[6..11]
        .try_into()
        .expect("fixed PDB line buffer contains the atom serial field");
    let serial = read_serial(serial_field)
        .map(PdbAtomSerial::new)
        .ok_or(PdbAtomFieldError::SerialOutsideSourceDefinedRange)?;

    let position = [
        read_double(&source_line_buffer[30..38]),
        read_double(&source_line_buffer[38..46]),
        read_double(&source_line_buffer[46..54]),
    ];
    let occupancy = if len > 58 {
        f64::from(read_double(&source_line_buffer[54..60]) as f32)
    } else {
        1.0
    };
    let b_iso = if len > 64 {
        f64::from(read_double(&source_line_buffer[60..66]) as f32)
    } else {
        20.0
    };
    let (element, isotope_mass_number) = resolve_pdb_atom_element(&name_bytes, line);
    let formal_charge = if len > 78 {
        read_charge(source_line_buffer[78], source_line_buffer[79])
            .map_err(PdbAtomFieldError::Charge)?
    } else {
        0
    };

    // Behavior review: the decoder preserves source assignment order and
    // initializes all six ANISOU components to the source Atom default zero;
    // a following ANISOU record updates the current row through the dedicated
    // PdbModelTransition helper. It also preserves both
    // strict byte-count gates (`> 58`, `> 64`, `> 76`, `> 78`), Gemmi's
    // occupancy/B float rounding, exact PDB source serial, raw four-column
    // AtomName mapping, explicit-element precedence, isotope state, and the
    // source's typed charge error. Only source-undefined signed overflow and
    // the canonical AtomName ASCII representation boundary return typed
    // errors; neither is converted into a fabricated field value.
    // Complexity review: fixed-offset indexing and fixed-size conversions are
    // O(1), with no per-atom allocation on success; source helper scans are
    // bounded by their declared field widths.
    Ok(PdbAtomFields {
        coordinate_text: source_line_buffer[30..54]
            .try_into()
            .expect("fixed source field"),
        serial,
        name,
        altloc: read_altloc(source_line_buffer[16]),
        position,
        occupancy,
        b_iso,
        anisou: [0.0; 6],
        element,
        isotope_mass_number,
        formal_charge,
    })
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum PdbInputFormatError {
    Cif { source: String },
    Mmjson { source: String },
}

impl std::fmt::Display for PdbInputFormatError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Cif { source } => write!(
                formatter,
                "Incorrect file format (perhaps it is cif not pdb?): {source}"
            ),
            Self::Mmjson { source } => write!(
                formatter,
                "Incorrect file format (perhaps it is mmJSON not pdb?): {source}"
            ),
        }
    }
}

impl std::error::Error for PdbInputFormatError {}

fn gemmi_pdb_path_basename(path: &str) -> String {
    // Gemmi source: third_party/gemmi/include/gemmi/fileutil.hpp,
    // path_basename.
    // Gemmi❗✔️: inline std::string path_basename(const std::string& path,
    // Gemmi❗✔️:                                  std::initializer_list<const char*> exts) {
    // Gemmi❗✔️:   size_t pos = path.find_last_of("\\/");
    // Gemmi❗✔️:   std::string basename = pos == std::string::npos ? path : path.substr(pos + 1);
    // Gemmi❗✔️:   for (const char* ext : exts) {
    // Gemmi❗✔️:     size_t len = std::strlen(ext);
    // Gemmi❗✔️:     if (basename.size() > len &&
    // Gemmi❗✔️:         basename.compare(basename.length() - len, len, ext, len) == 0)
    // Gemmi❗✔️:       basename.resize(basename.length() - len);
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   return basename;
    // Gemmi❗✔️: }
    // Behavior review: both path separators are recognized; the ordered,
    // case-sensitive suffix tests require a nonempty stem, exactly as source.
    // Complexity review: one reverse separator scan, two suffix checks, and
    // one owned basename allocation; no full-path component normalization.
    let basename = path.rsplit(['\\', '/']).next().unwrap_or(path);
    let mut name = basename.to_owned();
    for extension in [".gz", ".pdb"] {
        if name.len() > extension.len() && name.ends_with(extension) {
            name.truncate(name.len() - extension.len());
        }
    }
    name
}

fn gemmi_record_type4(line: &[u8], offset: usize, record: [u8; 4]) -> bool {
    // Gemmi source: third_party/gemmi/include/gemmi/pdb.hpp and util.hpp.
    // Gemmi❗✔️: inline bool is_record_type4(const char* s, const char* record) {
    // Gemmi❗✔️:   return ialpha4_id(s) == ialpha4_id(record);
    // Gemmi❗✔️: }
    // Gemmi❗✔️: constexpr int ialpha4_id(const char* s) {
    // Gemmi❗✔️:   return (s[0] << 24 | s[1] << 16 | s[2] << 8 | s[3]) & ~0x20202020;
    // Gemmi❗✔️: }
    // The modeled ASCII record prefixes avoid the source's implementation-
    // defined/undefined signed-char shifts for high-bit input bytes.
    let actual = [
        line.get(offset).copied().unwrap_or(0),
        line.get(offset + 1).copied().unwrap_or(0),
        line.get(offset + 2).copied().unwrap_or(0),
        line.get(offset + 3).copied().unwrap_or(0),
    ];
    gemmi_ialpha4_id(&actual) == gemmi_ialpha4_id(&record)
}

fn gemmi_record_type3(line: &[u8], offset: usize, record: [u8; 4]) -> bool {
    // Gemmi source: third_party/gemmi/include/gemmi/pdb.hpp.
    // Gemmi❗✔️: inline bool is_record_type3(const char* s, const char* record) {
    // Gemmi❗✔️:   return (ialpha4_id(s) & ~0xf) == ialpha4_id(record);
    // Gemmi❗✔️: }
    // Behavior review: matching clears only the low nibble of the fourth
    // source byte, preserving Gemmi's accepted TER padding/control suffixes
    // without accepting other fourth-byte values. Missing Rust slice bytes
    // read as the same NUL-terminated C buffer initialized by the PDB reader.
    // Complexity review: four bounded indexed reads and two integer masks.
    let actual = [
        line.get(offset).copied().unwrap_or(0),
        line.get(offset + 1).copied().unwrap_or(0),
        line.get(offset + 2).copied().unwrap_or(0),
        line.get(offset + 3).copied().unwrap_or(0),
    ];
    gemmi_ialpha4_id(&actual) & !0x0000_000f == gemmi_ialpha4_id(&record)
}

fn gemmi_ialpha4_id(source: &[u8; 4]) -> u32 {
    // Gemmi source: third_party/gemmi/include/gemmi/util.hpp.
    // Gemmi❗✔️: constexpr int ialpha4_id(const char* s) {
    // Gemmi❗✔️:   return (s[0] << 24 | s[1] << 16 | s[2] << 8 | s[3]) & ~0x20202020;
    // Gemmi❗✔️: }
    // Behavior review: source PDB and residue-name inputs are ASCII bytes;
    // widening each byte before shifting matches the defined source domain
    // and avoids relying on signed-char overflow for high-bit input.
    // Complexity review: four fixed byte reads and constant-time bitwise ops.
    ((u32::from(source[0]) << 24)
        | (u32::from(source[1]) << 16)
        | (u32::from(source[2]) << 8)
        | u32::from(source[3]))
        & !0x2020_2020
}

pub(super) fn gemmi_pdb_residue_is_water(name: ResidueName) -> bool {
    // Gemmi source: third_party/gemmi/include/gemmi/model.hpp, Residue::is_water.
    // Gemmi❗✔️: bool is_water() const {
    // Gemmi❗✔️:   if (name.length() != 3)
    // Gemmi❗✔️:     return false;
    // Gemmi❗✔️:   int id = ialpha4_id(name.c_str());
    // Gemmi❗✔️:   return id == ialpha4_id("HOH") || id == ialpha4_id("DOD") ||
    // Gemmi❗✔️:          id == ialpha4_id("WAT") || id == ialpha4_id("H2O");
    // Gemmi❗✔️: }
    // Behavior review: the approved three-byte ASCII residue-name value is
    // NUL-padded exactly as its source string; the same four source aliases
    // and case-insensitive identifier are tested without chemical inference.
    // Complexity review: at most three source bytes and four constant IDs.
    let bytes = name.as_bytes();
    if bytes.len() != 3 {
        return false;
    }
    let mut padded = [0; 4];
    padded[..3].copy_from_slice(bytes);
    let id = gemmi_ialpha4_id(&padded);
    [*b"HOH\0", *b"DOD\0", *b"WAT\0", *b"H2O\0"]
        .into_iter()
        .any(|water_name| id == gemmi_ialpha4_id(&water_name))
}

fn gemmi_ialpha3_id(line: &[u8], offset: usize) -> u32 {
    // Gemmi source: third_party/gemmi/include/gemmi/util.hpp.
    // Gemmi❗✔️: constexpr int ialpha3_id(const char* s) {
    // Gemmi❗✔️:   return (s[0] << 16 | s[1] << 8 | s[2]) & ~0x202020;
    // Gemmi❗✔️: }
    let first = line.get(offset).copied().unwrap_or(0) as u32;
    let second = line.get(offset + 1).copied().unwrap_or(0) as u32;
    let third = line.get(offset + 2).copied().unwrap_or(0) as u32;
    ((first << 16) | (second << 8) | third) & !0x0020_2020
}

pub(crate) fn gemmi_is_space(byte: u8) -> bool {
    // Gemmi source: third_party/gemmi/include/gemmi/atox.hpp
    // Gemmi❗✔️: inline bool is_space(char c) {
    // Gemmi❗✔️:   static const std::uint8_t table[256] = { // 1 for 9-13 and 32
    // Gemmi❗✔️:     0,0,0,0,0,0,0,0, 0,1,1,1,1,1,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
    // Gemmi❗✔️:     1,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
    // Gemmi❗✔️:     0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
    // Gemmi❗✔️:     0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
    // Gemmi❗✔️:     0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
    // Gemmi❗✔️:     0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
    // Gemmi❗✔️:     0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
    // Gemmi❗✔️:     0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,
    // Gemmi❗✔️:   };
    // Gemmi❗✔️:   return table[(std::uint8_t)c] != 0;
    // Gemmi❗✔️: }
    matches!(byte, b'\t'..=b'\r' | b' ')
}

fn gemmi_is_digit(byte: u8) -> bool {
    // Gemmi source: third_party/gemmi/include/gemmi/atox.hpp
    // Gemmi❗✔️: inline bool is_digit(char c) {
    // Gemmi❗✔️:   return c >= '0' && c <= '9';
    // Gemmi❗✔️: }
    (b'0'..=b'9').contains(&byte)
}

fn gemmi_pdb_skip_blank(bytes: &[u8], start: usize) -> usize {
    // Gemmi source: third_party/gemmi/include/gemmi/atox.hpp.
    // Gemmi❗✔️: inline const char* skip_blank(const char* p) {
    // Gemmi❗✔️:   if (p)
    // Gemmi❗✔️:     while (is_blank(*p))
    // Gemmi❗✔️:       ++p;
    // Gemmi❗✔️:   return p;
    // Gemmi❗✔️: }
    // Behavior review: `is_blank` is space or tab only; the input slice ends
    // at the source C-string terminator, so the cursor never consumes past it.
    // Complexity review: one forward pass, O(n) time and O(1) auxiliary state.
    let mut index = start.min(bytes.len());
    while index < bytes.len() && matches!(bytes[index], b' ' | b'\t') {
        index += 1;
    }
    index
}

fn gemmi_pdb_rtrim_cstr_end(bytes: &[u8], start: usize, end: Option<usize>) -> usize {
    // Gemmi source: third_party/gemmi/include/gemmi/util.hpp.
    // Gemmi❗✔️: inline const char* rtrim_cstr(const char* start, const char* end=nullptr) {
    // Gemmi❗✔️:   if (!start)
    // Gemmi❗✔️:     return nullptr;
    // Gemmi❗✔️:   if (!end) {
    // Gemmi❗✔️:     end = start;
    // Gemmi❗✔️:     while (*end != '\0')
    // Gemmi❗✔️:       ++end;
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   while (end > start && std::isspace(end[-1]))
    // Gemmi❗✔️:     --end;
    // Gemmi❗✔️:   return end;
    // Gemmi❗✔️: }
    // Behavior review: None scans to the already-isolated C-string end; Some
    // preserves the explicit colon boundary. Trailing bytes use Gemmi's
    // pinned C-locale whitespace table, distinct from skip_blank's isblank.
    // Complexity review: one backwards scan, O(n) time and O(1) storage.
    let mut end = end.unwrap_or(bytes.len()).min(bytes.len());
    while end > start && gemmi_is_space(bytes[end - 1]) {
        end -= 1;
    }
    end
}

fn gemmi_pdb_is_double(value: &[u8]) -> Result<bool, PdbRemark200Error> {
    // Gemmi source: third_party/gemmi/src/pdb.cpp.
    // Gemmi❗✔️: bool is_double(const char* p) {
    // Gemmi❗✔️:   while (is_space(*p)) ++p;
    // Gemmi❗✔️:   if (*p == '-' || *p == '+') ++p;
    // Gemmi❗✔️:   while (is_digit(*p)) ++p;
    // Gemmi❗✔️:   if (*p == '.') {
    // Gemmi❗✔️:     ++p;
    // Gemmi❗✔️:     while (is_digit(*++p)) ++p;
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   while (is_space(*p)) ++p;
    // Gemmi❗✔️:   return *p == '\0';
    // Gemmi❗✔️: }
    // Behavior review: preserve the source's pre-increment in the fractional
    // digit condition and its additional loop-body increment. If that sequence
    // would increment beyond the C terminator, return an explicit boundary
    // error rather than execute undefined pointer access; no parity is claimed
    // for that source-undefined pH screening input.
    // Complexity review: one cursor pass, O(n) time and O(1) auxiliary state.
    let c_end = value
        .iter()
        .position(|byte| *byte == 0)
        .unwrap_or(value.len());
    let value = &value[..c_end];
    let mut index = 0;

    while index < value.len() && gemmi_is_space(value[index]) {
        index += 1;
    }
    if value
        .get(index)
        .is_some_and(|byte| matches!(byte, b'-' | b'+'))
    {
        index += 1;
    }
    while value.get(index).is_some_and(|byte| gemmi_is_digit(*byte)) {
        index += 1;
    }
    if value.get(index) == Some(&b'.') {
        index += 1;
        loop {
            if index == value.len() {
                return Err(PdbRemark200Error::UndefinedPhNumericScreening);
            }
            index += 1;
            if !gemmi_is_digit(value.get(index).copied().unwrap_or(0)) {
                break;
            }
            index += 1;
        }
    }
    while index < value.len() && gemmi_is_space(value[index]) {
        index += 1;
    }
    Ok(index == value.len())
}

fn gemmi_alpha_up(byte: u8) -> u8 {
    // Gemmi source: third_party/gemmi/include/gemmi/util.hpp
    // Gemmi❗✔️: inline char alpha_up(char c) { return c & ~0x20; }
    // Behavior review: the source is used only with the four ASCII atom-name
    // columns required by `AtomName`; the bit mask therefore uppercases ASCII
    // lowercase letters and leaves the other modeled bytes unchanged.
    // Complexity review: one constant-time byte mask, with no allocation.
    byte & !0x20
}

fn pdb_date_format_to_iso(date: &[u8]) -> String {
    // Gemmi source: third_party/gemmi/src/pdb.cpp and its included helpers.
    // Gemmi❗✔️: inline bool is_digit(char c) {
    // Gemmi❗✔️:   return c >= '0' && c <= '9';
    // Gemmi❗✔️: }
    // Gemmi❗✔️: inline char alpha_up(char c) { return c & ~0x20; }
    // Gemmi❗✔️: // "28-MAR-07" -> "2007-03-28"
    // Gemmi❗✔️: // (we also accept less standard format "28-Mar-2007" as used by BUSTER)
    // Gemmi❗✔️: // We do not check if the date is correct.
    // Gemmi❗✔️: // The returned value is one of:
    // Gemmi❗✔️: //   DDDD-DD-DD - possibly correct date,
    // Gemmi❗✔️: //   DDDD-xx-DD - unrecognized month,
    // Gemmi❗✔️: //   empty string - the digits were not there.
    // Gemmi❗✔️: std::string pdb_date_format_to_iso(const std::string& date) {
    // Gemmi❗✔️:   const char months[] = "JAN01FEB02MAR03APR04MAY05JUN06"
    // Gemmi❗✔️:                         "JUL07AUG08SEP09OCT10NOV11DEC122222";
    // Gemmi❗✔️:   if (date.size() < 9 || !is_digit(date[0]) || !is_digit(date[1]) ||
    // Gemmi❗✔️:                          !is_digit(date[7]) || !is_digit(date[8]))
    // Gemmi❗✔️:     return std::string();
    // Gemmi❗✔️:   std::string iso = "xxxx-xx-xx";
    // Gemmi❗✔️:   if (date.size() >= 11 && is_digit(date[9]) && is_digit(date[10])) {
    // Gemmi❗✔️:     std::memcpy(&iso[0], &date[7], 4);
    // Gemmi❗✔️:   } else {
    // Gemmi❗✔️:     std::memcpy(&iso[0], (date[7] > '6' ? "19" : "20"), 2);
    // Gemmi❗✔️:     std::memcpy(&iso[2], &date[7], 2);
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   char month[4] = {alpha_up(date[3]), alpha_up(date[4]), alpha_up(date[5]), '\0'};
    // Gemmi❗✔️:   if (const char* m = std::strstr(months, month))
    // Gemmi❗✔️:     std::memcpy(&iso[5], m + 3, 2);
    // Gemmi❗✔️:   std::memcpy(&iso[8], &date[0], 2);
    // Gemmi❗✔️:   return iso;
    // Gemmi❗✔️: }
    //
    // Behavior review: the fixed character positions, short-year pivot, and
    // raw substring month lookup are preserved; the source explicitly skips
    // calendar validation. C-string needle termination is also preserved.
    // Complexity review: the month table has fixed size, so its single window
    // scan is bounded; the result is built once with constant-size storage.
    if date.len() < 9
        || !gemmi_is_digit(date[0])
        || !gemmi_is_digit(date[1])
        || !gemmi_is_digit(date[7])
        || !gemmi_is_digit(date[8])
    {
        return String::new();
    }

    const MONTHS: &[u8] = b"JAN01FEB02MAR03APR04MAY05JUN06JUL07AUG08SEP09OCT10NOV11DEC122222";
    let mut iso = *b"xxxx-xx-xx";
    if date.len() >= 11 && gemmi_is_digit(date[9]) && gemmi_is_digit(date[10]) {
        iso[..4].copy_from_slice(&date[7..11]);
    } else {
        iso[..2].copy_from_slice(if date[7] > b'6' { b"19" } else { b"20" });
        iso[2..4].copy_from_slice(&date[7..9]);
    }

    let month = [
        gemmi_alpha_up(date[3]),
        gemmi_alpha_up(date[4]),
        gemmi_alpha_up(date[5]),
        0,
    ];
    let needle_len = month.iter().position(|byte| *byte == 0).unwrap_or(3);
    let month_start = if needle_len == 0 {
        Some(0)
    } else {
        MONTHS
            .windows(needle_len)
            .position(|candidate| candidate == &month[..needle_len])
    };
    if let Some(month_start) = month_start {
        iso[5..7].copy_from_slice(&MONTHS[month_start + 3..month_start + 5]);
    }
    iso[8..10].copy_from_slice(&date[..2]);
    String::from_utf8(iso.to_vec()).expect("Gemmi PDB date output is ASCII")
}

fn change_author_name_format_to_mmcif(name: &mut String) {
    // Gemmi source: third_party/gemmi/src/pdb.cpp.
    // Gemmi❗✔️: // move initials after comma, as in mmCIF (A.-B.DOE -> DOE, A.-B.), see
    // Gemmi❗✔️: // https://www.wwpdb.org/documentation/file-format-content/format33/sect2.html#AUTHOR
    // Gemmi❗✔️: void change_author_name_format_to_mmcif(std::string& name) {
    // Gemmi❗✔️:   // If the AUTHOR record has comma followed by space we get leading space here
    // Gemmi❗✔️:   while (name[0] == ' ')
    // Gemmi❗✔️:     name.erase(name.begin());
    // Gemmi❗✔️:   size_t pos = 0;
    // Gemmi❗✔️:   // Initials may have multiple letters (e.g. JU. or PON.)
    // Gemmi❗✔️:   // but should not have space after dot.
    // Gemmi❗✔️:   for (size_t i = 1; i < pos+4 && i+1 < name.size(); ++i)
    // Gemmi❗✔️:     if (name[i] == '.' && name[i+1] != ' ')
    // Gemmi❗✔️:       pos = i+1;
    // Gemmi❗✔️:   if (pos > 0)
    // Gemmi❗✔️:     name = name.substr(pos) + ", " + name.substr(0, pos);
    // Gemmi❗✔️: }
    //
    // Behavior review: only ASCII space is stripped; the dot scan starts at
    // byte 1, its bound extends whenever another initial is found, and empty
    // names leave the loop without indexing or mutation.
    // Complexity review: repeated front removal intentionally retains the
    // source's worst-case quadratic shifting; the initials scan is linear and
    // the final reorder uses one bounded allocation.
    while name.as_bytes().first() == Some(&b' ') {
        name.remove(0);
    }

    let mut pos = 0;
    let mut i = 1;
    while i < pos + 4 && i + 1 < name.len() {
        let bytes = name.as_bytes();
        if bytes[i] == b'.' && bytes[i + 1] != b' ' {
            pos = i + 1;
        }
        i += 1;
    }
    if pos > 0 {
        let mut reordered = String::with_capacity(name.len() + 2);
        reordered.push_str(&name[pos..]);
        reordered.push_str(", ");
        reordered.push_str(&name[..pos]);
        *name = reordered;
    }
}

fn pdb_header_text(bytes: &[u8], field_offset: usize) -> Result<String, PdbHeaderError> {
    String::from_utf8(bytes.to_vec()).map_err(|error| PdbHeaderError::TextFieldNotUtf8 {
        field_offset,
        valid_up_to: error.utf8_error().valid_up_to(),
    })
}

fn gemmi_pdb_rtrim_text(bytes: &[u8]) -> &[u8] {
    // Gemmi source: third_party/gemmi/include/gemmi/util.hpp.
    // Gemmi❗✔️: inline std::string rtrim_str(const std::string& str) {
    // Gemmi❗✔️:   std::string::size_type last = str.find_last_not_of(" \r\n\t");
    // Gemmi❗✔️:   return str.substr(0, last == std::string::npos ? 0 : last + 1);
    // Gemmi❗✔️: }
    // Behavior review: trimming is limited to the exact four source bytes;
    // vertical tab and form feed are not treated as whitespace here.
    // Complexity review: one reverse scan and a borrowed slice; no extra
    // allocation beyond the later canonical String construction.
    let Some(end) = bytes
        .iter()
        .rposition(|byte| !matches!(byte, b' ' | b'\r' | b'\n' | b'\t'))
    else {
        return &bytes[..0];
    };
    &bytes[..end + 1]
}

fn gemmi_pdb_trim_text(bytes: &[u8]) -> &[u8] {
    // Gemmi source: third_party/gemmi/include/gemmi/util.hpp.
    // Gemmi❗✔️: inline std::string trim_str(const std::string& str) {
    // Gemmi❗✔️:   const std::string ws = " \r\n\t";
    // Gemmi❗✔️:   std::string::size_type first = str.find_first_not_of(ws);
    // Gemmi❗✔️:   if (first == std::string::npos)
    // Gemmi❗✔️:     return std::string{};
    // Gemmi❗✔️:   std::string::size_type last = str.find_last_not_of(ws);
    // Gemmi❗✔️:   return str.substr(first, last - first + 1);
    // Gemmi❗✔️: }
    // Behavior review: both ends use exactly space, CR, LF, and tab, with an
    // all-trim-byte field becoming empty; other byte values remain unchanged.
    // Complexity review: two bounded scans and one borrowed slice.
    let Some(start) = bytes
        .iter()
        .position(|byte| !matches!(byte, b' ' | b'\r' | b'\n' | b'\t'))
    else {
        return &bytes[..0];
    };
    let end = bytes
        .iter()
        .rposition(|byte| !matches!(byte, b' ' | b'\r' | b'\n' | b'\t'))
        .expect("a non-trim byte was found");
    &bytes[start..end + 1]
}

fn gemmi_pdb_author_fields(
    line: &[u8; 122],
    line_len: usize,
) -> Result<Vec<String>, PdbHeaderError> {
    // Gemmi source helper: third_party/gemmi/include/gemmi/atox.hpp.
    // Gemmi❗✔️: inline bool is_blank(char c) {
    // Gemmi❗✔️:   return c == ' ' || c == '\t';
    // Gemmi❗✔️: }
    // Gemmi❗✔️: inline const char* skip_blank(const char* p) {
    // Gemmi❗✔️:   if (p)
    // Gemmi❗✔️:     while (is_blank(*p))
    // Gemmi❗✔️:       ++p;
    // Gemmi❗✔️:   return p;
    // Gemmi❗✔️: }
    // Gemmi source helper: third_party/gemmi/include/gemmi/util.hpp.
    // Gemmi❗✔️: inline const char* rtrim_cstr(const char* start, const char* end=nullptr) {
    // Gemmi❗✔️:   if (!start)
    // Gemmi❗✔️:     return nullptr;
    // Gemmi❗✔️:   if (!end) {
    // Gemmi❗✔️:     end = start;
    // Gemmi❗✔️:     while (*end != '\0')
    // Gemmi❗✔️:       ++end;
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   while (end > start && std::isspace(end[-1]))
    // Gemmi❗✔️:     --end;
    // Gemmi❗✔️:   return end;
    // Gemmi❗✔️: }
    // Gemmi source helper: third_party/gemmi/include/gemmi/util.hpp.
    // Gemmi❗✔️: template<typename S>
    // Gemmi❗✔️: void split_str_into(const std::string& str, S sep,
    // Gemmi❗✔️:                     std::vector<std::string>& result) {
    // Gemmi❗✔️:   std::size_t start = 0, end;
    // Gemmi❗✔️:   while ((end = str.find(sep, start)) != std::string::npos) {
    // Gemmi❗✔️:     result.emplace_back(str, start, end - start);
    // Gemmi❗✔️:     start = end + impl::length(sep);
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   result.emplace_back(str, start);
    // Gemmi❗✔️: }
    // Behavior review: the leading scan removes only space/tab; the trailing
    // scan uses the pinned C-locale isspace set; comma splitting retains empty
    // leading, interior, and trailing components like split_str_into.
    // Complexity review: one linear leading/trailing scan, one linear split,
    // and one owned String per retained component; no repeated rescans.
    let mut start = 10;
    while matches!(line[start], b' ' | b'\t') {
        start += 1;
    }
    let mut end = line_len;
    while end > start && matches!(line[end - 1], b' ' | b'\t'..=b'\r') {
        end -= 1;
    }

    let value = &line[start..end];
    let mut fields = Vec::new();
    let mut field_start = 0;
    for (index, byte) in value.iter().enumerate() {
        if *byte == b',' {
            fields.push(pdb_header_text(
                &value[field_start..index],
                start + field_start,
            )?);
            field_start = index + 1;
        }
    }
    fields.push(pdb_header_text(&value[field_start..], start + field_start)?);
    Ok(fields)
}

fn gemmi_single_letter_element(byte: u8) -> (Element, Option<u16>) {
    // Gemmi source: third_party/gemmi/include/gemmi/elem.hpp
    // Gemmi❗✔️: inline El find_single_letter_element(char c) {
    // Gemmi❗✔️:   switch (c) {
    // Gemmi❗✔️:     case 'H': return El::H;
    // Gemmi❗✔️:     case 'B': return El::B;
    // Gemmi❗✔️:     case 'C': return El::C;
    // Gemmi❗✔️:     case 'N': return El::N;
    // Gemmi❗✔️:     case 'O': return El::O;
    // Gemmi❗✔️:     case 'F': return El::F;
    // Gemmi❗✔️:     case 'P': return El::P;
    // Gemmi❗✔️:     case 'S': return El::S;
    // Gemmi❗✔️:     case 'K': return El::K;
    // Gemmi❗✔️:     case 'V': return El::V;
    // Gemmi❗✔️:     case 'Y': return El::Y;
    // Gemmi❗✔️:     case 'I': return El::I;
    // Gemmi❗✔️:     case 'W': return El::W;
    // Gemmi❗✔️:     case 'U': return El::U;
    // Gemmi❗✔️:     case 'D': return El::D;
    // Gemmi❗✔️:     default: return El::X;
    // Gemmi❗✔️:   }
    // Gemmi❗✔️: }
    // Behavior review: Gemmi's closed one-letter dispatch is preserved;
    // source `El::D` is projected to the canonical BIO pair H+Some(2), while
    // source `El::X` becomes `Element::DUMMY` with unspecified isotope.
    // Complexity review: one byte match and constant-size value construction.
    match byte {
        b'H' => (Element::H, None),
        b'B' => (Element::B, None),
        b'C' => (Element::C, None),
        b'N' => (Element::N, None),
        b'O' => (Element::O, None),
        b'F' => (Element::F, None),
        b'P' => (Element::P, None),
        b'S' => (Element::S, None),
        b'K' => (Element::K, None),
        b'V' => (Element::V, None),
        b'Y' => (Element::Y, None),
        b'I' => (Element::I, None),
        b'W' => (Element::W, None),
        b'U' => (Element::U, None),
        b'D' => (Element::H, Some(2)),
        _ => (Element::DUMMY, None),
    }
}

pub(super) fn gemmi_find_element(field: [u8; 2]) -> (Element, Option<u16>) {
    // Gemmi source: third_party/gemmi/include/gemmi/elem.hpp
    // Gemmi❗✔️: inline elname_t& element_uppercase_name(El el) {
    // Gemmi❗✔️:   static constexpr elname_t names[] = {
    // Gemmi❗✔️:     "X",  "H",  "HE", "LI", "BE", "B",  "C",  "N",  "O", "F", "NE",
    // Gemmi❗✔️:     "NA", "MG", "AL", "SI", "P",  "S",  "CL", "AR",
    // Gemmi❗✔️:     "K",  "CA", "SC", "TI", "V",  "CR", "MN", "FE", "CO",
    // Gemmi❗✔️:     "NI", "CU", "ZN", "GA", "GE", "AS", "SE", "BR", "KR",
    // Gemmi❗✔️:     "RB", "SR", "Y",  "ZR", "NB", "MO", "TC", "RU", "RH",
    // Gemmi❗✔️:     "PD", "AG", "CD", "IN", "SN", "SB", "TE", "I", "XE",
    // Gemmi❗✔️:     "CS", "BA", "LA", "CE", "PR", "ND", "PM", "SM", "EU",
    // Gemmi❗✔️:     "GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU",
    // Gemmi❗✔️:     "HF", "TA", "W",  "RE", "OS", "IR", "PT", "AU", "HG",
    // Gemmi❗✔️:     "TL", "PB", "BI", "PO", "AT", "RN",
    // Gemmi❗✔️:     "FR", "RA", "AC", "TH", "PA", "U",  "NP", "PU", "AM",
    // Gemmi❗✔️:     "CM", "BK", "CF", "ES", "FM", "MD", "NO", "LR",
    // Gemmi❗✔️:     "RF", "DB", "SG", "BH", "HS", "MT", "DS", "RG", "CN",
    // Gemmi❗✔️:     "NH", "FL", "MC", "LV", "TS", "OG",
    // Gemmi❗✔️:     "D", "", ""
    // Gemmi❗✔️:   };
    // Gemmi❗✔️:   static_assert(sizeof(names) / sizeof(names[0]) == 122, "not 122");
    // Gemmi❗✔️:   return names[static_cast<int>(el)];
    // Gemmi❗✔️: }
    // Gemmi❗✔️: namespace impl {
    // Gemmi❗✔️: inline El find_single_letter_element(char c) {
    // Gemmi❗✔️:   switch (c) {
    // Gemmi❗✔️:     case 'H': return El::H;
    // Gemmi❗✔️:     case 'B': return El::B;
    // Gemmi❗✔️:     case 'C': return El::C;
    // Gemmi❗✔️:     case 'N': return El::N;
    // Gemmi❗✔️:     case 'O': return El::O;
    // Gemmi❗✔️:     case 'F': return El::F;
    // Gemmi❗✔️:     case 'P': return El::P;
    // Gemmi❗✔️:     case 'S': return El::S;
    // Gemmi❗✔️:     case 'K': return El::K;
    // Gemmi❗✔️:     case 'V': return El::V;
    // Gemmi❗✔️:     case 'Y': return El::Y;
    // Gemmi❗✔️:     case 'I': return El::I;
    // Gemmi❗✔️:     case 'W': return El::W;
    // Gemmi❗✔️:     case 'U': return El::U;
    // Gemmi❗✔️:     case 'D': return El::D;
    // Gemmi❗✔️:     default: return El::X;
    // Gemmi❗✔️:   }
    // Gemmi❗✔️: }
    // Gemmi❗✔️: } // namespace impl
    // Gemmi❗✔️: inline El find_element(const char* symbol) {
    // Gemmi❗✔️:   if (symbol == nullptr || symbol[0] == '\0')
    // Gemmi❗✔️:     return El::X;
    // Gemmi❗✔️:   char first = symbol[0] & ~0x20;  // lower -> upper, space -> NUL
    // Gemmi❗✔️:   char second = symbol[1] & ~0x20;
    // Gemmi❗✔️:   if (first == '\0')
    // Gemmi❗✔️:     return impl::find_single_letter_element(second);
    // Gemmi❗✔️:   // To handle symbol being "X\n" we have the condition below.
    // Gemmi❗✔️:   // In addition to \t, \v, \r and \n it catches also !"#$%&'()*+,- and
    // Gemmi❗✔️:   // some control characters - inconsistent but not necessarily bad.
    // Gemmi❗✔️:   if (second < 14)
    // Gemmi❗✔️:     return impl::find_single_letter_element(first);
    // Gemmi❗✔️:   elname_t* names = &element_uppercase_name(El::X);
    // Gemmi❗✔️:   for (int i = 0; i != 120; ++i)
    // Gemmi❗✔️:     if (names[i][0] == first && names[i][1] == second)
    // Gemmi❗✔️:       return static_cast<El>(i);
    // Gemmi❗✔️:   return El::X;
    // Gemmi❗✔️: }
    // Behavior review: Gemmi's table search is reproduced over the canonical
    // type symbols (with its source `X` spelling for `Element::DUMMY`); the
    // one-letter fallback, byte mask, and `< 14` second-byte branch remain
    // distinct. Its `D` enum result maps to H+Some(2), and all unknown rows to
    // DUMMY+None. Input bytes are the source's ASCII two-byte symbol domain.
    // Complexity review: the linear lookup checks 119 canonical entries plus
    // the source's D case, matching Gemmi's fixed 120-entry scan; no allocation.
    if field[0] == 0 {
        return (Element::DUMMY, None);
    }

    let first = gemmi_alpha_up(field[0]);
    let second = gemmi_alpha_up(field[1]);
    if first == 0 {
        return gemmi_single_letter_element(second);
    }
    if second < 14 {
        return gemmi_single_letter_element(first);
    }

    for element in Element::iter_with_dummy() {
        let symbol = if element == Element::DUMMY {
            "X"
        } else {
            element.symbol()
        };
        let symbol = symbol.as_bytes();
        let source_first = gemmi_alpha_up(symbol[0]);
        let source_second = symbol.get(1).copied().map(gemmi_alpha_up).unwrap_or(0);
        if first == source_first && second == source_second {
            return (element, None);
        }
    }

    if first == b'D' && second == 0 {
        return (Element::H, Some(2));
    }

    (Element::DUMMY, None)
}

fn infer_element_from_padded_name(name: &[u8; 4]) -> (Element, Option<u16>) {
    // Gemmi source: third_party/gemmi/src/pdb.cpp
    // Gemmi❗✔️: El infer_element_from_padded_name(const char* name) {
    // Gemmi❗✔️:   // Old versions of the PDB format had hydrogen names such as "1HB ".
    // Gemmi❗✔️:   // Some MD files use similar names for other elements ("1C4A" -> C).
    // Gemmi❗✔️:   if (name[0] == ' ' || is_digit(name[0]))
    // Gemmi❗✔️:     return impl::find_single_letter_element(name[1]);
    // Gemmi❗✔️:   // ... or it can be "C210"
    // Gemmi❗✔️:   if (is_digit(name[1]))
    // Gemmi❗✔️:     return impl::find_single_letter_element(name[0]);
    // Gemmi❗✔️:   if (name[3] != ' ') {
    // Gemmi❗✔️:     // Atom names HXXX are ambiguous, but Hg, He, Hf, Ho and Hs (almost)
    // Gemmi❗✔️:     // never have 4-character names, so H is assumed.
    // Gemmi❗✔️:     if (alpha_up(name[0]) == 'H')
    // Gemmi❗✔️:       return El::H;
    // Gemmi❗✔️:     // Similarly Deuterium (DXXX), but here alternatives are Dy, Db and Ds.
    // Gemmi❗✔️:     // Only Dysprosium is present in the PDB - in a single entry as of 2022.
    // Gemmi❗✔️:     if (alpha_up(name[0]) == 'D')
    // Gemmi❗✔️:       return El::D;
    // Gemmi❗✔️:     // Don't try harder for now. We don't recognize names such as CG11 as C
    // Gemmi❗✔️:     // (which we could; there is no Cg in the periodic table), but
    // Gemmi❗✔️:     // a name such as CL20 can be either Cl (in WGW) or C (in WQH) ¯\_(ツ)_/¯
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   return find_element(name);
    // Gemmi❗✔️: }
    // Behavior review: branch order is preserved over the exact four-byte
    // PDB source field. Name-derived H is H+None, Gemmi D is H+Some(2), and
    // the final lookup preserves the source's ambiguity decisions (including
    // CL20/CG11) rather than interpreting atom names heuristically.
    // Complexity review: two fixed-width digit checks and at most one fixed
    // 120-entry source-shaped symbol lookup; no allocation.
    if name[0] == b' ' || gemmi_is_digit(name[0]) {
        return gemmi_single_letter_element(name[1]);
    }
    if gemmi_is_digit(name[1]) {
        return gemmi_single_letter_element(name[0]);
    }
    if name[3] != b' ' {
        if gemmi_alpha_up(name[0]) == b'H' {
            return (Element::H, None);
        }
        if gemmi_alpha_up(name[0]) == b'D' {
            return (Element::H, Some(2));
        }
    }
    gemmi_find_element([name[0], name[1]])
}

fn resolve_pdb_atom_element(name: &[u8; 4], line: &[u8]) -> (Element, Option<u16>) {
    // Gemmi source: third_party/gemmi/src/pdb.cpp, ATOM/HETATM element selection
    // Gemmi❗✔️:       if (len > 76 && (std::isalpha(line[76]) || std::isalpha(line[77])))
    // Gemmi❗✔️:         atom.element = Element(line + 76);
    // Gemmi❗✔️:       else
    // Gemmi❗✔️:         atom.element = infer_element_from_padded_name(line+12);
    // Gemmi source: third_party/gemmi/include/gemmi/elem.hpp
    // Gemmi❗✔️:   explicit Element(const char* str) noexcept : elem(find_element(str)) {}
    // Behavior review: an explicitly recognized ASCII element field overrides
    // the atom-name inference, including conflicts. A missing byte at column
    // 78 is the source line buffer's NUL terminator. Non-ASCII `isalpha` input
    // is outside the modeled source-defined PDB text domain.
    // Complexity review: constant-time field selection followed by at most one
    // fixed 120-entry lookup; no allocation.
    let first = line.get(76).copied().unwrap_or(0);
    let second = line.get(77).copied().unwrap_or(0);
    if line.len() > 76 && (first.is_ascii_alphabetic() || second.is_ascii_alphabetic()) {
        gemmi_find_element([first, second])
    } else {
        infer_element_from_padded_name(name)
    }
}

pub(super) fn read_int(field: &[u8]) -> Option<i32> {
    // Gemmi source: third_party/gemmi/src/pdb.cpp
    // Gemmi❗✔️: int read_int(const char* p, int field_length) {
    // Gemmi❗✔️:   return string_to_int(p, false, field_length);
    // Gemmi❗✔️: }
    // Gemmi source: third_party/gemmi/include/gemmi/atox.hpp
    // Gemmi❗✔️: // no checking for overflow
    // Gemmi❗✔️: inline int string_to_int(const char* p, bool checked, size_t length=0) {
    // Gemmi❗✔️:   int mult = -1;
    // Gemmi❗✔️:   int n = 0;
    // Gemmi❗✔️:   size_t i = 0;
    // Gemmi❗✔️:   while ((length == 0 || i < length) && is_space(p[i]))
    // Gemmi❗✔️:     ++i;
    // Gemmi❗✔️:   if (p[i] == '-') {
    // Gemmi❗✔️:     mult = 1;
    // Gemmi❗✔️:     ++i;
    // Gemmi❗✔️:   } else if (p[i] == '+') {
    // Gemmi❗✔️:     ++i;
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   bool has_digits = false;
    // Gemmi❗✔️:   // use negative numbers because INT_MIN < -INT_MAX
    // Gemmi❗✔️:   for (; (length == 0 || i < length) && is_digit(p[i]); ++i) {
    // Gemmi❗✔️:     n = n * 10 - (p[i] - '0');
    // Gemmi❗✔️:     has_digits = true;
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   if (checked) {
    // Gemmi❗✔️:     while ((length == 0 || i < length) && is_space(p[i]))
    // Gemmi❗✔️:       ++i;
    // Gemmi❗✔️:     if (!has_digits || p[i] != '\0')
    // Gemmi❗✔️:       throw std::invalid_argument("not an integer: " +
    // Gemmi❗✔️:                            std::string(p, length ? length : i+1));
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   return mult * n;
    // Gemmi❗✔️: }
    // Behavior review: this PDB call fixes `checked` to false, so the parser
    // returns zero for no digits, accepts an optional sign and the initial
    // decimal digit prefix, and ignores all following bytes. `None` is reserved
    // for inputs whose source signed-int arithmetic overflows and is undefined;
    // it is not a replacement numeric value or a parity claim for those inputs.
    // Complexity review: leading whitespace and the digit prefix are disjoint
    // single forward scans, O(field length) time and O(1) auxiliary storage,
    // matching the source loops without allocation.
    let mut index = 0;
    while index < field.len() && gemmi_is_space(field[index]) {
        index += 1;
    }

    let negative = match field.get(index) {
        Some(b'-') => {
            index += 1;
            true
        }
        Some(b'+') => {
            index += 1;
            false
        }
        _ => false,
    };

    let mut value = 0_i32;
    while let Some(&digit) = field.get(index).filter(|&&byte| gemmi_is_digit(byte)) {
        value = value
            .checked_mul(10)?
            .checked_sub(i32::from(digit - b'0'))?;
        index += 1;
    }

    if negative {
        Some(value)
    } else {
        value.checked_neg()
    }
}

fn read_double(field: &[u8]) -> f64 {
    // Gemmi source: third_party/gemmi/src/pdb.cpp
    // Gemmi❗❗: double read_double(const char* p, int field_length) {
    // Gemmi❗❗:   double d = 0.;
    // Gemmi❗❗:   // we don't check for errors here
    // Gemmi❗❗:   fast_from_chars(p, p + field_length, d);
    // Gemmi❗❗:   return d;
    // Gemmi❗❗: }
    // Gemmi source: third_party/gemmi/include/gemmi/atof.hpp
    // Gemmi❗❗:   while (start < end && is_space(*start))
    // Gemmi❗❗:     ++start;
    // Gemmi❗❗:   if (start < end && *start == '+')
    // Gemmi❗❗:     ++start;
    // Gemmi❗❗:   return fast_float::from_chars(start, end, d);
    // fast_float source: third_party/gemmi/include/gemmi/third_party/fast_float.h
    // fast_float❗❗:   general = fixed | scientific,
    // fast_float❗❗: from_chars(UC const *first, UC const *last, T &value,
    // fast_float❗❗:           chars_format fmt /*= chars_format::general*/) noexcept {
    // fast_float❗❗:   to_float(pns.negative, am, value);
    // fast_float❗❗:   // Test for over/underflow.
    // fast_float❗❗:   if ((pns.mantissa != 0 && am.mantissa == 0 && am.power2 == 0) ||
    // fast_float❗❗:       am.power2 == binary_format<T>::infinite_power()) {
    // fast_float❗❗:     answer.ec = std::errc::result_out_of_range;
    // fast_float❗❗:   }
    // fast_float❗❗:   return answer;
    // Behavior review: preserve Gemmi's C-locale prefix whitespace and its
    // one-leading-plus removal before parsing. An invalid token leaves the
    // initialized +0.0 untouched; a valid numeric prefix is converted even
    // when the source status reports range error, because read_double ignores
    // that status. The prefix scanner also models the source general-format
    // incomplete-exponent rewind. Parity remains unverified pending P06 tests.
    // Complexity review: fixed PDB fields are at most 50 bytes. This path has
    // no heap allocation and O(n) time/O(1) auxiliary state, but Rust first
    // scans for the source prefix and then `FromStr` parses that prefix, so it
    // performs an additional linear pass relative to fast_float. The short,
    // fixed field bounds show no clear material cost difference by inspection;
    // the performance marker remains unresolved after that review.
    // BIO-CID-NUM N02: delegated to the canonical bio-gated numeric owner.
    crate::bio_numeric::gemmi_fast_atof_with_end(field).0
}

fn read_matrix(transform: &mut BioTransform, line: &[u8]) -> i32 {
    // Gemmi source: third_party/gemmi/src/pdb.cpp
    // Gemmi❗✔️: int read_matrix(Transform& t, const char* line, size_t len) {
    // Gemmi❗✔️:   if (len < 46)
    // Gemmi❗✔️:     return 0;
    // Gemmi❗✔️:   char n = line[5] - '0';
    // Gemmi❗✔️:   if (n >= 1 && n <= 3) {
    // Gemmi❗✔️:     t.mat[n-1][0] = read_double(line+10, 10);
    // Gemmi❗✔️:     t.mat[n-1][1] = read_double(line+20, 10);
    // Gemmi❗✔️:     t.mat[n-1][2] = read_double(line+30, 10);
    // Gemmi❗✔️:     t.vec.at(n-1) = read_double(line+45, 10);
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   return n;
    // Gemmi❗✔️: }
    // Behavior review: the short-line guard returns before mutation; only
    // suffix values 1..=3 assign the corresponding matrix row and translation
    // component, in source order. Other suffixes are returned without changing
    // the transform. Each numeric field reuses the existing Gemmi-shaped
    // `read_double`; a short final field ends at the line's NUL boundary.
    // Complexity review: four fixed-width numeric reads and a fixed 3x3 plus
    // three-value copy/reconstruction; no allocation or input-size-dependent
    // scan beyond the four bounded fields.
    if line.len() < 46 {
        return 0;
    }

    // The PDB row suffix is an ASCII byte in the source format. Promote before
    // subtraction so the returned value matches `char n = line[5] - '0'` over
    // the defined ASCII record-label domain without Rust u8 underflow.
    let row_number = i32::from(line[5]) - i32::from(b'0');
    if (1..=3).contains(&row_number) {
        let row = (row_number - 1) as usize;
        let mut matrix = *transform.matrix();
        let mut translation = *transform.translation();
        matrix[row][0] = read_double(&line[10..20]);
        matrix[row][1] = read_double(&line[20..30]);
        matrix[row][2] = read_double(&line[30..40]);
        translation[row] = read_double(&line[45..line.len().min(55)]);
        *transform = BioTransform::new(matrix, translation);
    }
    row_number
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct PdbChargeError {
    digit: u8,
    sign: u8,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum PdbResidueKeyError {
    SequenceIdNotRepresentable,
    ResidueNameTooWide { width: usize },
    ResidueNameNotAscii { byte: u8 },
    SegmentTooWide { width: usize },
    SegmentNotAscii { byte: u8 },
    SegmentNotRepresentable { width: usize },
}

impl std::fmt::Display for PdbResidueKeyError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::SequenceIdNotRepresentable => {
                formatter.write_str("PDB residue sequence identifier is not representable")
            }
            Self::ResidueNameTooWide { width } => write!(
                formatter,
                "residue name width {width} exceeds the BIO representation"
            ),
            Self::ResidueNameNotAscii { byte } => {
                write!(
                    formatter,
                    "residue name contains non-ASCII byte {byte:#04x}"
                )
            }
            Self::SegmentTooWide { width } => {
                write!(
                    formatter,
                    "residue segment width {width} exceeds the BIO representation"
                )
            }
            Self::SegmentNotAscii { byte } => {
                write!(
                    formatter,
                    "residue segment contains non-ASCII byte {byte:#04x}"
                )
            }
            Self::SegmentNotRepresentable { width } => write!(
                formatter,
                "residue segment width {width} is not representable"
            ),
        }
    }
}

impl std::error::Error for PdbResidueKeyError {}

pub(super) fn mapped_ccd_residue_address(
    address: ResidueAddress,
    aliases: &[PdbCcdAlias],
) -> Result<ResidueAddress, PdbResidueKeyError> {
    // Gemmi source helper: third_party/gemmi/include/gemmi/modify.hpp.
    // Gemmi❗✔️:   auto update = [&](ResidueId& rid) {
    // Gemmi❗✔️:     if (rid.name == old_name)
    // Gemmi❗✔️:       rid.name = new_name;
    // Gemmi❗✔️:   };
    // Behavior review: apply the retained aliases in encounter order, matching
    // each source residue name exactly. Prevalidation and commit call this same
    // deterministic mapper; a full code outside the approved four-byte ASCII
    // name representation returns its existing typed boundary error.
    // Complexity review: at most A alias comparisons over a fixed-size name,
    // O(A) time and O(1) auxiliary space.
    let original = address.name();
    let mut resolved = original.as_bytes();
    for alias in aliases {
        if resolved == alias.short_code.as_slice() {
            resolved = &alias.full_code;
        }
    }
    if resolved == original.as_bytes() {
        return Ok(address);
    }

    let name = match ResidueName::from_ascii(resolved) {
        Some(name) => name,
        None => {
            if let Some(byte) = resolved.iter().copied().find(|byte| !byte.is_ascii()) {
                return Err(PdbResidueKeyError::ResidueNameNotAscii { byte });
            }
            return Err(PdbResidueKeyError::ResidueNameTooWide {
                width: resolved.len(),
            });
        }
    };
    let segment = address.segment();
    ResidueAddress::new(
        address.sequence_number(),
        address.insertion_code(),
        segment.as_bytes(),
        name,
    )
    .ok_or(PdbResidueKeyError::SegmentNotRepresentable {
        width: segment.len(),
    })
}

pub(super) fn mapped_ccd_atom_address(
    address: &AtomAddress,
    aliases: &[PdbCcdAlias],
) -> Result<AtomAddress, PdbResidueKeyError> {
    // Gemmi source helper: third_party/gemmi/include/gemmi/modify.hpp.
    // Gemmi❗✔️:   process_addresses(st, [&](AtomAddress& aa) { update(aa.res_id); });
    // Behavior review: rename only the residue-address member; chain, logical
    // atom name and altloc are retained without reinterpretation. The same
    // alias sequence and existing fixed-width boundary apply to every source
    // address family.
    // Complexity review: residue mapping scans A aliases, then clones only the
    // logical atom name, O(A + atom-name bytes) time and O(atom-name bytes)
    // owned storage as the source address replacement requires.
    let residue = mapped_ccd_residue_address(address.residue(), aliases)?;
    Ok(AtomAddress::new(
        address.chain_name(),
        residue,
        address.logical_atom_name().to_owned(),
        Some(address.altloc()),
    ))
}

pub(super) fn rename_ccd_sequence_tokens(sequence: &mut Vec<u8>, old_name: &[u8], new_name: &[u8]) {
    // Gemmi source helper: third_party/gemmi/include/gemmi/modify.hpp.
    // Gemmi✔️✔️:       for (size_t start = 0;;) {
    // Gemmi✔️✔️:         size_t end = mon_ids.find(',', start);
    // Gemmi✔️✔️:         if (mon_ids.compare(start, end-start, old_name) == 0) {
    // Gemmi✔️✔️:           mon_ids.replace(start, end-start, new_name);
    // Gemmi✔️✔️:           if (end != std::string::npos)
    // Gemmi✔️✔️:             end = start + new_name.size();
    // Gemmi✔️✔️:         }
    // Gemmi✔️✔️:         if (end == std::string::npos)
    // Gemmi✔️✔️:           break;
    // Gemmi✔️✔️:         start = end + 1;
    // Gemmi✔️✔️:       }
    // Behavior review: compare and replace exact comma-delimited byte tokens;
    // after replacement, resume after the replaced token and original delimiter
    // exactly as the source does, including empty and repeated tokens.
    // Complexity review: one forward scan per alias; byte movement from Vec
    // splice is linear in the affected suffix, like std::string::replace.
    let mut start = 0;
    loop {
        let delimiter = sequence[start..]
            .iter()
            .position(|byte| *byte == b',')
            .map(|offset| start + offset);
        let end = delimiter.unwrap_or(sequence.len());
        if &sequence[start..end] == old_name {
            sequence
                .splice(start..end, new_name.iter().copied())
                .for_each(drop);
            if delimiter.is_some() {
                let end = start + new_name.len();
                start = end + 1;
                continue;
            }
        }
        let Some(end) = delimiter else {
            break;
        };
        start = end + 1;
    }
}

fn read_altloc(field: u8) -> Option<AltLocLabel> {
    // Gemmi source: third_party/gemmi/src/pdb.cpp
    // Gemmi❗✔️: char read_altloc(char c) { return c == ' ' ? '\0' : c; }
    // Behavior review: translate the source's returned NUL into the canonical
    // optional BIO label; all other field bytes remain unchanged.
    // Complexity review: constant-time byte comparison/construction, no allocation.
    let value = if field == b' ' { 0 } else { field };
    (value != 0).then(|| AltLocLabel::new(value))
}

fn read_res_id(
    sequence_field: &[u8; 5],
    residue_name_field: &[u8; 3],
) -> Result<ResidueAddress, PdbResidueKeyError> {
    // Gemmi source: third_party/gemmi/src/pdb.cpp
    // Gemmi❗✔️: ResidueId read_res_id(const char* seq_id, const char* name) {
    // Gemmi❗✔️:   return {read_seq_id(seq_id), {}, read_string(name, 3)};
    // Gemmi❗✔️: }
    // Behavior review: project the fixed five-byte SeqId and trimmed three-byte
    // residue name into the canonical address. The source's empty segment is
    // retained here; its ATOM/HETATM caller adds a segment only in a separate
    // length-guarded action. INT_MIN and space sentinels are canonicalized by
    // the existing ResidueAddress constructor.
    // Complexity review: all reads and conversion inputs are fixed-width;
    // constant auxiliary space and no allocation.
    let sequence =
        read_seq_id(sequence_field).ok_or(PdbResidueKeyError::SequenceIdNotRepresentable)?;
    let raw_name = read_string(residue_name_field);
    let name = match ResidueName::from_ascii(raw_name) {
        Some(name) => name,
        None => {
            if let Some(byte) = raw_name.iter().copied().find(|byte| !byte.is_ascii()) {
                return Err(PdbResidueKeyError::ResidueNameNotAscii { byte });
            }
            return Err(PdbResidueKeyError::ResidueNameTooWide {
                width: raw_name.len(),
            });
        }
    };

    ResidueAddress::new(Some(sequence.seq_num()), sequence.ins_code(), b"", name)
        .ok_or(PdbResidueKeyError::SegmentNotRepresentable { width: 0 })
}

fn read_pdb_residue_segment(
    residue: ResidueAddress,
    line: &[u8],
) -> Result<ResidueAddress, PdbResidueKeyError> {
    // Gemmi source: third_party/gemmi/src/pdb.cpp, populate_structure_from_pdb_stream
    // Gemmi❗✔️:       if (len > 72)
    // Gemmi❗✔️:         rid.segment = read_string(line+72, 4);
    // Behavior review: keep the caller's strict `len > 72` condition and its
    // four-byte source field. A short line retains the empty segment created
    // by read_res_id; selected bytes pass through the existing source-shaped
    // read_string trimming before the canonical ASCII/width boundary check.
    // Complexity review: at most four field bytes are inspected, so O(1) time
    // and auxiliary space; no heap allocation.
    if line.len() <= 72 {
        return Ok(residue);
    }

    let segment_end = line.len().min(76);
    let segment = read_string(&line[72..segment_end]);
    if segment.len() > 4 {
        return Err(PdbResidueKeyError::SegmentTooWide {
            width: segment.len(),
        });
    }
    if let Some(byte) = segment.iter().copied().find(|byte| !byte.is_ascii()) {
        return Err(PdbResidueKeyError::SegmentNotAscii { byte });
    }

    ResidueAddress::new(
        residue.sequence_number(),
        residue.insertion_code(),
        segment,
        residue.name(),
    )
    .ok_or(PdbResidueKeyError::SegmentNotRepresentable {
        width: segment.len(),
    })
}

fn read_charge(digit: u8, sign: u8) -> Result<i8, PdbChargeError> {
    // Gemmi source: third_party/gemmi/src/pdb.cpp
    // Gemmi❗✔️: // The standard charge format is 2+, but some files have +2.
    // Gemmi❗✔️: signed char read_charge(char digit, char sign) {
    // Gemmi❗✔️:   if (sign == ' ' && digit == ' ')  // by far the most common case
    // Gemmi❗✔️:     return 0;
    // Gemmi❗✔️:   if (sign >= '0' && sign <= '9')
    // Gemmi❗✔️:     std::swap(digit, sign);
    // Gemmi❗✔️:   if (digit >= '0' && digit <= '9') {
    // Gemmi❗✔️:     if (sign != '+' && sign != '-' && sign != '\0' && !is_space(sign))
    // Gemmi❗✔️:       fail("Wrong format for charge: " +
    // Gemmi❗✔️:            std::string(1, digit) + std::string(1, sign));
    // Gemmi❗✔️:     return (digit - '0') * (sign == '-' ? -1 : 1);
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   // if we are here the field should be blank, but maybe better not to check
    // Gemmi❗✔️:   return 0;
    // Gemmi❗✔️: }
    // Behavior review: preserve the source's order-sensitive second-byte digit
    // swap, digit gate, and sign validation only within that gate. A non-digit
    // first byte returns zero even if the second byte is otherwise malformed;
    // an invalid sign after a digit returns both original bytes in a typed
    // error. The result is explicitly i8, matching Gemmi's signed-char
    // destination without constraining the source-accepted one-digit range.
    // Complexity review: one constant-time optional swap and bounded byte
    // predicates, O(1) time and space; no success-path allocation. The typed
    // error retains both bytes without the source's temporary string creation.
    if sign == b' ' && digit == b' ' {
        return Ok(0);
    }

    let (mut digit, mut sign) = (digit, sign);
    if gemmi_is_digit(sign) {
        std::mem::swap(&mut digit, &mut sign);
    }
    if gemmi_is_digit(digit) {
        if sign != b'+' && sign != b'-' && sign != b'\0' && !gemmi_is_space(sign) {
            return Err(PdbChargeError { digit, sign });
        }
        let magnitude = (digit - b'0') as i8;
        return Ok(if sign == b'-' { -magnitude } else { magnitude });
    }
    Ok(0)
}

fn read_base36<const N: usize>(field: &[u8; N]) -> Option<i32> {
    // Gemmi source: third_party/gemmi/src/pdb.cpp
    // Gemmi❗✔️: template<int N> int read_base36(const char* p) {
    // Gemmi❗✔️:   char zstr[N+1] = {0};
    // Gemmi❗✔️:   std::memcpy(zstr, p, N);
    // Gemmi❗✔️:   return std::strtol(zstr, nullptr, 36);
    // Gemmi❗✔️: }
    // Behavior review: for the only PDB field widths used here (four or five
    // bytes), `strtol`'s base-36 digit prefix fits both its `long` and Gemmi's
    // returned `int`. It accepts 0-9/A-Z/a-z case-insensitively, stops at the
    // first non-digit or NUL, and yields zero when no digit converts. `None`
    // only marks a wider private use that cannot be represented here; no
    // such width is used by the PDB source paths.
    // Complexity review: one bounded forward scan and constant auxiliary
    // storage, O(N) time/O(1) space like `strtol`; unlike the C++ template it
    // does not copy the fixed field into a temporary NUL-terminated array.
    let mut index = 0;
    while index < N && gemmi_is_space(field[index]) {
        index += 1;
    }

    let negative = match field.get(index) {
        Some(b'-') => {
            index += 1;
            true
        }
        Some(b'+') => {
            index += 1;
            false
        }
        _ => false,
    };

    let mut value = 0_i64;
    let mut has_digits = false;
    while let Some(&byte) = field.get(index).filter(|&&byte| byte != 0) {
        let digit = match byte {
            b'0'..=b'9' => i64::from(byte - b'0'),
            b'A'..=b'Z' => i64::from(byte - b'A' + 10),
            b'a'..=b'z' => i64::from(byte - b'a' + 10),
            _ => break,
        };
        value = value.checked_mul(36)?.checked_add(digit)?;
        has_digits = true;
        index += 1;
    }

    if !has_digits {
        return Some(0);
    }

    let signed = if negative { -value } else { value };
    i32::try_from(signed).ok()
}

fn read_serial(field: &[u8; 5]) -> Option<i32> {
    // Gemmi source: third_party/gemmi/src/pdb.cpp
    // Gemmi❗✔️: int read_serial(const char* ptr) {
    // Gemmi❗✔️:   return ptr[0] < 'A' ? read_int(ptr, 5)
    // Gemmi❗✔️:                       : read_base36<5>(ptr) - 16796160 + 100000;
    // Gemmi❗✔️: }
    // The source compiler uses signed plain `char`; reproduce its comparison
    // for arbitrary bytes while keeping ASCII serial behavior unchanged.
    let first = i8::from_ne_bytes([field[0]]);
    if first < b'A' as i8 {
        read_int(field.as_slice())
    } else {
        read_base36(field)?
            .checked_sub(16_796_160)?
            .checked_add(100_000)
    }
}

fn read_seq_id(field: &[u8; 5]) -> Option<PdbSeqId> {
    // Gemmi source: third_party/gemmi/src/pdb.cpp
    // Gemmi❗✔️: SeqId read_seq_id(const char* str) {
    // Gemmi❗✔️:   SeqId seqid;
    // Gemmi❗✔️:   if (str[4] != '\r' && str[4] != '\n')
    // Gemmi❗✔️:     seqid.icode = str[4];
    // Gemmi❗✔️:   // We support hybrid-36 extension, although it is never used in practice
    // Gemmi❗✔️:   // as 9999 residues per chain are enough.
    // Gemmi❗✔️:   if (str[0] < 'A') {
    // Gemmi❗✔️:     for (int i = 4; i != 0; --i, ++str)
    // Gemmi❗✔️:       if (!is_space(*str)) {
    // Gemmi❗✔️:         seqid.num = read_int(str, i);
    // Gemmi❗✔️:         break;
    // Gemmi❗✔️:       }
    // Gemmi❗✔️:   } else {
    // Gemmi❗✔️:     seqid.num = read_base36<4>(str) - 466560 + 10000;
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   return seqid;
    // Gemmi❗✔️: }
    // Behavior review: the fixed five-byte field stores four sequence bytes
    // and one insertion byte. An all-whitespace decimal prefix preserves
    // Gemmi's INT_MIN absent-number sentinel; a non-space invalid prefix is
    // passed to read_int and therefore yields zero. The fifth byte is raw
    // except CR/LF, while blank space maps to the canonical absent insertion.
    // Complexity review: both source loops are bounded by the four-byte
    // sequence field; this implementation performs no allocation and O(1)
    // auxiliary work.
    let insertion = match field[4] {
        b'\r' | b'\n' => b' ',
        byte => byte,
    };
    let insertion_code = (insertion != b' ').then_some(insertion);

    // Gemmi compares plain `char`; the pinned Linux source build uses signed
    // char, so preserve its first-byte dispatch for arbitrary PDB bytes.
    let sequence_number = if i8::from_ne_bytes([field[0]]) < b'A' as i8 {
        let mut sequence_number = i32::MIN;
        for width in (1..=4).rev() {
            let offset = 4 - width;
            if !gemmi_is_space(field[offset]) {
                sequence_number = read_int(&field[offset..offset + width])?;
                break;
            }
        }
        sequence_number
    } else {
        let sequence_field: [u8; 4] = field[..4].try_into().ok()?;
        read_base36(&sequence_field)?
            .checked_sub(466_560)?
            .checked_add(10_000)?
    };

    Some(PdbSeqId::new(sequence_number, insertion_code))
}

fn read_string(field: &[u8]) -> &[u8] {
    // Gemmi source: third_party/gemmi/src/pdb.cpp
    // Gemmi❗✔️: std::string read_string(const char* p, int field_length) {
    // Gemmi❗✔️:   // left trim
    // Gemmi❗✔️:   while (field_length != 0 && is_space(*p)) {
    // Gemmi❗✔️:     ++p;
    // Gemmi❗✔️:     --field_length;
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   // EOL/EOF ends the string
    // Gemmi❗✔️:   for (int i = 0; i < field_length; ++i)
    // Gemmi❗✔️:     if (p[i] == '\n' || p[i] == '\r' || p[i] == '\0') {
    // Gemmi❗✔️:       field_length = i;
    // Gemmi❗✔️:       break;
    // Gemmi❗✔️:     }
    // Gemmi❗✔️:   // right trim
    // Gemmi❗✔️:   while (field_length != 0 && is_space(p[field_length-1]))
    // Gemmi❗✔️:     --field_length;
    // Gemmi❗✔️:   return std::string(p, field_length);
    // Gemmi❗✔️: }
    let mut start = 0;
    while start < field.len() && gemmi_is_space(field[start]) {
        start += 1;
    }

    let remaining = &field[start..];
    let end = remaining
        .iter()
        .position(|byte| matches!(*byte, b'\n' | b'\r' | 0))
        .unwrap_or(remaining.len());
    let terminated = &remaining[..end];
    let trimmed_end = terminated
        .iter()
        .rposition(|byte| !gemmi_is_space(*byte))
        .map_or(0, |index| index + 1);
    &terminated[..trimmed_end]
}

fn raw_atom_name(field: &[u8]) -> Option<AtomName> {
    let bytes: [u8; 4] = field.try_into().ok()?;
    AtomName::from_ascii(&bytes)
}

struct PdbLineCursor<'a> {
    input: &'a [u8],
    position: usize,
    max_line_length: usize,
    line_buffer: [u8; 122],
}

impl<'a> PdbLineCursor<'a> {
    fn new(input: &'a str, options: PdbInputOptions) -> Self {
        Self {
            input: input.as_bytes(),
            position: 0,
            max_line_length: options.effective_max_line_length(),
            // Gemmi source: third_party/gemmi/src/pdb.cpp,
            // populate_structure_from_pdb_stream.
            // Gemmi❗✔️:   char line[122] = {0};
            line_buffer: [0; 122],
        }
    }

    fn copy_line(&mut self) -> Option<&[u8]> {
        // Gemmi source: third_party/gemmi/include/gemmi/input.hpp
        // Gemmi❗✔️:   size_t copy_line(char* line, int size) {
        // Gemmi❗✔️:     if (!gets(line, size))
        // Gemmi❗✔️:       return 0;
        // Gemmi❗✔️:     size_t len = std::strlen(line);
        // Gemmi❗✔️:     // If a line is longer than size we discard the rest of it.
        // Gemmi❗✔️:     if (len > 0 && line[len-1] != '\n')
        // Gemmi❗✔️:       for (int c = getc(); c > 0 /* not 0 nor EOF */ && c != '\n'; c = getc())
        // Gemmi❗✔️:         continue;
        // Gemmi❗✔️:     return len;
        // Gemmi❗✔️:   };
        // Gemmi❗✔️:   char* gets(char* line, int size) override {
        // Gemmi❗✔️:     --size; // fgets reads in at most one less than size characters
        // Gemmi❗✔️:     if (cur >= end)
        // Gemmi❗✔️:       return nullptr;
        // Gemmi❗✔️:     if (size > end - cur)
        // Gemmi❗✔️:       size = int(end - cur);
        // Gemmi❗✔️:     const char* nl = (const char*) std::memchr(cur, '\n', size);
        // Gemmi❗✔️:     size_t len = nl ? nl - cur + 1 : size;
        // Gemmi❗✔️:     std::memcpy(line, cur, len);
        // Gemmi❗✔️:     line[len] = '\0';
        // Gemmi❗✔️:     cur += len;
        // Gemmi❗✔️:     return line;
        // Gemmi❗✔️:   }
        // Gemmi❗✔️:   int getc() override { return cur < end ? *cur++ : EOF; }
        //
        // Rust mirrors the reusable fixed buffer: a current line overwrites
        // only its copied prefix and its terminating NUL, leaving any suffix
        // from the previous line untouched just as `fgets` does. The cursor
        // performs the same bounded prefix copy and overlong-line drain with
        // O(1) storage and no per-record allocation.

        if self.position == self.input.len() {
            return None;
        }

        let start = self.position;
        let buffer_end = start + self.max_line_length.min(self.input.len() - start);
        let fetched = &self.input[start..buffer_end];
        let copied_end = fetched
            .iter()
            .position(|byte| *byte == b'\n')
            .map_or(buffer_end, |offset| start + offset + 1);
        self.position = copied_end;

        let copied_len = copied_end - start;
        self.line_buffer[..copied_len].copy_from_slice(&self.input[start..copied_end]);
        self.line_buffer[copied_len] = 0;
        let visible_len = self.line_buffer[..=copied_len]
            .iter()
            .position(|byte| *byte == 0)
            .unwrap_or(copied_len);
        if visible_len == 0 {
            return None;
        }

        if self.line_buffer[visible_len - 1] != b'\n' {
            while self.position < self.input.len() {
                // MemoryStream::getc returns `char` promoted to int. The pinned
                // Linux source environment uses signed char, so high-bit bytes
                // terminate this `c > 0` drain condition just as in Gemmi.
                let c = self.input[self.position] as i8;
                self.position += 1;
                if c <= 0 || c == b'\n' as i8 {
                    break;
                }
            }
        }

        Some(&self.line_buffer[..visible_len])
    }

    fn source_line_buffer(&self) -> &[u8; 122] {
        &self.line_buffer
    }
}

#[cfg(test)]
mod bio_cid_num_n01_tests {
    use crate::bio_numeric::parse_fast_float_prefix;

    // Expectations derived from pinned fast_float.h:4469-4522
    // (parse_infnan) and :4760-4774 (from_chars_advanced fallback), not
    // from the implementation under test.
    #[test]
    fn bio_cid_num_n01_inf_spellings_and_offsets() {
        // +inf / -inf: sign consumed inside parse_infnan, value signed.
        assert_eq!(parse_fast_float_prefix(b"inf"), Some((f64::INFINITY, 3)));
        assert_eq!(
            parse_fast_float_prefix(b"-inf"),
            Some((f64::NEG_INFINITY, 4))
        );
        // Case-insensitive 3-letter and full 8-letter spellings.
        assert_eq!(parse_fast_float_prefix(b"INF"), Some((f64::INFINITY, 3)));
        assert_eq!(parse_fast_float_prefix(b"iNf"), Some((f64::INFINITY, 3)));
        assert_eq!(
            parse_fast_float_prefix(b"Infinity"),
            Some((f64::INFINITY, 8))
        );
        assert_eq!(
            parse_fast_float_prefix(b"-INFINITY"),
            Some((f64::NEG_INFINITY, 9))
        );
        // Partial "infinit" (7 bytes) takes the short form: ptr = +3.
        assert_eq!(
            parse_fast_float_prefix(b"infinit"),
            Some((f64::INFINITY, 3))
        );
        // "infin" is 5 bytes: "ini" does not match -> the inf branch
        // fails (not strncasecmp "inf"? it does match inf...) -> 3 bytes
        // then trailing 'i','n' left unconsumed.
        assert_eq!(parse_fast_float_prefix(b"infin"), Some((f64::INFINITY, 3)));
    }

    #[test]
    fn bio_cid_num_n01_nan_spellings_and_payloads() {
        // Plain nan, all case classes, with sign.
        assert!(matches!(parse_fast_float_prefix(b"nan"), Some((v, 3)) if v.is_nan()));
        assert!(matches!(parse_fast_float_prefix(b"NAN"), Some((v, 3)) if v.is_nan()));
        assert!(matches!(parse_fast_float_prefix(b"-NaN"), Some((v, 4)) if v.is_nan()));
        // Closed alnum/underscore payload consumes through ')'.
        assert!(matches!(parse_fast_float_prefix(b"nan(ind)"), Some((v, 8)) if v.is_nan()));
        assert!(matches!(parse_fast_float_prefix(b"nan(x_9)"), Some((v, 8)) if v.is_nan()));
        // Unclosed payload: end stays at +3 (payload not consumed), the
        // NaN is still assigned with success status.
        assert!(matches!(parse_fast_float_prefix(b"nan(x"), Some((v, 3)) if v.is_nan()));
        // Forbidden payload character stops the scan at the same +3.
        assert!(matches!(parse_fast_float_prefix(b"nan(x!)"), Some((v, 3)) if v.is_nan()));
        // Empty payload "nan()" closes immediately at 5.
        assert!(matches!(parse_fast_float_prefix(b"nan()"), Some((v, 5)) if v.is_nan()));
        // "na" alone (fewer than 3 bytes) is not special: decimal scan.
        assert_eq!(parse_fast_float_prefix(b"na"), None);
    }

    #[test]
    fn bio_cid_num_n01_decimal_prefix_still_default() {
        // Default flags keep the ordinary decimal grammar ahead of the
        // inf/nan fallback: a valid decimal prefix consumes exactly the
        // digits; the fallback is reached only when it fails.
        assert_eq!(parse_fast_float_prefix(b"1.5e2x"), Some((150.0, 5)));
        assert_eq!(parse_fast_float_prefix(b"-.25"), Some((-0.25, 4)));
        // Incomplete exponent rolls back to the 'e'.
        assert_eq!(parse_fast_float_prefix(b"1e+"), Some((1.0, 1)));
        // Nothing decimal and nothing special -> None (source reports
        // invalid_argument through the wrapper).
        assert_eq!(parse_fast_float_prefix(b"x"), None);
    }
}

#[cfg(test)]
mod bio_cid_num_n09_tests {
    use super::read_double;

    // Fixed-field read_double regressions through the typed numeric
    // owner: the pinned source ignores the error status and returns d
    // (pdb.cpp:26-31 "we don't check for errors here"), so a range error
    // still yields the ASSIGNED extreme and a blank/invalid field the
    // initial +0.0.
    #[test]
    fn bio_cid_num_n09_fixed_field_statuses() {
        // Success.
        assert_eq!(read_double(b"  1.5    "), 1.5);
        assert_eq!(read_double(b"-2.25e2  "), -225.0);
        // Blank and invalid fields keep the initialized +0.0.
        assert_eq!(read_double(b"         "), 0.0);
        assert_eq!(read_double(b"  junk   "), 0.0);
        // Prefix conversion inside a fixed field stops at the field.
        assert_eq!(read_double(b"  1.5x   "), 1.5);
        // Range error with assigned infinity (status deliberately
        // ignored by read_double, exactly as the source).
        assert_eq!(read_double(b" 1.8e308 "), f64::INFINITY);
        assert_eq!(read_double(b"-1.8e308 "), f64::NEG_INFINITY);
        // Nonzero underflow: assigned signed zero, status ignored.
        assert_eq!(read_double(b" 1e-400  "), 0.0);
        // Literal specials inside a fixed field.
        assert_eq!(read_double(b"  inf    "), f64::INFINITY);
        assert!(read_double(b"  nan    ").is_nan());
    }
}

#[cfg(test)]
mod tests {
    use super::{
        PdbAtomFieldError, PdbBioReader, PdbCcdAlias, PdbChargeError, PdbCisPepRecordError,
        PdbCrystalError, PdbDbRefAction, PdbEntity, PdbEntityState, PdbGroupLocation,
        PdbHeaderError, PdbHelixError, PdbHierarchyGrouping, PdbInputFormatError, PdbInputOptions,
        PdbLineCursor, PdbLinkNameCounters, PdbModResError, PdbModelGroup, PdbModelTransition,
        PdbModelTransitionError, PdbMtrixError, PdbReaderCause, PdbReaderError, PdbReaderStage,
        PdbReaderState, PdbRemark3Continuation, PdbRemark3Error, PdbRemark350Error, PdbRemarkError,
        PdbRemarkMetadataError, PdbResidueKeyError, PdbSeqId, PdbSsbondResolutionError,
        assign_pdb_subchain_names, assign_pdb_subchains, change_author_name_format_to_mmcif,
        complete_ssbond, decode_pdb_atom_fields, decode_pdb_cispep_record, decode_pdb_link_address,
        decode_pdb_link_record, infer_element_from_padded_name, pdb_date_format_to_iso,
        raw_atom_name, read_altloc, read_charge, read_double, read_int, read_matrix,
        read_pdb_residue_segment, read_res_id, read_seq_id, read_serial, read_string,
        resolve_pdb_atom_element,
    };
    use cosmolkit_bio::{
        AtomAddress, AtomName, BioAssemblySpecialKind, BioAsu, BioCisPep, BioConnection,
        BioConnectionKind, BioCoordinateFormat, BioCrystalCell, BioCrystalInfo, BioHelix,
        BioHelixClass, BioModRes, BioSoftwareClassification, BioStructureError,
        BioStructureSourceState, BioTlsGroup, BioTlsSelection, BioTransform, EntityKind,
        PdbChainId, ResidueAddress, ResidueName,
    };
    use cosmolkit_types::Element;

    fn feed_remark3(state: &mut PdbReaderState, indentation: usize, payload: &str) {
        let line = format!("REMARK   3 {}{payload}", " ".repeat(indentation));
        state
            .remark3_record(&line)
            .expect("source-backed test remark is representable");
    }

    fn remark350_blank_line() -> Vec<u8> {
        let mut line = vec![b' '; 80];
        line[..11].copy_from_slice(b"REMARK 350 ");
        line
    }

    fn remark350_biomolecule(name: &str) -> String {
        assert!(name.len() <= 20);
        let mut line = remark350_blank_line();
        line[11..23].copy_from_slice(b"BIOMOLECULE:");
        line[23..23 + name.len()].copy_from_slice(name.as_bytes());
        String::from_utf8(line).expect("test REMARK text is ASCII")
    }

    fn remark350_apply(chains: &str, first: bool) -> String {
        assert!(chains.len() <= 39);
        let mut line = remark350_blank_line();
        let prefix = if first {
            b"APPLY THE FOLLOWING TO CHAINS".as_slice()
        } else {
            b"                   AND CHAINS".as_slice()
        };
        assert_eq!(prefix.len(), 29);
        line[11..40].copy_from_slice(prefix);
        line[40] = b':';
        line[41..41 + chains.len()].copy_from_slice(chains.as_bytes());
        String::from_utf8(line).expect("test REMARK text is ASCII")
    }

    fn remark350_metadata_row(
        key: &str,
        colon_offset: usize,
        value_offset: usize,
        value: &str,
    ) -> String {
        assert!(11 + key.len() <= colon_offset);
        assert!(value_offset + value.len() <= 80);
        let mut line = remark350_blank_line();
        line[11..11 + key.len()].copy_from_slice(key.as_bytes());
        line[colon_offset] = b':';
        line[value_offset..value_offset + value.len()].copy_from_slice(value.as_bytes());
        String::from_utf8(line).expect("test REMARK text is ASCII")
    }

    fn remark350_biomt(row: u8, operator_id: &str, values: [f64; 4]) -> String {
        assert!((b'1'..=b'3').contains(&row));
        assert_eq!(operator_id.len(), 3);
        let mut line = remark350_blank_line();
        line[11..18].copy_from_slice(b"  BIOMT");
        line[18] = row;
        line[20..23].copy_from_slice(operator_id.as_bytes());
        for (offset, value) in [23, 33, 43, 58].into_iter().zip(values) {
            let field = format!("{value:>10.6}");
            assert_eq!(field.len(), 10);
            line[offset..offset + 10].copy_from_slice(field.as_bytes());
        }
        String::from_utf8(line).expect("test REMARK text is ASCII")
    }

    fn feed_remark200(state: &mut PdbReaderState, number: u16, payload: &str) {
        let line = format!("REMARK {number:3} {payload}");
        state
            .remark_200_230_240_record(&line)
            .expect("source-backed experimental remark is representable");
    }

    fn dbref_source_line(variant: u8, chain: [u8; 2]) -> [u8; 122] {
        let mut line = [0; 122];
        line[..5].copy_from_slice(b"DBREF");
        line[5] = variant;
        line[11..13].copy_from_slice(&chain);
        line
    }

    fn modres_source_line(
        record: [u8; 6],
        chain: [u8; 2],
        sequence: [u8; 5],
        residue: [u8; 3],
        parent: [u8; 3],
        details: &[u8],
        mod_id: &[u8],
        line_len: usize,
        extension_guard: [u8; 2],
    ) -> [u8; 122] {
        assert!(line_len <= 80);
        assert!(details.len() <= 41);
        assert!(mod_id.len() <= 8);
        let mut line = [0; 122];
        line[..80].fill(b' ');
        line[..6].copy_from_slice(&record);
        line[7..11].copy_from_slice(b"1ABC");
        line[12..15].copy_from_slice(&residue);
        line[15..17].copy_from_slice(&chain);
        line[18..23].copy_from_slice(&sequence);
        line[24..27].copy_from_slice(&parent);
        line[29..29 + details.len()].copy_from_slice(details);
        line[70..72].copy_from_slice(&extension_guard);
        line[72..72 + mod_id.len()].copy_from_slice(mod_id);
        line[line_len..].fill(0);
        line
    }

    fn helix_source_line(
        record: [u8; 6],
        start_chain: [u8; 2],
        start_sequence: [u8; 5],
        start_residue: [u8; 3],
        end_chain: [u8; 2],
        end_sequence: [u8; 5],
        end_residue: [u8; 3],
        class: [u8; 2],
        length: [u8; 5],
        line_len: usize,
    ) -> [u8; 122] {
        assert!(line_len <= 80);
        let mut line = [0; 122];
        line[..80].fill(b' ');
        line[..6].copy_from_slice(&record);
        line[15..18].copy_from_slice(&start_residue);
        line[18..20].copy_from_slice(&start_chain);
        line[21..26].copy_from_slice(&start_sequence);
        line[27..30].copy_from_slice(&end_residue);
        line[30..32].copy_from_slice(&end_chain);
        line[33..38].copy_from_slice(&end_sequence);
        line[38..40].copy_from_slice(&class);
        line[72..77].copy_from_slice(&length);
        line[line_len..].fill(0);
        line
    }

    fn sheet_source_line(record: [u8; 6], line_len: usize, fields: &[(usize, &[u8])]) -> [u8; 122] {
        assert!(line_len <= 80);
        let mut line = [0; 122];
        line[..80].fill(b' ');
        line[..6].copy_from_slice(&record);
        for (offset, value) in fields {
            line[*offset..*offset + value.len()].copy_from_slice(value);
        }
        line[line_len..].fill(0);
        line
    }

    fn expected_helix_address(
        chain: &[u8],
        sequence: i32,
        insertion: u8,
        residue_name: &[u8],
    ) -> AtomAddress {
        let chain = PdbChainId::from_ascii(chain).unwrap();
        let name = ResidueName::from_ascii(residue_name).unwrap();
        let residue =
            cosmolkit_bio::ResidueAddress::new(Some(sequence), Some(insertion), b"", name).unwrap();
        AtomAddress::new(chain, residue, "", None)
    }

    fn expected_sheet_address(
        chain: &[u8],
        sequence: i32,
        insertion: u8,
        residue_name: &[u8],
        atom_name: &str,
    ) -> AtomAddress {
        let chain = PdbChainId::from_ascii(chain).unwrap();
        let name = ResidueName::from_ascii(residue_name).unwrap();
        let residue =
            cosmolkit_bio::ResidueAddress::new(Some(sequence), Some(insertion), b"", name).unwrap();
        AtomAddress::new(chain, residue, atom_name, None)
    }

    fn hetnam_source_line(
        record: [u8; 6],
        short_code: [u8; 3],
        full_code: &[u8],
        line_len: usize,
        column_70: u8,
    ) -> [u8; 122] {
        assert!(line_len < 122);
        assert!(full_code.len() <= 8);
        let mut line = [0; 122];
        line[..80].fill(b' ');
        line[..6].copy_from_slice(&record);
        line[11..14].copy_from_slice(&short_code);
        line[70] = column_70;
        let available = line_len.saturating_sub(71).min(8).min(full_code.len());
        line[71..71 + available].copy_from_slice(&full_code[..available]);
        line[line_len] = 0;
        line
    }

    fn set_dbref_source_field(line: &mut [u8; 122], start: usize, width: usize, value: &[u8]) {
        assert!(value.len() <= width);
        line[start..start + width].fill(b' ');
        line[start..start + value.len()].copy_from_slice(value);
    }

    fn conect_source_line(source_serial: i32, targets: &[i32], line_len: usize) -> [u8; 122] {
        assert!(line_len <= 122);
        let mut line = [0; 122];
        line[..6].copy_from_slice(b"CONECT");
        let source_field = format!("{source_serial:>5}");
        assert_eq!(source_field.len(), 5);
        line[6..11].copy_from_slice(source_field.as_bytes());
        for (index, target) in targets.iter().enumerate() {
            let offset = 11 + index * 5;
            let target_field = format!("{target:>5}");
            assert_eq!(target_field.len(), 5);
            line[offset..offset + 5].copy_from_slice(target_field.as_bytes());
        }
        line[line_len..].fill(0);
        line
    }

    #[test]
    fn gemmi_pdb_read_altloc_maps_only_space_and_nul_to_absence() {
        // Gemmi returns NUL only for a blank field; the canonical optional
        // label represents that source absence without normalizing other bytes.
        assert_eq!(read_altloc(b' '), None);
        assert_eq!(read_altloc(0), None);
        assert_eq!(read_altloc(b'A').map(|label| label.value()), Some(b'A'));
        assert_eq!(read_altloc(0x80).map(|label| label.value()), Some(0x80));
    }

    #[test]
    fn gemmi_pdb_read_res_id_preserves_seqid_insertion_and_trimmed_name() {
        let residue = read_res_id(b"  12A", b" GL").unwrap();
        assert_eq!(residue.sequence_number(), Some(12));
        assert_eq!(residue.insertion_code(), Some(b'A'));
        assert_eq!(residue.segment(), "");
        assert_eq!(residue.name().as_str(), "GL");

        // Gemmi's SeqId uses INT_MIN and space as absence sentinels, while
        // CR/LF suppress insertion but an explicit NUL remains a byte value.
        let blank = read_res_id(b"     ", b"ALA").unwrap();
        assert_eq!(blank.sequence_number(), None);
        assert_eq!(blank.insertion_code(), None);
        assert_eq!(
            read_res_id(b"  12\n", b"ALA").unwrap().insertion_code(),
            None
        );
        assert_eq!(
            read_res_id(&[b' ', b' ', b'1', b'2', 0], b"ALA")
                .unwrap()
                .insertion_code(),
            Some(0)
        );
        assert_eq!(
            read_res_id(b"  12A", &[b'A', 0x80, b' ']).unwrap_err(),
            PdbResidueKeyError::ResidueNameNotAscii { byte: 0x80 }
        );
    }

    #[test]
    fn gemmi_pdb_segment_uses_strict_length_gate_and_exact_four_byte_field() {
        let residue = read_res_id(b"  12A", b"GLY").unwrap();
        for (length, expected) in [
            (72, ""),
            (73, "X"),
            (74, "XY"),
            (75, "XYZ"),
            (76, "XYZW"),
            (77, "XYZW"),
        ] {
            let mut line = vec![b' '; length];
            let segment_bytes = expected.as_bytes();
            let copied = segment_bytes.len().min(length.saturating_sub(72));
            if copied != 0 {
                line[72..72 + copied].copy_from_slice(&segment_bytes[..copied]);
            }
            if length == 77 {
                line[76] = b'Q'; // outside the source's four-byte field
            }
            let actual = read_pdb_residue_segment(residue, &line).unwrap();
            assert_eq!(actual.segment(), expected, "line length {length}");
            assert_eq!(actual.sequence_number(), Some(12));
            assert_eq!(actual.insertion_code(), Some(b'A'));
            assert_eq!(actual.name().as_str(), "GLY");
        }

        let mut middle_space = vec![b' '; 76];
        middle_space[72..76].copy_from_slice(b"A B ");
        assert_eq!(
            read_pdb_residue_segment(residue, &middle_space)
                .unwrap()
                .segment(),
            "A B"
        );

        let mut non_ascii = vec![b' '; 73];
        non_ascii[72] = 0x80;
        assert_eq!(
            read_pdb_residue_segment(residue, &non_ascii).unwrap_err(),
            PdbResidueKeyError::SegmentNotAscii { byte: 0x80 }
        );
    }

    #[test]
    fn pdb_line_cursor_normalizes_gemmi_max_line_options_and_truncates_at_120() {
        for max_line_length in [0, -4, 121] {
            let options = PdbInputOptions {
                max_line_length,
                ..PdbInputOptions::default()
            };
            assert_eq!(options.effective_max_line_length(), 120);

            let input = format!("{}\nNEXT\n", "x".repeat(121));
            let mut cursor = PdbLineCursor::new(&input, options);
            assert_eq!(cursor.copy_line(), Some(&b"x".repeat(120)[..]));
            assert_eq!(cursor.copy_line(), Some(&b"NEXT\n"[..]));
            assert_eq!(cursor.copy_line(), None);
        }

        assert_eq!(PdbInputOptions::default().effective_max_line_length(), 120);
        assert_eq!(
            PdbInputOptions {
                max_line_length: 1,
                ..PdbInputOptions::default()
            }
            .effective_max_line_length(),
            1
        );
        assert_eq!(
            PdbInputOptions {
                max_line_length: 120,
                ..PdbInputOptions::default()
            }
            .effective_max_line_length(),
            120
        );
    }

    #[test]
    fn gemmi_pdb_reader_initializes_canonical_state_and_tracks_first_non_ascii_line() {
        let mut state = PdbReaderState::new(
            r"C:\fixtures\sample.pdb.gz",
            true,
            BioStructureSourceState::default(),
        );
        assert_eq!(state.input_format, BioCoordinateFormat::Pdb);
        assert_eq!(state.source_state.name, "sample");
        assert_eq!(state.line_number, 0);
        assert_eq!(state.source_state.non_ascii_line, 0);

        // The pinned basename strips suffixes in the supplied order, only
        // when a nonempty stem remains, and case-sensitively.
        for (source, expected) in [
            ("/fixtures/.pdb", ".pdb"),
            ("/fixtures/name.pdb.gz", "name"),
            ("/fixtures/NAME.PDB", "NAME.PDB"),
            ("/fixtures/name.pdb.gz.extra", "name.pdb.gz.extra"),
        ] {
            assert_eq!(
                PdbReaderState::new(source, false, BioStructureSourceState::default())
                    .source_state
                    .name,
                expected,
                "source {source:?}"
            );
        }

        let mut prior = BioStructureSourceState::default();
        prior.name = "old".to_owned();
        prior.non_ascii_line = 17;
        let mut prepopulated = PdbReaderState::new("new.pdb", true, prior);
        assert_eq!(prepopulated.source_state.name, "new");
        prepopulated.record_line(b"ASCII");
        prepopulated.record_line(&[b'R', 0x80]);
        assert_eq!(prepopulated.line_number, 2);
        assert_eq!(prepopulated.source_state.non_ascii_line, 17);

        let mut cursor = PdbLineCursor::new(
            "REMARK ASCII\nREMARK \u{00e9}\nREMARK \u{2603}\n",
            PdbInputOptions::default(),
        );
        let mut checked =
            PdbReaderState::new("nonascii.pdb", true, BioStructureSourceState::default());
        while let Some(line) = cursor.copy_line() {
            checked.record_line(line);
        }
        assert_eq!(checked.line_number, 3);
        assert_eq!(checked.source_state.non_ascii_line, 2);

        let mut unchecked =
            PdbReaderState::new("nonascii.pdb", false, BioStructureSourceState::default());
        unchecked.record_line("REMARK \u{00e9}".as_bytes());
        assert_eq!(unchecked.line_number, 1);
        assert_eq!(unchecked.source_state.non_ascii_line, 0);
    }

    #[test]
    fn gemmi_pdb_seqres_trims_entity_names_and_appends_nonempty_slots_in_source_order() {
        let mut state =
            PdbReaderState::new("seqres.pdb", false, BioStructureSourceState::default());

        let first = seqres_source_line(*b"sEqReS", *b"A ", &[(19, *b"ALA"), (27, *b"GLY")], 68);
        assert!(state.seqres_record(&first));

        // Both fixed-width spellings trim to the same source chain name. The
        // blank slot at 23 is skipped; it does not stop slot 27 or later rows.
        let second = seqres_source_line(*b"SEQRES", *b" A", &[(19, *b"SER")], 68);
        assert!(state.seqres_record(&second));

        // Entity identity remains byte/case-sensitive even though record
        // matching is case-insensitive.
        let lower_case_chain = seqres_source_line(*b"SEQRES", *b"a ", &[(19, *b"VAL")], 68);
        assert!(state.seqres_record(&lower_case_chain));

        let other_record = seqres_source_line(*b"HETATM", *b"Z ", &[(19, *b"THR")], 68);
        assert!(!state.seqres_record(&other_record));

        let entities = &state.entities.entities;
        assert_eq!(entities.len(), 2);
        assert_eq!(&entities[0].source_name[..], b"A");
        assert_eq!(entities[0].entity_kind, EntityKind::Polymer);
        assert_eq!(
            entities[0].full_sequence,
            vec![b"ALA".to_vec(), b"GLY".to_vec(), b"SER".to_vec()]
        );
        assert_eq!(&entities[1].source_name[..], b"a");
        assert_eq!(entities[1].entity_kind, EntityKind::Polymer);
        assert_eq!(entities[1].full_sequence, vec![b"VAL".to_vec()]);
    }

    #[test]
    fn gemmi_pdb_seqres_creates_empty_entities_and_promotes_existing_kind() {
        let mut state =
            PdbReaderState::new("seqres.pdb", false, BioStructureSourceState::default());
        state.entities.entities.push(PdbEntity {
            source_name: b"Q".to_vec(),
            entity_kind: EntityKind::NonPolymer,
            full_sequence: vec![b"OLD".to_vec()],
            dbrefs: Vec::new(),
            subchains: Vec::new(),
        });

        let promotion = seqres_source_line(*b"SEQRES", *b"Q ", &[(19, *b"THR")], 68);
        assert!(state.seqres_record(&promotion));

        // The source creates/promotes an entity even when no sequence slot is
        // present because it finds/adds and assigns Polymer before the loop.
        let empty = seqres_source_line(*b"SEQRES", *b"C ", &[], 19);
        assert!(state.seqres_record(&empty));

        let entities = &state.entities.entities;
        assert_eq!(entities.len(), 2);
        assert_eq!(&entities[0].source_name[..], b"Q");
        assert_eq!(entities[0].entity_kind, EntityKind::Polymer);
        assert_eq!(
            entities[0].full_sequence,
            vec![b"OLD".to_vec(), b"THR".to_vec()]
        );
        assert_eq!(&entities[1].source_name[..], b"C");
        assert_eq!(entities[1].entity_kind, EntityKind::Polymer);
        assert!(entities[1].full_sequence.is_empty());
    }

    #[test]
    fn gemmi_pdb_dbref_populates_blank_record_fields_and_ranges() {
        let mut state = PdbReaderState::new("dbref.pdb", false, BioStructureSourceState::default());
        let mut line = dbref_source_line(b' ', *b"A ");
        set_dbref_source_field(&mut line, 14, 5, b"  12A");
        set_dbref_source_field(&mut line, 20, 5, b"  20B");
        set_dbref_source_field(&mut line, 26, 6, b"UNP");
        set_dbref_source_field(&mut line, 33, 8, b"P12345");
        set_dbref_source_field(&mut line, 42, 12, b"PROT_HUMAN");
        set_dbref_source_field(&mut line, 55, 5, b"  101");
        line[60] = b'C';
        set_dbref_source_field(&mut line, 62, 5, b"  202");
        line[67] = b'D';

        assert_eq!(state.dbref_record(&line), Ok(PdbDbRefAction::Continue));
        assert_eq!(state.entities.entities.len(), 1);
        let entity = &state.entities.entities[0];
        assert_eq!(entity.source_name, b"A");
        assert_eq!(entity.entity_kind, EntityKind::Polymer);
        assert_eq!(entity.dbrefs.len(), 1);
        let dbref = &entity.dbrefs[0];
        assert_eq!(dbref.db_name, b"UNP");
        assert_eq!(dbref.accession_code, b"P12345");
        assert_eq!(dbref.id_code, b"PROT_HUMAN");
        assert!(dbref.isoform.is_empty());
        assert_eq!(dbref.seq_begin, PdbSeqId::new(12, Some(b'A')));
        assert_eq!(dbref.seq_end, PdbSeqId::new(20, Some(b'B')));
        assert_eq!(dbref.db_begin, PdbSeqId::new(101, Some(b'C')));
        assert_eq!(dbref.db_end, PdbSeqId::new(202, Some(b'D')));
        assert_eq!(dbref.label_seq_begin, None);
        assert_eq!(dbref.label_seq_end, None);
    }

    #[test]
    fn gemmi_pdb_dbref1_dbref2_pair_preserves_and_completes_the_same_row() {
        let mut state =
            PdbReaderState::new("dbref-pair.pdb", false, BioStructureSourceState::default());
        let mut first = dbref_source_line(b'1', *b"B ");
        set_dbref_source_field(&mut first, 14, 5, b"   1 ");
        set_dbref_source_field(&mut first, 20, 5, b"  99 ");
        set_dbref_source_field(&mut first, 26, 6, b"UNP");
        set_dbref_source_field(&mut first, 47, 20, b"Q98765_HUMAN");
        assert_eq!(state.dbref_record(&first), Ok(PdbDbRefAction::Continue));

        let mut continuation = dbref_source_line(b'2', *b"B ");
        set_dbref_source_field(&mut continuation, 18, 22, b"Q98765");
        set_dbref_source_field(&mut continuation, 45, 10, b"        10");
        set_dbref_source_field(&mut continuation, 57, 10, b"       100");
        assert_eq!(
            state.dbref_record(&continuation),
            Ok(PdbDbRefAction::Continue)
        );

        let entity = &state.entities.entities[0];
        assert_eq!(entity.source_name, b"B");
        assert_eq!(entity.entity_kind, EntityKind::Polymer);
        assert_eq!(entity.dbrefs.len(), 1);
        let dbref = &entity.dbrefs[0];
        assert_eq!(dbref.db_name, b"UNP");
        assert_eq!(dbref.accession_code, b"Q98765");
        assert_eq!(dbref.id_code, b"Q98765_HUMAN");
        assert_eq!(dbref.seq_begin, PdbSeqId::new(1, None));
        assert_eq!(dbref.seq_end, PdbSeqId::new(99, None));
        assert_eq!(dbref.db_begin, PdbSeqId::new(10, None));
        assert_eq!(dbref.db_end, PdbSeqId::new(100, None));
    }

    #[test]
    fn gemmi_pdb_dbref_rows_append_in_order_and_unknown_suffix_does_not_rewrite() {
        let mut state = PdbReaderState::new(
            "dbref-multiple.pdb",
            false,
            BioStructureSourceState::default(),
        );
        let mut first = dbref_source_line(b' ', *b"A ");
        set_dbref_source_field(&mut first, 14, 5, b"   1 ");
        assert_eq!(state.dbref_record(&first), Ok(PdbDbRefAction::Continue));

        let mut second = dbref_source_line(b' ', *b"A ");
        set_dbref_source_field(&mut second, 14, 5, b"   2 ");
        assert_eq!(state.dbref_record(&second), Ok(PdbDbRefAction::Continue));
        assert_eq!(state.entities.entities[0].dbrefs.len(), 2);
        assert_eq!(
            state.entities.entities[0].dbrefs[0].seq_begin,
            PdbSeqId::new(1, None)
        );
        assert_eq!(
            state.entities.entities[0].dbrefs[1].seq_begin,
            PdbSeqId::new(2, None)
        );

        let before = state.entities.entities[0].dbrefs.clone();
        let unknown = dbref_source_line(b'3', *b"A ");
        assert_eq!(state.dbref_record(&unknown), Ok(PdbDbRefAction::Continue));
        assert_eq!(state.entities.entities[0].dbrefs, before);
    }

    #[test]
    fn gemmi_pdb_dbref2_without_predecessor_stops_after_entity_promotion() {
        let mut state = PdbReaderState::new(
            "dbref-orphan.pdb",
            false,
            BioStructureSourceState::default(),
        );
        let orphan = dbref_source_line(b'2', *b"C ");
        assert_eq!(state.dbref_record(&orphan), Ok(PdbDbRefAction::Stop));
        assert_eq!(state.entities.entities.len(), 1);
        assert_eq!(state.entities.entities[0].source_name, b"C");
        assert_eq!(state.entities.entities[0].entity_kind, EntityKind::Polymer);
        assert!(state.entities.entities[0].dbrefs.is_empty());

        let mut unknown_state = PdbReaderState::new(
            "dbref-unknown.pdb",
            false,
            BioStructureSourceState::default(),
        );
        let unknown_without_predecessor = dbref_source_line(b'3', *b"  ");
        assert_eq!(
            unknown_state.dbref_record(&unknown_without_predecessor),
            Ok(PdbDbRefAction::Stop)
        );
        assert_eq!(unknown_state.entities.entities.len(), 1);
        assert!(unknown_state.entities.entities[0].source_name.is_empty());
        assert_eq!(
            unknown_state.entities.entities[0].entity_kind,
            EntityKind::Polymer
        );
    }

    #[test]
    fn gemmi_pdb_modres_retains_canonical_fields_and_record_order() {
        let mut state =
            PdbReaderState::new("modres.pdb", false, BioStructureSourceState::default());
        let first = modres_source_line(
            *b"MODRES",
            *b" A",
            *b"  12B",
            *b"MSE",
            *b"MET",
            b"SELENOMETHIONINE",
            b"REFMAC1",
            80,
            *b"  ",
        );
        let second = modres_source_line(
            *b"modres", *b" B", *b"   7C", *b"TPO", *b"THR", b"D", b"RFX", 30, *b"  ",
        );

        assert_eq!(state.modres_record(&first, 80), Ok(true));
        assert_eq!(state.modres_record(&second, 30), Ok(true));
        assert_eq!(state.modres_record(&first, 80), Ok(true));
        assert_eq!(
            state.mod_residues,
            vec![
                BioModRes {
                    chain_name: PdbChainId::from_ascii(b"A").unwrap(),
                    res_id: read_res_id(b"  12B", b"MSE").unwrap(),
                    parent_comp_id: "MET".to_owned(),
                    mod_id: "REFMAC1".to_owned(),
                    details: "SELENOMETHIONINE".to_owned(),
                },
                BioModRes {
                    chain_name: PdbChainId::from_ascii(b"B").unwrap(),
                    res_id: read_res_id(b"   7C", b"TPO").unwrap(),
                    parent_comp_id: "THR".to_owned(),
                    mod_id: String::new(),
                    details: "D".to_owned(),
                },
                BioModRes {
                    chain_name: PdbChainId::from_ascii(b"A").unwrap(),
                    res_id: read_res_id(b"  12B", b"MSE").unwrap(),
                    parent_comp_id: "MET".to_owned(),
                    mod_id: "REFMAC1".to_owned(),
                    details: "SELENOMETHIONINE".to_owned(),
                },
            ]
        );
    }

    #[test]
    fn gemmi_pdb_modres_applies_exact_details_and_refmac_gates() {
        let cases = [
            (29, *b"  ", "", ""),
            (30, *b"  ", "D", ""),
            (72, *b"  ", "DETAIL", ""),
            (73, *b"  ", "DETAIL", "R"),
            (80, *b"X ", "DETAIL", ""),
            (80, *b" X", "DETAIL", ""),
        ];
        let mut state = PdbReaderState::new(
            "modres-gates.pdb",
            false,
            BioStructureSourceState::default(),
        );

        for (line_len, guard, expected_details, expected_mod_id) in cases {
            let line = modres_source_line(
                *b"MODRES", *b" A", *b"  12B", *b"MSE", *b"MET", b"DETAIL", b"REFMAC1", line_len,
                guard,
            );
            assert_eq!(state.modres_record(&line, line_len), Ok(true));
            let record = state.mod_residues.last().unwrap();
            assert_eq!(record.details, expected_details, "line length {line_len}");
            assert_eq!(record.mod_id, expected_mod_id, "line length {line_len}");
        }
        assert_eq!(state.mod_residues.len(), cases.len());
    }

    #[test]
    fn gemmi_pdb_modres_reports_unrepresentable_fields_without_partial_append() {
        let mut state = PdbReaderState::new(
            "modres-invalid.pdb",
            false,
            BioStructureSourceState::default(),
        );
        let valid = |chain, residue, parent: [u8; 3], details: &[u8], mod_id: &[u8]| {
            modres_source_line(
                *b"MODRES", chain, *b"  12B", residue, parent, details, mod_id, 80, *b"  ",
            )
        };

        let bad_chain = valid(*b" \x80", *b"MSE", *b"MET", b"", b"");
        assert_eq!(
            state.modres_record(&bad_chain, 80),
            Err(PdbModResError::ChainNameNotRepresentable {
                width: 1,
                non_ascii: Some(0x80),
            })
        );
        let bad_residue = valid(*b" A", [0x80, b'S', b'E'], *b"MET", b"", b"");
        assert_eq!(
            state.modres_record(&bad_residue, 80),
            Err(PdbModResError::ResidueAddress(
                PdbResidueKeyError::ResidueNameNotAscii { byte: 0x80 }
            ))
        );
        let bad_parent = valid(*b" A", *b"MSE", [0xC3, b' ', b' '], b"", b"");
        assert_eq!(
            state.modres_record(&bad_parent, 80),
            Err(PdbModResError::TextFieldNotUtf8 {
                field_offset: 24,
                valid_up_to: 0,
            })
        );
        let bad_details = valid(*b" A", *b"MSE", *b"MET", &[0xC3], b"");
        assert_eq!(
            state.modres_record(&bad_details, 80),
            Err(PdbModResError::TextFieldNotUtf8 {
                field_offset: 29,
                valid_up_to: 0,
            })
        );
        let bad_mod_id = valid(*b" A", *b"MSE", *b"MET", b"", &[0xC3]);
        assert_eq!(
            state.modres_record(&bad_mod_id, 80),
            Err(PdbModResError::TextFieldNotUtf8 {
                field_offset: 72,
                valid_up_to: 0,
            })
        );
        assert!(state.mod_residues.is_empty());

        let other_record = [0; 122];
        assert_eq!(state.modres_record(&other_record, 80), Ok(false));
        assert!(state.mod_residues.is_empty());
    }

    #[test]
    fn gemmi_pdb_helix_short_lines_and_record_dispatch_match_source() {
        let mut state =
            PdbReaderState::new("helix-short.pdb", false, BioStructureSourceState::default());
        let short = helix_source_line(
            *b"HELIX ", *b" A", *b"  12B", *b"ALA", *b" B", *b"  34C", *b"GLY", *b" 5", *b"   17",
            39,
        );

        assert_eq!(state.helix_record(&short, 39), Ok(true));
        assert!(state.helices.is_empty());

        let other_record = helix_source_line(
            *b"OTHER ", *b" A", *b"  12B", *b"ALA", *b" B", *b"  34C", *b"GLY", *b" 5", *b"   17",
            80,
        );
        assert_eq!(state.helix_record(&other_record, 80), Ok(false));
        assert!(state.helices.is_empty());
    }

    #[test]
    fn gemmi_pdb_helix_endpoints_class_and_length_match_pinned_oracle() {
        let cases = [
            (40, *b"HELIX ", *b" 5", *b"   17", BioHelixClass::R310, -1),
            (72, *b"HELIX ", *b" 5", *b"   17", BioHelixClass::R310, -1),
            (73, *b"HELIX ", *b" 5", *b"7    ", BioHelixClass::R310, 7),
            (77, *b"HELIX ", *b" 5", *b"   17", BioHelixClass::R310, 17),
            (
                80,
                *b"HELIX ",
                *b"11",
                *b"   17",
                BioHelixClass::UnknownHelix,
                17,
            ),
            (
                80,
                *b"HELIX ",
                *b" X",
                *b"  X  ",
                BioHelixClass::UnknownHelix,
                0,
            ),
            (80, *b"helix ", *b" 1", *b"   17", BioHelixClass::RAlpha, 17),
        ];
        let mut state = PdbReaderState::new(
            "helix-branches.pdb",
            false,
            BioStructureSourceState::default(),
        );

        for (line_len, record, class, length_field, expected_class, expected_length) in cases {
            let line = helix_source_line(
                record,
                *b" A",
                *b"  12B",
                *b"ALA",
                *b" B",
                *b"  34C",
                *b"GLY",
                class,
                length_field,
                line_len,
            );
            assert_eq!(state.helix_record(&line, line_len), Ok(true));
            let helix = state.helices.last().unwrap();
            assert_eq!(
                helix.start,
                expected_helix_address(b"A", 12, b'B', b"ALA"),
                "start endpoint at physical line length {line_len}"
            );
            assert_eq!(
                helix.end,
                expected_helix_address(b"B", 34, b'C', b"GLY"),
                "end endpoint at physical line length {line_len}"
            );
            assert_eq!(helix.pdb_helix_class, expected_class);
            assert_eq!(helix.length, expected_length);
        }
        assert_eq!(state.helices.len(), cases.len());
    }

    #[test]
    fn gemmi_pdb_helix_unrepresentable_chain_is_typed_and_atomic() {
        let mut state = PdbReaderState::new(
            "helix-invalid.pdb",
            false,
            BioStructureSourceState::default(),
        );
        let invalid = helix_source_line(
            *b"HELIX ",
            [0x80, b' '],
            *b"  12B",
            *b"ALA",
            *b" B",
            *b"  34C",
            *b"GLY",
            *b" 5",
            *b"   17",
            80,
        );

        assert_eq!(
            state.helix_record(&invalid, 80),
            Err(PdbHelixError::ChainNameNotRepresentable {
                field_offset: 18,
                width: 1,
                non_ascii: Some(0x80),
            })
        );
        assert!(state.helices.is_empty());
    }

    #[test]
    fn gemmi_pdb_sheet_short_records_and_nonmatching_prefix_preserve_state() {
        let mut state = PdbReaderState::new(
            "sheet-dispatch.pdb",
            false,
            BioStructureSourceState::default(),
        );
        let short = sheet_source_line(*b"SHEET ", 39, &[]);
        assert_eq!(state.sheet_record(&short, 39), Ok(true));
        assert!(state.sheets.is_empty());

        let other = sheet_source_line(*b"OTHER ", 80, &[]);
        assert_eq!(state.sheet_record(&other, 80), Ok(false));
        assert!(state.sheets.is_empty());
    }

    #[test]
    fn gemmi_pdb_sheet_groups_first_seen_ids_and_preserves_strand_addresses() {
        let mut state =
            PdbReaderState::new("sheet-order.pdb", false, BioStructureSourceState::default());
        let first_aa = sheet_source_line(
            *b"SHEET ",
            67,
            &[
                (11, &b"AA "[..]),
                (17, &b"ALA"[..]),
                (20, &b" A"[..]),
                (22, &b"  12B"[..]),
                (28, &b"GLY"[..]),
                (31, &b" B"[..]),
                (33, &b"  20C"[..]),
                (38, &b" 0"[..]),
            ],
        );
        let first_bb = sheet_source_line(
            *b"SHEET ",
            80,
            &[
                (11, &b"BB "[..]),
                (17, &b"SER"[..]),
                (20, &b" C"[..]),
                (22, &b"  30 "[..]),
                (28, &b"THR"[..]),
                (31, &b" D"[..]),
                (33, &b"  40 "[..]),
                (38, &b" 1"[..]),
                (41, &b" O  "[..]),
                (45, &b"GLY"[..]),
                (48, &b" E"[..]),
                (50, &b"  33F"[..]),
                (56, &b" N  "[..]),
                (60, &b"ASN"[..]),
                (63, &b" F"[..]),
                (65, &b"  37G"[..]),
            ],
        );
        let second_aa = sheet_source_line(
            *b"SHEET ",
            80,
            &[
                (11, &b"AA "[..]),
                (17, &b"LYS"[..]),
                (20, &b" G"[..]),
                (22, &b"  50H"[..]),
                (28, &b"ASP"[..]),
                (31, &b" H"[..]),
                (33, &b"  60I"[..]),
                (38, &b"-1"[..]),
                (41, &b" O  "[..]),
                (45, &b"GLU"[..]),
                (48, &b" I"[..]),
                (50, &b"  55J"[..]),
                (56, &b" N  "[..]),
                (60, &b"VAL"[..]),
                (63, &b" J"[..]),
                (65, &b"  59K"[..]),
            ],
        );

        assert_eq!(state.sheet_record(&first_aa, 67), Ok(true));
        assert_eq!(state.sheet_record(&first_bb, 80), Ok(true));
        assert_eq!(state.sheet_record(&second_aa, 80), Ok(true));

        assert_eq!(
            state
                .sheets
                .iter()
                .map(|sheet| sheet.name.as_str())
                .collect::<Vec<_>>(),
            ["AA", "BB"]
        );
        assert_eq!(state.sheets[0].strands.len(), 2);
        assert_eq!(state.sheets[1].strands.len(), 1);

        let aa_first = &state.sheets[0].strands[0];
        assert_eq!(aa_first.sense, 0);
        assert_eq!(
            aa_first.start,
            expected_sheet_address(b"A", 12, b'B', b"ALA", "")
        );
        assert_eq!(
            aa_first.end,
            expected_sheet_address(b"B", 20, b'C', b"GLY", "")
        );
        assert_eq!(aa_first.hbond_atom2, AtomAddress::default());
        assert_eq!(aa_first.hbond_atom1, AtomAddress::default());
        assert!(aa_first.name.is_empty());

        let aa_second = &state.sheets[0].strands[1];
        assert_eq!(aa_second.sense, -1);
        assert_eq!(
            aa_second.start,
            expected_sheet_address(b"G", 50, b'H', b"LYS", "")
        );
        assert_eq!(
            aa_second.end,
            expected_sheet_address(b"H", 60, b'I', b"ASP", "")
        );
        assert_eq!(
            aa_second.hbond_atom2,
            expected_sheet_address(b"I", 55, b'J', b"GLU", "O")
        );
        assert_eq!(
            aa_second.hbond_atom1,
            expected_sheet_address(b"J", 59, b'K', b"VAL", "N")
        );

        let bb = &state.sheets[1].strands[0];
        assert_eq!(bb.sense, 1);
        assert_eq!(bb.start, expected_sheet_address(b"C", 30, b' ', b"SER", ""));
        assert_eq!(bb.end, expected_sheet_address(b"D", 40, b' ', b"THR", ""));
        assert_eq!(
            bb.hbond_atom2,
            expected_sheet_address(b"E", 33, b'F', b"GLY", "O")
        );
        assert_eq!(
            bb.hbond_atom1,
            expected_sheet_address(b"F", 37, b'G', b"ASN", "N")
        );
    }

    #[test]
    fn gemmi_pdb_sheet_hbond_gate_uses_source_length_and_reused_suffix() {
        let mut state = PdbReaderState::new(
            "sheet-reused-suffix.pdb",
            false,
            BioStructureSourceState::default(),
        );
        let previous_full_row = sheet_source_line(
            *b"SHEET ",
            80,
            &[
                (11, &b"OLD"[..]),
                (17, &b"SER"[..]),
                (20, &b" C"[..]),
                (22, &b"  30 "[..]),
                (28, &b"THR"[..]),
                (31, &b" D"[..]),
                (33, &b"  40 "[..]),
                (38, &b" 1"[..]),
                (41, &b" O  "[..]),
                (45, &b"GLY"[..]),
                (48, &b" E"[..]),
                (50, &b"  33F"[..]),
                (56, &b" N  "[..]),
                (60, &b"ASN"[..]),
                (63, &b" F"[..]),
                (65, &b"  37G"[..]),
            ],
        );
        assert_eq!(state.sheet_record(&previous_full_row, 80), Ok(true));

        let current_full_row = sheet_source_line(
            *b"SHEET ",
            80,
            &[
                (11, &b"EDG"[..]),
                (17, &b"LYS"[..]),
                (20, &b" G"[..]),
                (22, &b"  50H"[..]),
                (28, &b"ASP"[..]),
                (31, &b" H"[..]),
                (33, &b"  60I"[..]),
                (38, &b"-1"[..]),
                (41, &b" O  "[..]),
                (45, &b"GLU"[..]),
                (48, &b" I"[..]),
                (50, &b"  55J"[..]),
                (56, &b" N  "[..]),
                (60, &b"VAL"[..]),
                (63, &b" J"[..]),
                (65, &b"  59K"[..]),
            ],
        );
        let mut reused_line_buffer = previous_full_row;
        reused_line_buffer[..67].copy_from_slice(&current_full_row[..67]);
        reused_line_buffer[67] = b'\n';
        reused_line_buffer[68] = 0;

        assert_eq!(state.sheet_record(&reused_line_buffer, 68), Ok(true));
        assert_eq!(state.sheets[1].name, "EDG");
        let edge = &state.sheets[1].strands[0];
        assert_eq!(edge.sense, -1);
        assert_eq!(
            edge.hbond_atom2,
            expected_sheet_address(b"I", 55, b'J', b"GLU", "O")
        );
        assert_eq!(
            edge.hbond_atom1,
            expected_sheet_address(b"J", 0, b'G', b"VAL", "N")
        );
    }

    #[test]
    fn gemmi_pdb_hetnam_capture_matches_record_length_column_and_trim_gates() {
        let mut state = PdbReaderState::new(
            "hetnam-gates.pdb",
            false,
            BioStructureSourceState::default(),
        );

        // Pinned Gemmi oracle: lower-case record prefix with a trimmed short
        // field and full code is accepted; `len == 71` does not enter the arm;
        // at length 72, read_string sees the one available `B` before NUL.
        let mixed_case = hetnam_source_line(*b"hEtNAM", *b" A ", b" BBB    ", 80, b' ');
        assert!(state.hetnam_record(&mixed_case, 80));
        assert_eq!(
            state.shortened_ccd_codes,
            vec![PdbCcdAlias {
                full_code: b"BBB".to_vec(),
                short_code: b"A".to_vec(),
            }]
        );

        let length_71 = hetnam_source_line(*b"HETNAM", *b"AAA", b"BBBB", 71, b' ');
        assert!(state.hetnam_record(&length_71, 71));
        assert_eq!(state.shortened_ccd_codes.len(), 1);

        let length_72 = hetnam_source_line(*b"HETNAM", *b"AAA", b"BBBB", 72, b' ');
        assert!(state.hetnam_record(&length_72, 72));
        assert_eq!(
            state.shortened_ccd_codes.last(),
            Some(&PdbCcdAlias {
                full_code: b"B".to_vec(),
                short_code: b"AAA".to_vec(),
            })
        );

        let wrong_column = hetnam_source_line(*b"HETNAM", *b"AAA", b"CCCC", 80, b'X');
        assert!(state.hetnam_record(&wrong_column, 80));
        let blank_full_code = hetnam_source_line(*b"HETNAM", *b"AAA", b"        ", 80, b' ');
        assert!(state.hetnam_record(&blank_full_code, 80));
        assert_eq!(state.shortened_ccd_codes.len(), 2);

        let other_record = hetnam_source_line(*b"HETX  ", *b"AAA", b"DDDD", 80, b' ');
        assert!(!state.hetnam_record(&other_record, 80));
        assert_eq!(state.shortened_ccd_codes.len(), 2);
    }

    #[test]
    fn gemmi_pdb_hetnam_restores_modres_sequences_and_residues_in_alias_order() {
        let mut state = PdbReaderState::new(
            "hetnam-restore.pdb",
            false,
            BioStructureSourceState::default(),
        );
        let first_alias = hetnam_source_line(*b"HETNAM", *b"AAA", b"BBB", 80, b' ');
        let second_alias = hetnam_source_line(*b"HETNAM", *b"BBB", b"CCC", 80, b' ');
        assert!(state.hetnam_record(&first_alias, 80));
        assert!(state.hetnam_record(&second_alias, 80));

        let sequence = seqres_source_line(
            *b"SEQRES",
            *b"A ",
            &[(19, *b"AAA"), (23, *b"BBB"), (27, *b"AAA")],
            30,
        );
        assert!(state.seqres_record(&sequence));

        let modres = modres_source_line(
            *b"MODRES", *b" A", *b"   1 ", *b"AAA", *b"ALA", b"DETAIL", b"", 80, *b"  ",
        );
        assert_eq!(state.modres_record(&modres, 80), Ok(true));

        let mut grouping = PdbHierarchyGrouping::default();
        let model_id = grouping.add_model(Some(1)).unwrap();
        for (serial, name, sequence_number) in [(1, *b"AAA", 1), (2, *b"XXX", 2)] {
            let line = grouped_pdb_atom_line(
                *b"ATOM  ",
                serial,
                *b" CA ",
                name,
                b'A',
                sequence_number,
                b' ',
                *b"SEG1",
            );
            add_grouped_pdb_atom(&mut grouping, model_id, &line);
        }

        state
            .restore_full_ccd_codes(&mut grouping, &mut [], &mut [])
            .unwrap();

        assert!(state.shortened_ccd_codes.is_empty());
        assert_eq!(
            state.mod_residues[0].res_id.name(),
            ResidueName::from_ascii(b"CCC").unwrap()
        );
        assert_eq!(
            state.entities.entities[0].full_sequence,
            [b"CCC".to_vec(), b"CCC".to_vec(), b"CCC".to_vec()]
        );
        let residues = &grouping.models[model_id.index()].chains[0].residues;
        assert_eq!(
            residues[0].address.name(),
            ResidueName::from_ascii(b"CCC").unwrap()
        );
        assert_eq!(residues[0].address.sequence_number(), Some(1));
        assert_eq!(residues[0].address.segment().as_bytes(), b"SEG1");
        assert_eq!(
            residues[1].address.name(),
            ResidueName::from_ascii(b"XXX").unwrap()
        );
        assert_eq!(residues[1].address.sequence_number(), Some(2));
    }

    #[test]
    fn gemmi_pdb_hetnam_overwide_restore_is_typed_and_failure_atomic() {
        let mut state = PdbReaderState::new(
            "hetnam-width.pdb",
            false,
            BioStructureSourceState::default(),
        );
        let alias = hetnam_source_line(*b"HETNAM", *b"AAA", b"ABCDE", 80, b' ');
        assert!(state.hetnam_record(&alias, 80));
        let sequence = seqres_source_line(*b"SEQRES", *b"A ", &[(19, *b"AAA")], 22);
        assert!(state.seqres_record(&sequence));

        let mut grouping = PdbHierarchyGrouping::default();
        let model_id = grouping.add_model(Some(1)).unwrap();
        let line = grouped_pdb_atom_line(*b"ATOM  ", 1, *b" CA ", *b"AAA", b'A', 1, b' ', *b"SEG1");
        add_grouped_pdb_atom(&mut grouping, model_id, &line);
        let state_before = state.clone();
        let grouping_before = grouping.clone();

        assert_eq!(
            state.restore_full_ccd_codes(&mut grouping, &mut [], &mut []),
            Err(PdbResidueKeyError::ResidueNameTooWide { width: 5 })
        );
        assert_eq!(state, state_before);
        assert_eq!(grouping, grouping_before);
    }

    #[test]
    fn gemmi_pdb_conect_accumulates_ordered_targets_without_lowering_to_bonds() {
        use std::collections::BTreeMap;

        let mut state =
            PdbReaderState::new("conect.pdb", false, BioStructureSourceState::default());
        for (source, targets, line_len) in [
            (9, &[20, 0, 20, 3, 99][..], 60),
            (2, &[7, 0, 8, 0][..], 31),
            (9, &[3, 2, 0, 4][..], 31),
        ] {
            let line = conect_source_line(source, targets, line_len);
            assert_eq!(state.conect_record(&line, line_len), Ok(true));
        }

        assert_eq!(
            state.source_state.conect_map,
            BTreeMap::from([(2, vec![7, 8]), (9, vec![20, 20, 3, 3, 2, 4]),])
        );
        assert!(state.entities.entities.is_empty());
    }

    #[test]
    fn gemmi_pdb_conect_observes_short_zero_empty_and_record_prefix_gates() {
        use std::collections::BTreeMap;

        let mut state = PdbReaderState::new(
            "conect-gates.pdb",
            false,
            BioStructureSourceState::default(),
        );

        let short = conect_source_line(12, &[1], 10);
        assert_eq!(state.conect_record(&short, 10), Ok(true));
        assert!(state.source_state.conect_map.is_empty());

        let zero_source = conect_source_line(0, &[1, 2, 3, 4], 31);
        assert_eq!(state.conect_record(&zero_source, 31), Ok(true));
        assert!(state.source_state.conect_map.is_empty());

        let empty = conect_source_line(13, &[], 11);
        assert_eq!(state.conect_record(&empty, 11), Ok(true));
        assert_eq!(
            state.source_state.conect_map,
            BTreeMap::from([(13, Vec::new())])
        );

        let mut lowercase_prefix = conect_source_line(14, &[5], 16);
        lowercase_prefix[..6].copy_from_slice(b"coneCT");
        assert_eq!(state.conect_record(&lowercase_prefix, 16), Ok(true));

        let mut other_record = conect_source_line(15, &[6], 16);
        other_record[..6].copy_from_slice(b"ATOM  ");
        assert_eq!(state.conect_record(&other_record, 16), Ok(false));
        assert_eq!(
            state.source_state.conect_map,
            BTreeMap::from([(13, Vec::new()), (14, vec![5])])
        );
        assert!(state.entities.entities.is_empty());
    }

    #[test]
    fn gemmi_pdb_wrong_input_format_errors_preserve_prefix_state_and_source() {
        let cases = [
            (
                "DaTa_model\n",
                "folder/model.cif",
                PdbInputFormatError::Cif {
                    source: "folder/model.cif".to_owned(),
                },
                "Incorrect file format (perhaps it is cif not pdb?): folder/model.cif",
            ),
            (
                "{\"DaTa_\n",
                "folder/model.json",
                PdbInputFormatError::Mmjson {
                    source: "folder/model.json".to_owned(),
                },
                "Incorrect file format (perhaps it is mmJSON not pdb?): folder/model.json",
            ),
        ];

        for (input, source, expected, display) in cases {
            let mut cursor = PdbLineCursor::new(input, PdbInputOptions::default());
            let mut state = PdbReaderState::new(source, false, BioStructureSourceState::default());
            let line = cursor.copy_line().expect("input has a line").to_vec();
            state.record_line(&line);
            let error = state
                .wrong_input_format(cursor.source_line_buffer(), false, source)
                .expect("source format mismatch");
            assert_eq!(error, expected);
            assert_eq!(error.to_string(), display);
            assert_eq!(state.line_number, 1);
            assert!(
                state
                    .wrong_input_format(cursor.source_line_buffer(), true, source)
                    .is_none()
            );
        }

        let mut ordinary = PdbLineCursor::new("data\n", PdbInputOptions::default());
        let line = ordinary.copy_line().expect("input has a line").to_vec();
        let state = PdbReaderState::new("ordinary.pdb", false, BioStructureSourceState::default());
        assert!(
            state
                .wrong_input_format(ordinary.source_line_buffer(), false, "ordinary.pdb")
                .is_none()
        );
        assert_eq!(line, b"data\n");

        // Gemmi checks the reusable C line buffer, not only its strlen prefix.
        // A NUL after the complete mmJSON discriminator leaves the discriminator
        // readable and the source reports the format error; a NUL before `_`
        // does not complete the source's `ta_` check.
        for (input, expected_error) in [("{\"data_\0\nNEXT\n", true), ("{\"data\0_\nNEXT\n", false)]
        {
            let mut cursor = PdbLineCursor::new(input, PdbInputOptions::default());
            let mut state =
                PdbReaderState::new("embedded.json", false, BioStructureSourceState::default());
            let visible = cursor.copy_line().expect("visible prefix exists").to_vec();
            state.record_line(&visible);
            let error =
                state.wrong_input_format(cursor.source_line_buffer(), false, "embedded.json");
            assert_eq!(error.is_some(), expected_error, "input {input:?}");
            if expected_error {
                assert_eq!(
                    error.unwrap().to_string(),
                    "Incorrect file format (perhaps it is mmJSON not pdb?): embedded.json"
                );
            }
            assert_eq!(cursor.copy_line(), None, "copy_line drains the next line");
        }
    }

    #[test]
    fn gemmi_pdb_model_transition_matches_model_and_implicit_atom_rules() {
        // Pinned Gemmi 5cc1c23 oracle in IO-bio_pdb.md Step 252: repeated
        // empty MODEL number is reused; populated duplicates and MODEL before
        // ENDMDL fail; an implicit model uses model_count + 1 and fails only
        // when that exact number already exists.
        let mut empty_duplicate = PdbModelTransition::default();
        let first = empty_duplicate
            .model_record(&pdb_model_source_line(1))
            .unwrap();
        empty_duplicate.end_model_record();
        let reused = empty_duplicate
            .model_record(&pdb_model_source_line(1))
            .unwrap();
        assert_eq!(first, reused);
        assert_eq!(empty_duplicate.grouping.models.len(), 1);

        let mut populated_duplicate = PdbModelTransition::default();
        populated_duplicate
            .model_record(&pdb_model_source_line(1))
            .unwrap();
        add_transition_pdb_atom(&mut populated_duplicate).unwrap();
        populated_duplicate.end_model_record();
        assert!(matches!(
            populated_duplicate.model_record(&pdb_model_source_line(1)),
            Err(PdbModelTransitionError::DuplicateModelNumber { number: 1 })
        ));

        let mut missing_endmdl = PdbModelTransition::default();
        missing_endmdl
            .model_record(&pdb_model_source_line(1))
            .unwrap();
        add_transition_pdb_atom(&mut missing_endmdl).unwrap();
        assert!(matches!(
            missing_endmdl.model_record(&pdb_model_source_line(2)),
            Err(PdbModelTransitionError::ModelWithoutEnd)
        ));

        let mut implicit_free = PdbModelTransition::default();
        let location = add_transition_pdb_atom(&mut implicit_free).unwrap();
        assert_eq!(location.model_id.value(), 0);
        assert_eq!(
            implicit_free.grouping.models[0].source_model_number,
            Some(1)
        );

        let mut implicit_next_free = PdbModelTransition::default();
        implicit_next_free
            .model_record(&pdb_model_source_line(1))
            .unwrap();
        implicit_next_free.end_model_record();
        let location = add_transition_pdb_atom(&mut implicit_next_free).unwrap();
        assert_eq!(location.model_id.value(), 1);
        assert_eq!(
            implicit_next_free.grouping.models[1].source_model_number,
            Some(2)
        );

        let mut implicit_collision = PdbModelTransition::default();
        implicit_collision
            .model_record(&pdb_model_source_line(2))
            .unwrap();
        implicit_collision.end_model_record();
        assert!(matches!(
            add_transition_pdb_atom(&mut implicit_collision),
            Err(PdbModelTransitionError::AtomBetweenModels)
        ));
    }

    #[test]
    fn gemmi_pdb_model_transition_matches_end_eof_and_wrong_format_gates() {
        // Gemmi's END breaks record dispatch; post-loop finalization adds one
        // empty model only if none was read. Format detection is disabled only
        // while its model pointer is active, including after MODEL before ATOM.
        let mut empty = PdbModelTransition::default();
        empty.finish_eof().unwrap();
        assert_eq!(empty.grouping.models.len(), 1);
        assert_eq!(empty.grouping.models[0].source_model_number, Some(1));

        let mut explicit = PdbModelTransition::default();
        explicit.model_record(&pdb_model_source_line(5)).unwrap();
        explicit.end_record();
        assert!(explicit.stopped_at_end_record());
        explicit.finish_eof().unwrap();
        assert_eq!(explicit.grouping.models.len(), 1);
        assert_eq!(explicit.grouping.models[0].source_model_number, Some(5));

        for (input, source, expected) in [
            (
                "DaTa_model\n",
                "folder/model.cif",
                PdbInputFormatError::Cif {
                    source: "folder/model.cif".to_owned(),
                },
            ),
            (
                "{\"DaTa_\n",
                "folder/model.json",
                PdbInputFormatError::Mmjson {
                    source: "folder/model.json".to_owned(),
                },
            ),
        ] {
            let mut no_model = PdbModelTransition::default();
            assert_eq!(
                transition_format_error(&no_model, input, source),
                Some(expected.clone())
            );

            let mut active_model = PdbModelTransition::default();
            active_model
                .model_record(&pdb_model_source_line(1))
                .unwrap();
            assert_eq!(transition_format_error(&active_model, input, source), None);
            active_model.end_model_record();
            assert_eq!(
                transition_format_error(&active_model, input, source),
                Some(expected)
            );
        }
    }

    #[test]
    fn gemmi_pdb_ter_type3_suffix_ignore_and_missing_chain_follow_source() {
        // Gemmi is_record_type3 accepts the source's low-nibble suffix set;
        // pinned source comments list space/LF/CR/TAB. The pinned P18 oracle
        // confirms LF is handled and TERE is not. ignore_ter and missing
        // current-chain paths leave the durable TER status untouched.
        for suffix in [b' ', b'\n', b'\r', b'\t', 0] {
            let mut transition = PdbModelTransition::default();
            let mut state = BioStructureSourceState::default();
            add_transition_pdb_atom_named(&mut transition, *b"ATOM  ", *b"ALA", b'A', 1).unwrap();
            let line = transition_ter_line(suffix);
            transition
                .ter_record(&line, PdbInputOptions::default(), &mut state)
                .unwrap();
            assert_eq!(state.ter_status, b'y', "suffix {suffix:?}");
            assert_eq!(
                transition.grouping.models[0].chains[0].residues[0].entity_kind,
                EntityKind::Polymer,
                "suffix {suffix:?}"
            );
        }

        for suffix in [b'E', b'1'] {
            let mut transition = PdbModelTransition::default();
            let mut state = BioStructureSourceState::default();
            add_transition_pdb_atom_named(&mut transition, *b"ATOM  ", *b"ALA", b'A', 1).unwrap();
            let line = transition_ter_line(suffix);
            transition
                .ter_record(&line, PdbInputOptions::default(), &mut state)
                .unwrap();
            assert_eq!(state.ter_status, 0, "suffix {suffix:?}");
            assert_eq!(
                transition.grouping.models[0].chains[0].residues[0].entity_kind,
                EntityKind::Unknown,
                "suffix {suffix:?}"
            );
        }

        let mut ignored = PdbModelTransition::default();
        let mut ignored_state = BioStructureSourceState::default();
        add_transition_pdb_atom_named(&mut ignored, *b"ATOM  ", *b"ALA", b'A', 1).unwrap();
        let line = transition_ter_line(b' ');
        ignored
            .ter_record(
                &line,
                PdbInputOptions {
                    ignore_ter: true,
                    ..PdbInputOptions::default()
                },
                &mut ignored_state,
            )
            .unwrap();
        assert_eq!(ignored_state.ter_status, 0);
        assert_eq!(
            ignored.grouping.models[0].chains[0].residues[0].entity_kind,
            EntityKind::Unknown
        );

        let mut no_chain = PdbModelTransition::default();
        let mut no_chain_state = BioStructureSourceState::default();
        no_chain
            .ter_record(&line, PdbInputOptions::default(), &mut no_chain_state)
            .unwrap();
        assert_eq!(no_chain_state.ter_status, 0);
    }

    #[test]
    fn gemmi_pdb_ter_marks_current_chain_and_classifies_later_same_chain_residues() {
        // Pinned Gemmi 5cc1c23 output in IO-bio_pdb.md, Step264:
        // normal_same_chain|ter=y|...|m1:[A:ALA1=1,GLY2=2,]. The first
        // TER promotes all then-current residues; a later new residue in the
        // same chain is NonPolymer because after_ter remains true.
        let mut transition = PdbModelTransition::default();
        let mut state = BioStructureSourceState::default();
        add_transition_pdb_atom_named(&mut transition, *b"ATOM  ", *b"ALA", b'A', 1).unwrap();
        let line = transition_ter_line(b' ');
        transition
            .ter_record(&line, PdbInputOptions::default(), &mut state)
            .unwrap();
        add_transition_pdb_atom_named(&mut transition, *b"ATOM  ", *b"GLY", b'A', 2).unwrap();

        assert_eq!(state.ter_status, b'y');
        let residues = &transition.grouping.models[0].chains[0].residues;
        assert_eq!(residues.len(), 2);
        assert_eq!(residues[0].entity_kind, EntityKind::Polymer);
        assert_eq!(residues[1].entity_kind, EntityKind::NonPolymer);
    }

    #[test]
    fn gemmi_pdb_ter_split_chain_preserves_unknown_entity_kinds() {
        // The pinned oracle prints two A chain parts and durable status y;
        // neither part is typed Polymer because split_chain_on_ter bypasses
        // entity inference.
        let mut transition = PdbModelTransition::default();
        let mut state = BioStructureSourceState::default();
        add_transition_pdb_atom_named(&mut transition, *b"ATOM  ", *b"ALA", b'A', 1).unwrap();
        let line = transition_ter_line(b' ');
        transition
            .ter_record(
                &line,
                PdbInputOptions {
                    split_chain_on_ter: true,
                    ..PdbInputOptions::default()
                },
                &mut state,
            )
            .unwrap();
        add_transition_pdb_atom_named(&mut transition, *b"ATOM  ", *b"GLY", b'A', 2).unwrap();

        assert_eq!(state.ter_status, b'y');
        let chains = &transition.grouping.models[0].chains;
        assert_eq!(chains.len(), 2);
        assert!(chains.iter().all(|chain| {
            chain.residues.len() == 1 && chain.residues[0].entity_kind == EntityKind::Unknown
        }));
    }

    #[test]
    fn gemmi_pdb_ter_reappearing_chain_uses_first_prior_part_and_water_kind() {
        // Pinned Gemmi output in IO-bio_pdb.md, Step264: reappearing A after
        // Polymer A / unseen B is NonPolymer, while DOD in that reappearing
        // part is Water. Both decisions use the first prior A part.
        let mut transition = PdbModelTransition::default();
        let mut state = BioStructureSourceState::default();
        add_transition_pdb_atom_named(&mut transition, *b"ATOM  ", *b"ALA", b'A', 1).unwrap();
        let line = transition_ter_line(b' ');
        transition
            .ter_record(&line, PdbInputOptions::default(), &mut state)
            .unwrap();
        add_transition_pdb_atom_named(&mut transition, *b"ATOM  ", *b"GLY", b'B', 1).unwrap();
        add_transition_pdb_atom_named(&mut transition, *b"HETATM", *b"DOD", b'A', 3).unwrap();

        assert_eq!(state.ter_status, b'y');
        let chains = &transition.grouping.models[0].chains;
        assert_eq!(chains.len(), 3);
        assert_eq!(chains[0].residues[0].entity_kind, EntityKind::Polymer);
        assert_eq!(chains[1].residues[0].entity_kind, EntityKind::Unknown);
        assert_eq!(chains[2].residues[0].entity_kind, EntityKind::Water);
    }

    #[test]
    fn gemmi_pdb_ter_water_aliases_and_repeated_ter_roll_back_all_models() {
        // Gemmi 5cc1c23 oracle: HOH/DOD/WAT/H2O/h2o each set ter_status=e
        // when present in the current chain at TER; OH is not water. A second
        // TER also sets e. remove_entity_types then clears every model.
        for water_name in [*b"HOH", *b"DOD", *b"WAT", *b"H2O", *b"h2o"] {
            let mut transition = PdbModelTransition::default();
            let mut state = BioStructureSourceState::default();
            add_transition_pdb_atom_named(&mut transition, *b"ATOM  ", *b"ALA", b'A', 1).unwrap();
            add_transition_pdb_atom_named(&mut transition, *b"HETATM", water_name, b'A', 2)
                .unwrap();
            let line = transition_ter_line(b' ');
            transition
                .ter_record(&line, PdbInputOptions::default(), &mut state)
                .unwrap();
            assert_eq!(state.ter_status, b'e', "water {water_name:?}");
            assert!(
                transition.grouping.models[0].chains[0]
                    .residues
                    .iter()
                    .all(|residue| residue.entity_kind == EntityKind::Polymer)
            );
            transition.finalize_ter_entity_types(&state);
            assert!(
                transition.grouping.models[0].chains[0]
                    .residues
                    .iter()
                    .all(|residue| residue.entity_kind == EntityKind::Unknown)
            );
        }

        let mut transition = PdbModelTransition::default();
        let mut state = BioStructureSourceState::default();
        transition.model_record(&pdb_model_source_line(1)).unwrap();
        add_transition_pdb_atom_named(&mut transition, *b"ATOM  ", *b"ALA", b'A', 1).unwrap();
        let line = transition_ter_line(b' ');
        transition
            .ter_record(&line, PdbInputOptions::default(), &mut state)
            .unwrap();
        transition
            .ter_record(&line, PdbInputOptions::default(), &mut state)
            .unwrap();
        assert_eq!(state.ter_status, b'e');
        assert_eq!(
            transition.grouping.models[0].chains[0].residues[0].entity_kind,
            EntityKind::Polymer
        );
        transition.end_model_record();
        transition.model_record(&pdb_model_source_line(2)).unwrap();
        add_transition_pdb_atom_named(&mut transition, *b"ATOM  ", *b"GLY", b'B', 2).unwrap();
        transition.finalize_ter_entity_types(&state);
        assert!(transition.grouping.models.iter().all(|model| {
            model
                .chains
                .iter()
                .flat_map(|chain| &chain.residues)
                .all(|residue| residue.entity_kind == EntityKind::Unknown)
        }));
    }

    #[test]
    fn gemmi_pdb_ter_finalization_assigns_and_appends_subchains_in_source_order() {
        // Fixed Gemmi 5cc1c23 read_pdb_string observation for SEQRES A:ALA,
        // SEQRES B:GLY, ATOM A/ALA/1, TER, HETATM A/LIG/2, ATOM B/GLY/1,
        // TER, HETATM A/LIG/3, END:
        // model=1; A=[ALA,Polymer,Axp][LIG,NonPolymer,Ax1];
        // B=[GLY,Polymer,Bxp]; A=[LIG,NonPolymer,Ax2];
        // entity A=[Axp]; entity B=[Bxp].
        let mut transition = PdbModelTransition::default();
        let mut source_state = BioStructureSourceState::default();
        add_transition_pdb_atom_named(&mut transition, *b"ATOM  ", *b"ALA", b'A', 1).unwrap();
        let ter = transition_ter_line(b' ');
        transition
            .ter_record(&ter, PdbInputOptions::default(), &mut source_state)
            .unwrap();
        add_transition_pdb_atom_named(&mut transition, *b"HETATM", *b"LIG", b'A', 2).unwrap();
        add_transition_pdb_atom_named(&mut transition, *b"ATOM  ", *b"GLY", b'B', 1).unwrap();
        transition
            .ter_record(&ter, PdbInputOptions::default(), &mut source_state)
            .unwrap();
        add_transition_pdb_atom_named(&mut transition, *b"HETATM", *b"LIG", b'A', 3).unwrap();

        assert_eq!(source_state.ter_status, b'y');
        let mut entities = PdbEntityState {
            entities: [b'A', b'B']
                .into_iter()
                .map(|chain| PdbEntity {
                    source_name: vec![chain],
                    entity_kind: EntityKind::Polymer,
                    full_sequence: Vec::new(),
                    dbrefs: Vec::new(),
                    subchains: Vec::new(),
                })
                .collect(),
        };

        transition
            .finalize_subchains_and_entities(&source_state, &mut entities)
            .unwrap();

        let chains = &transition.grouping.models[0].chains;
        assert_eq!(chains.len(), 3);
        assert_eq!(
            chains
                .iter()
                .map(|chain| {
                    chain
                        .residues
                        .iter()
                        .map(|residue| residue.subchain_id.as_str())
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>(),
            vec![vec!["Axp", "Ax1"], vec!["Bxp"], vec!["Ax2"]]
        );
        assert_eq!(entities.entities[0].subchains, vec!["Axp".to_owned()]);
        assert_eq!(entities.entities[1].subchains, vec!["Bxp".to_owned()]);
    }

    #[test]
    fn gemmi_pdb_subchain_kind_suffixes_and_model_local_counter_match_source() {
        // Fixed Gemmi 5cc1c23 assign_subchains(false,false) observations from
        // a first model with 46 A NonPolymer rows, B rows of Polymer/Water/
        // Branched, then another A NonPolymer row, plus one A NonPolymer row
        // in model 2: A[1,9,10,11,45,46]=Ax1,Ax9,Ax0,Ax01,Ax0Z,Ax10;
        // B= Bxp,Bxw,Bxb; repeated A=Ax11; model-2 A=Ax1.
        let mut transition = PdbModelTransition::default();
        for sequence in 1..=46 {
            add_transition_pdb_atom_named(&mut transition, *b"HETATM", *b"LIG", b'A', sequence)
                .unwrap();
        }
        add_transition_pdb_atom_named(&mut transition, *b"ATOM  ", *b"ALA", b'B', 1).unwrap();
        add_transition_pdb_atom_named(&mut transition, *b"HETATM", *b"LIG", b'B', 2).unwrap();
        add_transition_pdb_atom_named(&mut transition, *b"HETATM", *b"LIG", b'B', 3).unwrap();
        add_transition_pdb_atom_named(&mut transition, *b"HETATM", *b"LIG", b'A', 100).unwrap();

        let second_model = transition.grouping.add_model(Some(2)).unwrap();
        let second_model_line =
            grouped_pdb_atom_line(*b"HETATM", 101, *b" C1 ", *b"LIG", b'A', 1, b' ', *b"    ");
        add_grouped_pdb_atom(&mut transition.grouping, second_model, &second_model_line);

        for residue in &mut transition.grouping.models[0].chains[0].residues {
            residue.entity_kind = EntityKind::NonPolymer;
        }
        for (residue, kind) in transition.grouping.models[0].chains[1]
            .residues
            .iter_mut()
            .zip([EntityKind::Polymer, EntityKind::Water, EntityKind::Branched])
        {
            residue.entity_kind = kind;
        }
        transition.grouping.models[0].chains[2].residues[0].entity_kind = EntityKind::NonPolymer;
        transition.grouping.models[1].chains[0].residues[0].entity_kind = EntityKind::NonPolymer;

        assign_pdb_subchains(&mut transition.grouping.models).unwrap();

        let first_a = &transition.grouping.models[0].chains[0].residues;
        for (index, expected) in ["Ax1", "Ax9", "Ax0", "Ax01", "Ax0Z", "Ax10"]
            .into_iter()
            .enumerate()
        {
            let source_index = [0, 8, 9, 10, 44, 45][index];
            assert_eq!(first_a[source_index].subchain_id, expected);
        }
        assert_eq!(
            transition.grouping.models[0].chains[1]
                .residues
                .iter()
                .map(|residue| residue.subchain_id.as_str())
                .collect::<Vec<_>>(),
            ["Bxp", "Bxw", "Bxb"]
        );
        assert_eq!(
            transition.grouping.models[0].chains[2].residues[0].subchain_id,
            "Ax11"
        );
        assert_eq!(
            transition.grouping.models[1].chains[0].residues[0].subchain_id,
            "Ax1"
        );
    }

    #[test]
    fn gemmi_pdb_subchain_assignment_preserves_source_existing_and_unknown_branches() {
        // Fixed Gemmi 5cc1c23 assign_subchains(false,false) output:
        // fully-subchained C=[Unknown/old1,Polymer/old2] stays intact;
        // partial fully-typed D becomes Dxp,Dx1; partial D-like E containing
        // Unknown remains [Unknown/keep,Polymer/empty]; empty F is unchanged.
        let mut transition = PdbModelTransition::default();
        for (chain, sequence) in [
            (b'C', 1),
            (b'C', 2),
            (b'D', 1),
            (b'D', 2),
            (b'E', 1),
            (b'E', 2),
        ] {
            add_transition_pdb_atom_named(&mut transition, *b"ATOM  ", *b"ALA", chain, sequence)
                .unwrap();
        }
        let chains = &mut transition.grouping.models[0].chains;
        chains[0].residues[0].entity_kind = EntityKind::Unknown;
        chains[0].residues[0].subchain_id = "old1".to_owned();
        chains[0].residues[1].entity_kind = EntityKind::Polymer;
        chains[0].residues[1].subchain_id = "old2".to_owned();
        chains[1].residues[0].entity_kind = EntityKind::Polymer;
        chains[1].residues[0].subchain_id = "stale".to_owned();
        chains[1].residues[1].entity_kind = EntityKind::NonPolymer;
        chains[2].residues[0].entity_kind = EntityKind::Unknown;
        chains[2].residues[0].subchain_id = "keep".to_owned();
        chains[2].residues[1].entity_kind = EntityKind::Polymer;
        chains.push(super::PdbChainGroup {
            row_id: cosmolkit_bio::BioChainId::new(3),
            source_id: PdbChainId::from_ascii(b"F").unwrap(),
            residues: Vec::new(),
            residue_indices: std::collections::HashMap::new(),
            current_residue: None,
        });

        assign_pdb_subchains(&mut transition.grouping.models).unwrap();

        assert_eq!(
            transition.grouping.models[0].chains[0]
                .residues
                .iter()
                .map(|residue| (residue.entity_kind, residue.subchain_id.as_str()))
                .collect::<Vec<_>>(),
            [(EntityKind::Unknown, "old1"), (EntityKind::Polymer, "old2")]
        );
        assert_eq!(
            transition.grouping.models[0].chains[1]
                .residues
                .iter()
                .map(|residue| residue.subchain_id.as_str())
                .collect::<Vec<_>>(),
            ["Dxp", "Dx1"]
        );
        assert_eq!(
            transition.grouping.models[0].chains[2]
                .residues
                .iter()
                .map(|residue| (residue.entity_kind, residue.subchain_id.as_str()))
                .collect::<Vec<_>>(),
            [(EntityKind::Unknown, "keep"), (EntityKind::Polymer, "")]
        );
        assert!(transition.grouping.models[0].chains[3].residues.is_empty());
    }

    #[test]
    fn gemmi_pdb_subchain_assignment_keeps_chain_width_and_counter_error_boundaries() {
        assert!(PdbChainId::from_ascii(b"ABCD").is_some());
        assert!(PdbChainId::from_ascii(b"ABCDE").is_none());

        let mut transition = PdbModelTransition::default();
        add_transition_pdb_atom_named(&mut transition, *b"ATOM  ", *b"ALA", b'A', 1).unwrap();
        let chain = &mut transition.grouping.models[0].chains[0];
        chain.source_id = PdbChainId::from_ascii(b"ABCD").unwrap();
        chain.residues[0].entity_kind = EntityKind::Polymer;
        assign_pdb_subchains(&mut transition.grouping.models).unwrap();
        assert_eq!(
            transition.grouping.models[0].chains[0].residues[0].subchain_id,
            "ABCDxp"
        );

        let chain = &mut transition.grouping.models[0].chains[0];
        chain.residues[0].entity_kind = EntityKind::NonPolymer;
        chain.residues[0].subchain_id.clear();
        let chain_id = chain.source_id;
        let mut counter = i32::MAX;
        let error = assign_pdb_subchain_names(chain, &mut counter).unwrap_err();
        assert!(matches!(
            error,
            PdbModelTransitionError::NonpolymerSubchainCounterOutsideSourceIntRange { chain }
                if chain == chain_id
        ));
        assert!(
            transition.grouping.models[0].chains[0].residues[0]
                .subchain_id
                .is_empty()
        );
    }

    #[test]
    fn gemmi_pdb_anisou_associates_to_last_atom_and_rounds_six_components() {
        let mut transition = PdbModelTransition::default();
        add_transition_pdb_atom(&mut transition).unwrap();
        add_transition_pdb_atom(&mut transition).unwrap();

        // The pinned Gemmi input uses an ANISOU serial/name that do not match
        // either atom. `populate_structure_from_pdb_stream` associates it with
        // `resi->atoms.back()` without comparing those fields.
        let line =
            pdb_anisou_source_line(999, *b" XX ", [12345, -23456, 34567, -45678, 56789, -67890]);
        transition.anisou_record(&line).unwrap();

        let atoms = &transition.grouping.models[0].chains[0].residues[0].atoms;
        assert_eq!(atoms.len(), 2);
        assert_eq!(atoms[0].anisou, [0.0; 6]);
        let expected_source_f32_bits = [
            0x3f9e_0418_u32,
            0xc016_1e4f,
            0x405d_3a92,
            0xc092_2b6b,
            0x40b5_b98c,
            0xc0d9_3f7d,
        ];
        for (observed, source_bits) in atoms[1].anisou.iter().zip(expected_source_f32_bits) {
            assert_eq!(
                observed.to_bits(),
                f64::from(f32::from_bits(source_bits)).to_bits()
            );
        }
    }

    #[test]
    fn gemmi_pdb_anisou_requires_current_model_chain_residue_and_atom() {
        let line = pdb_anisou_source_line(1, *b" CA ", [1, 2, 3, 4, 5, 6]);

        let mut no_model = PdbModelTransition::default();
        let error = no_model.anisou_record(&line).unwrap_err();
        assert!(matches!(error, PdbModelTransitionError::AnisouWithoutAtom));
        assert_eq!(
            error.to_string(),
            "ANISOU record not directly after ATOM/HETATM."
        );

        let mut no_chain = PdbModelTransition::default();
        no_chain.model_record(&pdb_model_source_line(1)).unwrap();
        let error = no_chain.anisou_record(&line).unwrap_err();
        assert!(matches!(error, PdbModelTransitionError::AnisouWithoutAtom));
        assert_eq!(
            error.to_string(),
            "ANISOU record not directly after ATOM/HETATM."
        );

        let mut no_residue = PdbModelTransition::default();
        add_transition_pdb_atom(&mut no_residue).unwrap();
        no_residue.grouping.models[0].chains[0].current_residue = None;
        let error = no_residue.anisou_record(&line).unwrap_err();
        assert!(matches!(error, PdbModelTransitionError::AnisouWithoutAtom));
        assert_eq!(
            error.to_string(),
            "ANISOU record not directly after ATOM/HETATM."
        );

        let mut empty_residue = PdbModelTransition::default();
        add_transition_pdb_atom(&mut empty_residue).unwrap();
        empty_residue.grouping.models[0].chains[0].residues[0]
            .atoms
            .clear();
        let error = empty_residue.anisou_record(&line).unwrap_err();
        assert!(matches!(error, PdbModelTransitionError::AnisouWithoutAtom));
        assert_eq!(
            error.to_string(),
            "ANISOU record not directly after ATOM/HETATM."
        );
    }

    #[test]
    fn gemmi_pdb_anisou_uses_only_nonzero_u11_as_duplicate_sentinel() {
        let mut nonzero_u11 = PdbModelTransition::default();
        add_transition_pdb_atom(&mut nonzero_u11).unwrap();
        nonzero_u11
            .anisou_record(&pdb_anisou_source_line(
                1,
                *b" CA ",
                [10000, 20000, 30000, 40000, 50000, 60000],
            ))
            .unwrap();
        let error = nonzero_u11
            .anisou_record(&pdb_anisou_source_line(
                1,
                *b" CA ",
                [12345, -23456, 34567, -45678, 56789, -67890],
            ))
            .unwrap_err();
        assert!(matches!(error, PdbModelTransitionError::DuplicateAnisou));
        assert_eq!(
            error.to_string(),
            "Duplicated ANISOU record or not directly after ATOM/HETATM."
        );

        let mut zero_u11 = PdbModelTransition::default();
        add_transition_pdb_atom(&mut zero_u11).unwrap();
        zero_u11
            .anisou_record(&pdb_anisou_source_line(
                1,
                *b" CA ",
                [0, 111, 222, 333, 444, 555],
            ))
            .unwrap();
        zero_u11
            .anisou_record(&pdb_anisou_source_line(
                1,
                *b" CA ",
                [10000, 20000, 30000, 40000, 50000, 60000],
            ))
            .unwrap();
        let expected_source_f32_bits = [
            0x3f80_0000_u32,
            0x4000_0000,
            0x4040_0000,
            0x4080_0000,
            0x40a0_0000,
            0x40c0_0000,
        ];
        let observed = &zero_u11.grouping.models[0].chains[0].residues[0].atoms[0].anisou;
        for (observed, source_bits) in observed.iter().zip(expected_source_f32_bits) {
            assert_eq!(
                observed.to_bits(),
                f64::from(f32::from_bits(source_bits)).to_bits()
            );
        }
    }

    #[test]
    fn pdb_line_cursor_keeps_newline_and_returns_final_unterminated_line_once() {
        let mut cursor = PdbLineCursor::new("one\nlast", PdbInputOptions::default());
        assert_eq!(cursor.copy_line(), Some(&b"one\n"[..]));
        assert_eq!(cursor.copy_line(), Some(&b"last"[..]));
        assert_eq!(cursor.copy_line(), None);

        let mut one_byte = PdbLineCursor::new(
            "A\nB",
            PdbInputOptions {
                max_line_length: 1,
                ..PdbInputOptions::default()
            },
        );
        assert_eq!(one_byte.copy_line(), Some(&b"A"[..]));
        assert_eq!(one_byte.copy_line(), Some(&b"B"[..]));
        assert_eq!(one_byte.copy_line(), None);
    }

    #[test]
    fn pdb_line_cursor_discards_truncated_suffix_before_reading_next_line() {
        let mut cursor = PdbLineCursor::new(
            "abcdefgh\nSECOND\n",
            PdbInputOptions {
                max_line_length: 5,
                ..PdbInputOptions::default()
            },
        );
        assert_eq!(cursor.copy_line(), Some(&b"abcde"[..]));
        assert_eq!(cursor.copy_line(), Some(&b"SECON"[..]));
        assert_eq!(cursor.copy_line(), None);

        let exact_limit = format!("{}\nNEXT\n", "x".repeat(120));
        let mut cursor = PdbLineCursor::new(
            &exact_limit,
            PdbInputOptions {
                max_line_length: 120,
                ..PdbInputOptions::default()
            },
        );
        assert_eq!(cursor.copy_line(), Some(&b"x".repeat(120)[..]));
        assert_eq!(cursor.copy_line(), Some(&b"NEXT\n"[..]));
        assert_eq!(cursor.copy_line(), None);
    }

    #[test]
    fn pdb_line_cursor_matches_gemmi_strlen_and_signed_char_drain_edges() {
        let input = "\0first\nnext\n";
        let mut nul_prefix = PdbLineCursor::new(input, PdbInputOptions::default());
        assert_eq!(nul_prefix.copy_line(), None);
        assert_eq!(nul_prefix.position, 7);

        let mut signed_char = PdbLineCursor::new(
            "ABérest\nNEXT\n",
            PdbInputOptions {
                max_line_length: 2,
                ..PdbInputOptions::default()
            },
        );
        assert_eq!(signed_char.copy_line(), Some(&b"AB"[..]));
        assert_eq!(signed_char.copy_line(), Some(&[0xa9, b'r'][..]));
        assert_eq!(signed_char.copy_line(), Some(&b"NE"[..]));
        assert_eq!(signed_char.copy_line(), None);
    }

    #[test]
    fn gemmi_pdb_read_string_trims_c_locale_whitespace_in_source_order() {
        assert_eq!(read_string(b"\t\n\x0b\x0c\r CA \t\x0b\x0c"), b"CA");
        assert_eq!(read_string(b"\nCA"), b"CA");
        assert_eq!(read_string(b"CA \t\x0b\x0c"), b"CA");
        assert_eq!(read_string(b"   \t\x0b\x0c"), b"");
    }

    #[test]
    fn gemmi_pdb_read_string_terminates_at_cr_lf_and_nul() {
        assert_eq!(read_string(b"CA\rignored"), b"CA");
        assert_eq!(read_string(b"ABC\nignored"), b"ABC");
        assert_eq!(read_string(b"CA\0ignored"), b"CA");
        assert_eq!(read_string(b"  \rignored"), b"ignored");
    }

    #[test]
    fn raw_pdb_atom_name_preserves_four_columns_separately_from_read_string() {
        let leading_and_trailing_spaces = raw_atom_name(b" CA ").unwrap();
        assert_eq!(leading_and_trailing_spaces.as_bytes(), b" CA ");
        assert_eq!(read_string(leading_and_trailing_spaces.as_bytes()), b"CA");

        let trailing_spaces = raw_atom_name(b"CA  ").unwrap();
        assert_eq!(trailing_spaces.as_bytes(), b"CA  ");
        assert!(raw_atom_name(b" CA").is_none());
        assert!(raw_atom_name(b"CA   ").is_none());
        assert!(raw_atom_name(&[0xc3, b'C', b'A', b' ']).is_none());
    }

    #[test]
    fn pdb_atom_name_width_preserves_raw_columns_for_source_element_inference() {
        let mut line = pdb_atom_field_line();
        for (raw, element, isotope) in [
            (*b" CA ", Element::C, None),
            (*b"CA  ", Element::CA, None),
            (*b" D  ", Element::H, Some(2)),
            (*b"DXYZ", Element::H, Some(2)),
        ] {
            line[12..16].copy_from_slice(&raw);
            let fields = decode_pdb_atom_field_line(&line[..76]).unwrap();
            assert_eq!(fields.name.as_bytes(), raw);
            assert_eq!(
                (fields.element, fields.isotope_mass_number),
                (element, isotope)
            );
            assert_eq!(infer_element_from_padded_name(&raw), (element, isotope));
        }
    }

    #[test]
    fn pdb_atom_name_width_link_keeps_raw_fields_for_classification_and_logical_address() {
        let record = p39_link_record(
            *b"LINK  ",
            [*b" CA ", *b"CA  "],
            [b' ', b' '],
            Some([*b"1555  ", *b"1555  "]),
            None,
            None,
            80,
        );
        let (first_address, first_raw) = decode_pdb_link_address(&record, 1, 0).unwrap();
        let (second_address, second_raw) = decode_pdb_link_address(&record, 2, 30).unwrap();
        assert_eq!(first_raw, *b" CA ");
        assert_eq!(second_raw, *b"CA  ");
        assert_eq!(first_address.logical_atom_name(), "CA");
        assert_eq!(second_address.logical_atom_name(), "CA");
        let connection = decode_pdb_link_record(
            &record,
            &p38_model_group(&[]),
            &mut PdbLinkNameCounters::default(),
        )
        .unwrap()
        .unwrap();
        assert_eq!(connection.kind, BioConnectionKind::MetalC);
        assert_eq!(connection.partner1.logical_atom_name(), "CA");
        assert_eq!(connection.partner2.logical_atom_name(), "CA");
    }

    fn pdb_atom_field_line() -> [u8; 80] {
        let mut line = [b' '; 80];
        line[..6].copy_from_slice(b"ATOM  ");
        line[6..11].copy_from_slice(b"    1");
        line[12..16].copy_from_slice(b" CA ");
        line[16] = b'B';
        line[17..20].copy_from_slice(b"GLY");
        line[21] = b'A';
        line[22..26].copy_from_slice(b"   1");
        line[26] = b'A';
        line[30..38].copy_from_slice(b"   1.250");
        line[38..46].copy_from_slice(b"  -2.500");
        line[46..54].copy_from_slice(b"   3.750");
        line[54..60].copy_from_slice(b"  0.50");
        line[60..66].copy_from_slice(b" 12.25");
        line[76..78].copy_from_slice(b" C");
        line[78..80].copy_from_slice(b"2+");
        line
    }

    fn decode_pdb_atom_field_line(line: &[u8]) -> Result<super::PdbAtomFields, PdbAtomFieldError> {
        let mut source_line_buffer = [0; 122];
        source_line_buffer[..line.len()].copy_from_slice(line);
        decode_pdb_atom_fields(line, &source_line_buffer)
    }

    fn grouped_pdb_atom_line(
        record: [u8; 6],
        serial: i32,
        atom_name: [u8; 4],
        residue_name: [u8; 3],
        chain_name: u8,
        sequence_number: i32,
        insertion_code: u8,
        segment: [u8; 4],
    ) -> [u8; 80] {
        let mut line = pdb_atom_field_line();
        line[..6].copy_from_slice(&record);
        let serial_field = format!("{serial:>5}");
        assert_eq!(serial_field.len(), 5);
        line[6..11].copy_from_slice(serial_field.as_bytes());
        line[12..16].copy_from_slice(&atom_name);
        line[17..20].copy_from_slice(&residue_name);
        line[21] = chain_name;
        let sequence_field = format!("{sequence_number:>4}");
        assert_eq!(sequence_field.len(), 4);
        line[22..26].copy_from_slice(sequence_field.as_bytes());
        line[26] = insertion_code;
        line[72..76].copy_from_slice(&segment);
        line
    }

    fn pdb_model_source_line(number: i32) -> [u8; 122] {
        let mut line = [b' '; 122];
        line[..5].copy_from_slice(b"MODEL");
        let field = format!("{number:>8}");
        assert_eq!(field.len(), 8);
        line[6..14].copy_from_slice(field.as_bytes());
        line[14] = b'\n';
        line
    }

    fn pdb_anisou_source_line(serial: i32, atom_name: [u8; 4], values: [i32; 6]) -> [u8; 122] {
        let mut line = [b' '; 122];
        line[..6].copy_from_slice(b"ANISOU");
        let serial_field = format!("{serial:>5}");
        assert_eq!(serial_field.len(), 5);
        line[6..11].copy_from_slice(serial_field.as_bytes());
        line[12..16].copy_from_slice(&atom_name);
        line[17..20].copy_from_slice(b"GLY");
        line[21] = b'A';
        line[22..26].copy_from_slice(b"   1");
        for (value, offset) in values.into_iter().zip([28_usize, 35, 42, 49, 56, 63]) {
            let field = format!("{value:>7}");
            assert_eq!(field.len(), 7);
            line[offset..offset + 7].copy_from_slice(field.as_bytes());
        }
        line[80] = b'\n';
        line
    }

    fn add_transition_pdb_atom(
        transition: &mut PdbModelTransition,
    ) -> Result<PdbGroupLocation, PdbModelTransitionError> {
        add_transition_pdb_atom_named(transition, *b"ATOM  ", *b"GLY", b'A', 1)
    }

    fn add_transition_pdb_atom_named(
        transition: &mut PdbModelTransition,
        record: [u8; 6],
        residue_name: [u8; 3],
        chain_name: u8,
        sequence_number: i32,
    ) -> Result<PdbGroupLocation, PdbModelTransitionError> {
        let line = grouped_pdb_atom_line(
            record,
            sequence_number,
            *b" CA ",
            residue_name,
            chain_name,
            sequence_number,
            b' ',
            *b"    ",
        );
        let fields = decode_pdb_atom_field_line(&line).unwrap();
        let chain_name = PdbChainId::from_ascii(read_string(&line[20..22])).unwrap();
        let sequence_field: &[u8; 5] = line[22..27].try_into().unwrap();
        let residue_name_field: &[u8; 3] = line[17..20].try_into().unwrap();
        let residue = read_res_id(sequence_field, residue_name_field).unwrap();
        let residue = read_pdb_residue_segment(residue, &line).unwrap();
        transition.add_atom_record(chain_name, residue, line[0], fields)
    }

    fn transition_ter_line(suffix: u8) -> [u8; 4] {
        [b'T', b'E', b'R', suffix]
    }

    fn seqres_source_line(
        record: [u8; 6],
        chain: [u8; 2],
        slots: &[(usize, [u8; 3])],
        length: usize,
    ) -> Vec<u8> {
        let mut line = vec![b' '; 68];
        line[..6].copy_from_slice(&record);
        line[10..12].copy_from_slice(&chain);
        for (offset, residue_name) in slots {
            line[*offset..*offset + 3].copy_from_slice(residue_name);
        }
        line.truncate(length);
        line
    }

    fn transition_format_error(
        transition: &PdbModelTransition,
        input: &str,
        source: &str,
    ) -> Option<PdbInputFormatError> {
        let mut cursor = PdbLineCursor::new(input, PdbInputOptions::default());
        let line = cursor.copy_line()?;
        let mut reader = PdbReaderState::new(source, false, BioStructureSourceState::default());
        reader.record_line(line);
        transition.wrong_input_format(&reader, cursor.source_line_buffer(), source)
    }

    fn add_grouped_pdb_atom(
        grouping: &mut PdbHierarchyGrouping,
        model_id: cosmolkit_bio::BioModelId,
        line: &[u8; 80],
    ) -> PdbGroupLocation {
        let fields = decode_pdb_atom_field_line(line).unwrap();
        let chain_name = PdbChainId::from_ascii(read_string(&line[20..22])).unwrap();
        let sequence_field: &[u8; 5] = line[22..27].try_into().unwrap();
        let residue_name_field: &[u8; 3] = line[17..20].try_into().unwrap();
        let residue = read_res_id(sequence_field, residue_name_field).unwrap();
        let residue = read_pdb_residue_segment(residue, line).unwrap();
        grouping
            .add_atom_record(model_id, chain_name, residue, line[0], fields)
            .unwrap()
    }

    type P38AtomRow = (i32, [u8; 4], u8, u8, i32, [u8; 4], f64, [u8; 2]);

    fn p38_model_group(rows: &[P38AtomRow]) -> PdbModelGroup {
        let mut grouping = PdbHierarchyGrouping::default();
        let model_id = grouping.add_model(Some(1)).unwrap();
        for (serial, atom_name, altloc, chain, sequence, segment, x, element) in rows {
            let mut line = grouped_pdb_atom_line(
                *b"ATOM  ", *serial, *atom_name, *b"CYS", *chain, *sequence, b' ', *segment,
            );
            line[16] = *altloc;
            for (range, value) in [(30..38, *x), (38..46, 0.0), (46..54, 0.0)] {
                let field = format!("{value:8.3}");
                assert_eq!(field.len(), range.len());
                line[range].copy_from_slice(field.as_bytes());
            }
            line[76..78].copy_from_slice(element);
            add_grouped_pdb_atom(&mut grouping, model_id, &line);
        }
        grouping
            .models
            .pop()
            .expect("the P38 fixture has one model")
    }

    fn p38_residue(sequence: i32, segment: &[u8]) -> ResidueAddress {
        ResidueAddress::new(
            Some(sequence),
            None,
            segment,
            ResidueName::from_ascii(b"CYS").unwrap(),
        )
        .unwrap()
    }

    fn p38_address(chain: u8, sequence: i32) -> AtomAddress {
        AtomAddress::new(
            PdbChainId::from_ascii(&[chain]).unwrap(),
            p38_residue(sequence, b""),
            "",
            None,
        )
    }

    fn p38_connection(asu: BioAsu) -> BioConnection {
        BioConnection {
            name: "disulf1".to_owned(),
            kind: BioConnectionKind::Disulf,
            asu,
            partner1: p38_address(b'A', 1),
            partner2: p38_address(b'B', 2),
            ..BioConnection::default()
        }
    }

    fn p38_prepared_crystal() -> BioCrystalInfo {
        let mut crystal = BioCrystalInfo::new(
            BioCrystalCell::default(),
            None,
            None,
            BioTransform::identity(),
            BioTransform::identity(),
            false,
            0,
            Vec::new(),
        );
        cosmolkit_bio::set_crystal_cell(
            &mut crystal,
            BioCrystalCell {
                a: 10.0,
                b: 10.0,
                c: 10.0,
                alpha: 90.0,
                beta: 90.0,
                gamma: 90.0,
            },
        )
        .unwrap();
        cosmolkit_bio::set_crystal_space_group_hm(&mut crystal, "P 1".to_owned());
        cosmolkit_bio::setup_cell_images(&mut crystal, &[]);
        crystal
    }

    fn p39_link_record(
        record_name: [u8; 6],
        atom_names: [[u8; 4]; 2],
        altlocs: [u8; 2],
        symmetries: Option<[[u8; 6]; 2]>,
        distance: Option<[u8; 5]>,
        link_id: Option<[u8; 8]>,
        length: usize,
    ) -> Vec<u8> {
        let mut record = vec![b' '; 80];
        record[..6].copy_from_slice(&record_name);
        for (index, offset) in [0, 30].into_iter().enumerate() {
            record[offset + 12..offset + 16].copy_from_slice(&atom_names[index]);
            record[offset + 16] = altlocs[index];
            record[offset + 17..offset + 20].copy_from_slice(b"LIG");
            record[offset + 21] = if index == 0 { b'A' } else { b'B' };
            let sequence = if index == 0 { 1 } else { 12 };
            let sequence = format!("{sequence:>4}");
            assert_eq!(sequence.len(), 4);
            record[offset + 22..offset + 26].copy_from_slice(sequence.as_bytes());
            record[offset + 26] = b' ';
        }
        if let Some([first, second]) = symmetries {
            record[59..65].copy_from_slice(&first);
            record[66..72].copy_from_slice(&second);
        }
        if let Some(distance) = distance {
            record[73..78].copy_from_slice(&distance);
        }
        if let Some(link_id) = link_id {
            record[72..80].copy_from_slice(&link_id);
        }
        record.truncate(length);
        record
    }

    fn p40_cispep_record(
        record_name: [u8; 6],
        partner_c_chain: [u8; 2],
        partner_c_sequence: [u8; 5],
        partner_c_name: [u8; 3],
        partner_n_chain: [u8; 2],
        partner_n_sequence: [u8; 5],
        partner_n_name: [u8; 3],
        model_number: [u8; 3],
        angle: [u8; 6],
    ) -> Vec<u8> {
        let mut record = vec![b' '; 80];
        record[..6].copy_from_slice(&record_name);
        record[11..14].copy_from_slice(&partner_c_name);
        record[14..16].copy_from_slice(&partner_c_chain);
        record[17..22].copy_from_slice(&partner_c_sequence);
        record[25..28].copy_from_slice(&partner_n_name);
        record[28..30].copy_from_slice(&partner_n_chain);
        record[31..36].copy_from_slice(&partner_n_sequence);
        record[43..46].copy_from_slice(&model_number);
        record[53..59].copy_from_slice(&angle);
        record
    }

    fn p40_models(numbers: &[i32]) -> Vec<PdbModelGroup> {
        let mut grouping = PdbHierarchyGrouping::default();
        for number in numbers {
            grouping.add_model(Some(*number)).unwrap();
        }
        grouping.models
    }

    fn p40_address(
        chain: &[u8],
        residue_name: &[u8],
        sequence: i32,
        insertion: Option<u8>,
    ) -> AtomAddress {
        AtomAddress::new(
            PdbChainId::from_ascii(chain).unwrap(),
            ResidueAddress::new(
                Some(sequence),
                insertion,
                b"",
                ResidueName::from_ascii(residue_name).unwrap(),
            )
            .unwrap(),
            "",
            None,
        )
    }

    fn p39_model_group_with_explicit_carbon_ca11() -> PdbModelGroup {
        let mut grouping = PdbHierarchyGrouping::default();
        let model_id = grouping.add_model(Some(1)).unwrap();
        let mut line =
            grouped_pdb_atom_line(*b"ATOM  ", 1, *b"CA11", *b"LIG", b'A', 1, b' ', *b"    ");
        // Gemmi's fixed reference row uses a blank altloc. The common ATOM
        // fixture defaults to B, which must not match a blank LINK request.
        line[16] = b' ';
        line[76..78].copy_from_slice(b" C");
        add_grouped_pdb_atom(&mut grouping, model_id, &line);
        grouping
            .models
            .pop()
            .expect("the P39 fixture has one model")
    }

    #[test]
    fn pdb_link_decodes_both_addresses_asu_symmetry_and_reported_distance() {
        // Pinned Gemmi 0.7.5 read_pdb_string probe: LINK C1(A/1)->O1(B/12),
        // altloc B, symmetry 1555/2556 and distance 1.23 produces covale1,
        // Different, reported_sym [2, 0, 0, 1], and the two trimmed addresses.
        // The exact source probe output is recorded in the P39 Step390 report.
        let record = p39_link_record(
            *b"LINK  ",
            [*b" C1 ", *b" O1 "],
            [b' ', b'B'],
            Some([*b"1555  ", *b"2556  "]),
            Some(*b" 1.23"),
            None,
            80,
        );
        let model = p38_model_group(&[]);
        let connection =
            decode_pdb_link_record(&record, &model, &mut PdbLinkNameCounters::default())
                .unwrap()
                .unwrap();

        assert_eq!(connection.name, "covale1");
        assert_eq!(connection.kind, BioConnectionKind::Covale);
        assert_eq!(connection.asu, BioAsu::Different);
        assert_eq!(connection.reported_sym, [2, 0, 0, 1]);
        assert_eq!(connection.reported_distance, 1.23);
        assert!(connection.link_id.is_empty());

        assert_eq!(
            connection.partner1.chain_name(),
            PdbChainId::from_ascii(b"A").unwrap()
        );
        assert_eq!(connection.partner1.logical_atom_name(), "C1");
        assert_eq!(connection.partner1.altloc(), 0);
        assert_eq!(
            connection.partner1.residue(),
            ResidueAddress::new(Some(1), None, b"", ResidueName::from_ascii(b"LIG").unwrap())
                .unwrap()
        );
        assert_eq!(
            connection.partner2.chain_name(),
            PdbChainId::from_ascii(b"B").unwrap()
        );
        assert_eq!(connection.partner2.logical_atom_name(), "O1");
        assert_eq!(connection.partner2.altloc(), b'B');
        assert_eq!(
            connection.partner2.residue(),
            ResidueAddress::new(
                Some(12),
                None,
                b"",
                ResidueName::from_ascii(b"LIG").unwrap()
            )
            .unwrap()
        );
    }

    #[test]
    fn pdb_link_symmetry_comparison_preserves_any_same_and_unparsed_different_cases() {
        let model = p38_model_group(&[]);
        let mut counters = PdbLinkNameCounters::default();

        let short = p39_link_record(
            *b"LINK  ",
            [*b" C1 ", *b" O1 "],
            [b' ', b' '],
            None,
            None,
            None,
            57,
        );
        let short = decode_pdb_link_record(&short, &model, &mut counters)
            .unwrap()
            .unwrap();
        assert_eq!(short.asu, BioAsu::Any);
        assert_eq!(short.reported_sym, [0; 4]);

        let same = p39_link_record(
            *b"LINK  ",
            [*b" C1 ", *b" O1 "],
            [b' ', b' '],
            Some([*b"1555  ", *b"1555  "]),
            None,
            None,
            80,
        );
        let same = decode_pdb_link_record(&same, &model, &mut counters)
            .unwrap()
            .unwrap();
        assert_eq!(same.asu, BioAsu::Same);
        assert_eq!(same.reported_sym, [0; 4]);

        let unparsed_different = p39_link_record(
            *b"LINK  ",
            [*b" C1 ", *b" O1 "],
            [b' ', b' '],
            Some([*b"A     ", *b"BBBBBB"]),
            None,
            None,
            80,
        );
        let unparsed_different = decode_pdb_link_record(&unparsed_different, &model, &mut counters)
            .unwrap()
            .unwrap();
        assert_eq!(unparsed_different.asu, BioAsu::Different);
        assert_eq!(unparsed_different.reported_sym, [0; 4]);

        let too_short = p39_link_record(
            *b"LINK  ",
            [*b" C1 ", *b" O1 "],
            [b' ', b' '],
            None,
            None,
            None,
            56,
        );
        assert_eq!(
            decode_pdb_link_record(&too_short, &model, &mut counters),
            Ok(None)
        );
    }

    #[test]
    fn pdb_link_classification_linkr_and_type_counters_follow_source() {
        let model = p38_model_group(&[]);
        let mut counters = PdbLinkNameCounters::default();

        let covale_record = p39_link_record(
            *b"LINK  ",
            [*b" C1 ", *b" O1 "],
            [b' ', b' '],
            Some([*b"1555  ", *b"1555  "]),
            Some(*b" 1.23"),
            None,
            80,
        );
        let covale = decode_pdb_link_record(&covale_record, &model, &mut counters)
            .unwrap()
            .unwrap();
        assert_eq!(covale.name, "covale1");
        assert_eq!(covale.kind, BioConnectionKind::Covale);
        assert_eq!(covale.reported_distance, 1.23);

        let metal_record = p39_link_record(
            *b"LINK  ",
            [*b"ZN  ", *b" C1 "],
            [b' ', b' '],
            Some([*b"1555  ", *b"1555  "]),
            Some(*b" 2.50"),
            None,
            80,
        );
        let metal = decode_pdb_link_record(&metal_record, &model, &mut counters)
            .unwrap()
            .unwrap();
        assert_eq!(metal.name, "metalc1");
        assert_eq!(metal.kind, BioConnectionKind::MetalC);
        assert_eq!(metal.reported_distance, 2.5);

        let second_covale_record = p39_link_record(
            *b"LINK  ",
            [*b" C1 ", *b" O1 "],
            [b' ', b' '],
            Some([*b"1555  ", *b"1555  "]),
            None,
            None,
            80,
        );
        let second_covale = decode_pdb_link_record(&second_covale_record, &model, &mut counters)
            .unwrap()
            .unwrap();
        assert_eq!(second_covale.name, "covale2");

        let linkr_record = p39_link_record(
            *b"LINKR ",
            [*b" C1 ", *b" O1 "],
            [b' ', b' '],
            Some([*b"1555  ", *b"1555  "]),
            None,
            Some(*b"CHEM1234"),
            80,
        );
        let linkr = decode_pdb_link_record(&linkr_record, &model, &mut counters)
            .unwrap()
            .unwrap();
        assert_eq!(linkr.name, "covale3");
        assert_eq!(linkr.link_id, "CHEM1234");
        assert_eq!(linkr.reported_distance, 0.0);
    }

    #[test]
    fn pdb_link_ambiguous_padded_name_uses_first_model_element_before_fallback() {
        let record = p39_link_record(
            *b"LINK  ",
            [*b"CA11", *b" C1 "],
            [b' ', b' '],
            Some([*b"1555  ", *b"1555  "]),
            None,
            None,
            80,
        );
        let explicit_carbon_model = p39_model_group_with_explicit_carbon_ca11();
        let from_model = decode_pdb_link_record(
            &record,
            &explicit_carbon_model,
            &mut PdbLinkNameCounters::default(),
        )
        .unwrap()
        .unwrap();
        assert_eq!(from_model.kind, BioConnectionKind::Covale);
        assert_eq!(from_model.name, "covale1");

        let no_match_model = p38_model_group(&[]);
        let from_padded_name = decode_pdb_link_record(
            &record,
            &no_match_model,
            &mut PdbLinkNameCounters::default(),
        )
        .unwrap()
        .unwrap();
        assert_eq!(from_padded_name.kind, BioConnectionKind::MetalC);
        assert_eq!(from_padded_name.name, "metalc1");
    }

    #[test]
    fn p40_cispep_single_model_uses_structure_number_and_residue_only_addresses() {
        // Pinned Gemmi 0.7.5 probe: an implicit single model reports model 1;
        // explicit MODEL 7 reports 7 even when CISPEP modNum says 99. Both
        // preserve PHE A 58B / PRO BC 59 with empty atom names and NUL altloc.
        let record = p40_cispep_record(
            *b"CISPEP", *b"A ", *b"  58B", *b"PHE", *b"BC", *b"  59 ", *b"PRO", *b" 99", *b" -0.65",
        );
        let implicit = decode_pdb_cispep_record(&record, &p40_models(&[1]))
            .unwrap()
            .unwrap();
        let explicit = decode_pdb_cispep_record(&record, &p40_models(&[7]))
            .unwrap()
            .unwrap();

        for (cispep, expected_model) in [(&implicit, 1), (&explicit, 7)] {
            assert_eq!(cispep.model_num, expected_model);
            assert_eq!(cispep.partner_c, p40_address(b"A", b"PHE", 58, Some(b'B')));
            assert_eq!(cispep.partner_n, p40_address(b"BC", b"PRO", 59, None));
            assert_eq!(cispep.partner_c.logical_atom_name(), "");
            assert_eq!(cispep.partner_n.logical_atom_name(), "");
            assert_eq!(cispep.partner_c.altloc(), 0);
            assert_eq!(cispep.partner_n.altloc(), 0);
            assert_eq!(cispep.only_altloc, 0);
            assert_eq!(cispep.reported_angle.to_bits(), (-0.65_f64).to_bits());
        }
    }

    #[test]
    fn p40_cispep_multiple_models_read_record_number_and_angle_defaults() {
        // Pinned Gemmi 0.7.5 probe: multiple MODEL rows read modNum; a blank
        // three-byte integer and invalid six-byte angle both retain zero.
        let models = p40_models(&[7, 9]);
        let valid = p40_cispep_record(
            *b"CISPEP", *b"A ", *b"  58B", *b"PHE", *b"BC", *b"  59C", *b"PRO", *b"  4", *b" 12.50",
        );
        let parsed = decode_pdb_cispep_record(&valid, &models).unwrap().unwrap();
        assert_eq!(parsed.model_num, 4);
        assert_eq!(parsed.reported_angle.to_bits(), 12.5_f64.to_bits());
        assert_eq!(parsed.partner_n, p40_address(b"BC", b"PRO", 59, Some(b'C')));

        let source_zeroes = p40_cispep_record(
            *b"cispep", *b"A ", *b"  58 ", *b"PHE", *b"BC", *b"  59C", *b"PRO", *b"   ", *b"abcdef",
        );
        let parsed = decode_pdb_cispep_record(&source_zeroes, &models)
            .unwrap()
            .unwrap();
        assert_eq!(parsed.model_num, 0);
        assert_eq!(parsed.reported_angle.to_bits(), 0.0_f64.to_bits());
        assert_eq!(parsed.only_altloc, 0);
    }

    #[test]
    fn p40_cispep_source_skip_and_undefined_short_field_boundary_are_distinct() {
        // The <22-byte branch is source-defined and skips. For lengths 22..59,
        // Gemmi reads past the std::string character sequence; this assertion
        // verifies only the Rust safety boundary and is not a parity claim.
        let too_short = p40_cispep_record(
            *b"CISPEP", *b"A ", *b"  58 ", *b"PHE", *b"BC", *b"  59 ", *b"PRO", *b"  1", *b"  0.00",
        );
        assert_eq!(decode_pdb_cispep_record(&too_short[..21], &[]), Ok(None));
        assert_eq!(
            decode_pdb_cispep_record(&too_short[..58], &p40_models(&[1])),
            Err(PdbCisPepRecordError::RequiredFieldsPastRecord {
                record_len: 58,
                required_len: 59,
            })
        );
    }

    #[test]
    fn p38_ssbond_uses_periodic_nearest_image_for_different_asu_altlocs() {
        // Pinned Gemmi 0.7.5 read_pdb_string probe: a 10-A P1 cell and
        // Different-ASU SSBOND choose S1/A against S1/A. The Euclidean B/B
        // distance is 0.5 A, while A/A crosses the unit-cell boundary at 0.3 A.
        // Exact input and output are recorded in the P38 Step384 report.
        let model = p38_model_group(&[
            (1, *b" S1 ", b'B', b'A', 1, *b"    ", 0.8, *b" S"),
            (2, *b" CA ", b' ', b'A', 1, *b"    ", 0.0, *b" C"),
            (3, *b" S1 ", b'A', b'A', 1, *b"    ", 0.2, *b" S"),
            (4, *b" S1 ", b'B', b'B', 2, *b"    ", 1.3, *b" S"),
            (5, *b" CA ", b' ', b'B', 2, *b"    ", 0.0, *b" C"),
            (6, *b" S1 ", b'A', b'B', 2, *b"    ", 9.9, *b" S"),
        ]);
        let crystal = p38_prepared_crystal();
        let mut connection = p38_connection(BioAsu::Different);

        complete_ssbond(&mut connection, &model, &crystal).unwrap();

        assert_eq!(connection.asu, BioAsu::Different);
        assert_eq!(connection.partner1.logical_atom_name(), "S1");
        assert_eq!(connection.partner2.logical_atom_name(), "S1");
        assert_eq!(connection.partner1.altloc(), b'A');
        assert_eq!(connection.partner2.altloc(), b'A');
    }

    #[test]
    fn p38_ssbond_equal_distance_keeps_first_source_order_pair() {
        // Pinned Gemmi 0.7.5 read_pdb_string probe: equal Same-ASU distances
        // retain the first B/B pair because complete_ssbond replaces only on
        // strict `dist_sq < min_dist_sq`.
        let model = p38_model_group(&[
            (1, *b" S1 ", b'B', b'A', 1, *b"    ", 0.0, *b" S"),
            (2, *b" S1 ", b'A', b'A', 1, *b"    ", 2.0, *b" S"),
            (3, *b" S1 ", b'B', b'B', 2, *b"    ", 1.0, *b" S"),
            (4, *b" S1 ", b'A', b'B', 2, *b"    ", 3.0, *b" S"),
        ]);
        let crystal = p38_prepared_crystal();
        let mut connection = p38_connection(BioAsu::Same);

        complete_ssbond(&mut connection, &model, &crystal).unwrap();

        assert_eq!(connection.partner1.altloc(), b'B');
        assert_eq!(connection.partner2.altloc(), b'B');
    }

    #[test]
    fn p38_ssbond_falls_back_from_non_sulfur_sg_to_first_sulfur() {
        // Gemmi's complete_ssbond_atom selects the residue's first sulfur if
        // SG is absent or the address-resolved SG atom is not sulfur.
        let model = p38_model_group(&[
            (1, *b" SG ", b' ', b'A', 1, *b"    ", 0.0, *b" C"),
            (2, *b" S1 ", b'B', b'A', 1, *b"    ", 0.5, *b" S"),
            (3, *b" S1 ", b'B', b'B', 2, *b"    ", 0.7, *b" S"),
        ]);
        let crystal = p38_prepared_crystal();
        let mut connection = p38_connection(BioAsu::Same);

        complete_ssbond(&mut connection, &model, &crystal).unwrap();

        for partner in [&connection.partner1, &connection.partner2] {
            assert_eq!(partner.logical_atom_name(), "S1");
            assert_eq!(partner.altloc(), b'B');
        }
    }

    #[test]
    fn p38_ssbond_missing_groups_follow_altloc_gate_and_source_error_order() {
        // Gemmi retains an SSBOND with SG/NUL addresses when both residues lack
        // sulfur and no resolved altloc opens the group-lookup branch.
        let no_sulfur = p38_model_group(&[
            (1, *b" CA ", b' ', b'A', 1, *b"    ", 0.0, *b" C"),
            (2, *b" CA ", b' ', b'B', 2, *b"    ", 1.0, *b" C"),
        ]);
        let crystal = p38_prepared_crystal();
        let mut connection = p38_connection(BioAsu::Any);
        complete_ssbond(&mut connection, &no_sulfur, &crystal).unwrap();
        assert_eq!(connection.partner1.logical_atom_name(), "SG");
        assert_eq!(connection.partner2.logical_atom_name(), "SG");
        assert_eq!(connection.partner1.altloc(), 0);
        assert_eq!(connection.partner2.altloc(), 0);

        // The pinned reader throws `No such atom: SG` for either missing group
        // once the other resolved endpoint has a nonblank altloc. Rust reports
        // partner one first and preserves the caller's value on Err.
        let missing_partner1 = p38_model_group(&[
            (1, *b" CA ", b' ', b'A', 1, *b"    ", 0.0, *b" C"),
            (2, *b" S1 ", b'B', b'B', 2, *b"    ", 1.0, *b" S"),
        ]);
        let before = p38_connection(BioAsu::Any);
        let mut connection = before.clone();
        let error = complete_ssbond(&mut connection, &missing_partner1, &crystal).unwrap_err();
        assert_eq!(
            error,
            PdbSsbondResolutionError {
                partner: 1,
                atom_name: "SG".to_owned(),
            }
        );
        assert_eq!(error.to_string(), "No such atom: SG");
        assert_eq!(connection, before);

        let missing_partner2 = p38_model_group(&[
            (1, *b" S1 ", b'B', b'A', 1, *b"    ", 0.0, *b" S"),
            (2, *b" CA ", b' ', b'B', 2, *b"    ", 1.0, *b" C"),
        ]);
        let before = p38_connection(BioAsu::Any);
        let mut connection = before.clone();
        let error = complete_ssbond(&mut connection, &missing_partner2, &crystal).unwrap_err();
        assert_eq!(
            error,
            PdbSsbondResolutionError {
                partner: 2,
                atom_name: "SG".to_owned(),
            }
        );
        assert_eq!(error.to_string(), "No such atom: SG");
        assert_eq!(connection, before);
    }

    #[test]
    fn p38_ssbond_missing_residue_keeps_connection_address_fallbacks() {
        // Gemmi retains the connection when partner two's residue is absent;
        // partner one still resolves to its first sulfur, and partner two
        // remains the source SG/NUL address.
        let model = p38_model_group(&[(1, *b" S1 ", b'B', b'A', 1, *b"    ", 0.0, *b" S")]);
        let crystal = p38_prepared_crystal();
        let mut connection = p38_connection(BioAsu::Any);

        complete_ssbond(&mut connection, &model, &crystal).unwrap();

        assert_eq!(connection.partner1.logical_atom_name(), "S1");
        assert_eq!(connection.partner1.altloc(), b'B');
        assert_eq!(connection.partner2.logical_atom_name(), "SG");
        assert_eq!(connection.partner2.altloc(), 0);
    }

    #[test]
    fn p38_ssbond_address_requires_the_source_segment_match() {
        // SSBOND's read_res_id initializes segment to empty. Gemmi::find_cra
        // requires exact segment equality, so a same chain/seq/name row with
        // SEG1 is unresolved rather than implicitly matched.
        let model = p38_model_group(&[
            (1, *b" S1 ", b'B', b'A', 1, *b"SEG1", 0.0, *b" S"),
            (2, *b" S1 ", b'B', b'B', 2, *b"    ", 0.5, *b" S"),
        ]);
        let crystal = p38_prepared_crystal();
        let mut connection = p38_connection(BioAsu::Any);

        complete_ssbond(&mut connection, &model, &crystal).unwrap();

        assert_eq!(connection.partner1.logical_atom_name(), "SG");
        assert_eq!(connection.partner1.altloc(), 0);
        assert_eq!(connection.partner2.logical_atom_name(), "S1");
        assert_eq!(connection.partner2.altloc(), b'B');
    }

    #[test]
    fn gemmi_pdb_atom_hierarchy_groups_chain_parts_and_source_residue_keys() {
        // Pinned Gemmi 5cc1c23 probe of populate_structure_from_pdb_stream:
        // source output and full fixed-width input are in IO-bio_pdb.md's P15
        // oracle record. A/B/A becomes three chain parts; A/1A/SEG1/GLY and
        // A/1a/SEG1/GLY reuse the first residue, but segment/name differences
        // append distinct residues. Repeated atoms remain on the first-seen
        // residue in record order; its first ATOM/HETATM record owns het_flag.
        let mut grouping = PdbHierarchyGrouping::default();
        let model_id = grouping.add_model(Some(1)).unwrap();
        let rows = [
            (*b"ATOM  ", 1, *b" CA ", *b"GLY", b'A', 1, b'A', *b"SEG1"),
            (*b"ATOM  ", 2, *b" N  ", *b"GLY", b'A', 2, b' ', *b"SEG1"),
            (*b"ATOM  ", 3, *b" CB ", *b"GLY", b'A', 1, b'a', *b"SEG1"),
            (*b"ATOM  ", 4, *b" C  ", *b"GLY", b'A', 1, b'A', *b"OTHR"),
            (*b"ATOM  ", 5, *b" CA ", *b"ALA", b'A', 1, b'A', *b"SEG1"),
            (*b"ATOM  ", 6, *b" CA ", *b"GLY", b'B', 1, b'A', *b"SEG1"),
            (*b"ATOM  ", 7, *b" CA ", *b"GLY", b'A', 1, b'a', *b"SEG1"),
            (*b"HETATM", 8, *b" O  ", *b"GLY", b'A', 1, b'A', *b"SEG1"),
        ];
        let mut locations = Vec::new();
        for (record, serial, atom, residue, chain, sequence, insertion, segment) in rows {
            let line = grouped_pdb_atom_line(
                record, serial, atom, residue, chain, sequence, insertion, segment,
            );
            locations.push(add_grouped_pdb_atom(&mut grouping, model_id, &line));
        }

        let model = &grouping.models[model_id.index()];
        assert_eq!(model.row_id, model_id);
        assert_eq!(model.source_model_number, Some(1));
        assert_eq!(
            model
                .chains
                .iter()
                .map(|chain| chain.source_id.as_str())
                .collect::<Vec<_>>(),
            ["A", "B", "A"]
        );
        assert_eq!(
            model
                .chains
                .iter()
                .map(|chain| chain.row_id.value())
                .collect::<Vec<_>>(),
            [0, 1, 2]
        );

        let first_chain = &model.chains[0];
        assert_eq!(first_chain.residues.len(), 4);
        assert_eq!(
            first_chain
                .residues
                .iter()
                .map(|residue| residue.row_id.value())
                .collect::<Vec<_>>(),
            [0, 1, 2, 3]
        );
        assert_eq!(
            first_chain
                .residues
                .iter()
                .map(|residue| (
                    residue.address.sequence_number(),
                    residue.address.insertion_code(),
                    residue.address.segment(),
                    residue.address.name(),
                    residue.het_flag,
                ))
                .collect::<Vec<_>>(),
            [
                (
                    Some(1),
                    Some(b'A'),
                    "SEG1",
                    ResidueName::from_ascii(b"GLY").unwrap(),
                    b'A',
                ),
                (
                    Some(2),
                    None,
                    "SEG1",
                    ResidueName::from_ascii(b"GLY").unwrap(),
                    b'A',
                ),
                (
                    Some(1),
                    Some(b'A'),
                    "OTHR",
                    ResidueName::from_ascii(b"GLY").unwrap(),
                    b'A',
                ),
                (
                    Some(1),
                    Some(b'A'),
                    "SEG1",
                    ResidueName::from_ascii(b"ALA").unwrap(),
                    b'A',
                ),
            ]
        );
        assert_eq!(
            first_chain.residues[0]
                .atoms
                .iter()
                .map(|atom| atom.serial.value())
                .collect::<Vec<_>>(),
            [1, 3]
        );
        assert_eq!(
            first_chain.residues[0]
                .atoms
                .iter()
                .map(|atom| {
                    <[u8; 4]>::try_from(atom.name.as_bytes())
                        .expect("PDB atom name has four columns")
                })
                .collect::<Vec<_>>(),
            [*b" CA ", *b" CB "]
        );
        assert_eq!(model.chains[1].residues[0].atoms[0].serial.value(), 6);
        assert_eq!(model.chains[2].residues.len(), 1);
        assert_eq!(
            model.chains[2].residues[0]
                .atoms
                .iter()
                .map(|atom| atom.serial.value())
                .collect::<Vec<_>>(),
            [7, 8]
        );
        assert_eq!(model.chains[2].residues[0].het_flag, b'A');
        assert_eq!(
            locations
                .iter()
                .map(|location| (
                    location.model_id.value(),
                    location.chain_id.value(),
                    location.residue_id.value(),
                    location.atom_index,
                ))
                .collect::<Vec<_>>(),
            [
                (0, 0, 0, 0),
                (0, 0, 1, 0),
                (0, 0, 0, 1),
                (0, 0, 2, 0),
                (0, 0, 3, 0),
                (0, 1, 4, 0),
                (0, 2, 5, 0),
                (0, 2, 5, 1),
            ]
        );
    }

    #[test]
    fn gemmi_pdb_atom_field_decoder_preserves_minimum_xyz_serial_raw_name_and_altloc() {
        let line = pdb_atom_field_line();
        assert!(matches!(
            decode_pdb_atom_field_line(&line[..54]),
            Err(PdbAtomFieldError::TooShort { len: 54 })
        ));

        // Pinned Gemmi 5cc1c23 field oracle: 54 bytes fails; the same row at
        // 55 bytes returns serial 1, trimmed source Atom.name "CA", altloc B,
        // XYZ 1.25/-2.5/3.75, default occupancy/B=1/20, element C, charge 0.
        let fields = decode_pdb_atom_field_line(&line[..55]).unwrap();
        assert_eq!(fields.serial.value(), 1);
        assert_eq!(fields.name.as_bytes(), b" CA ");
        assert_eq!(fields.altloc.map(|value| value.value()), Some(b'B'));
        assert_eq!(fields.position, [1.25, -2.5, 3.75]);
        assert_eq!((fields.occupancy, fields.b_iso), (1.0, 20.0));
        assert_eq!(
            (fields.element, fields.isotope_mass_number),
            (Element::C, None)
        );
        assert_eq!(fields.formal_charge, 0);

        let mut base36 = line;
        base36[6..11].copy_from_slice(b"A0000");
        assert_eq!(
            decode_pdb_atom_field_line(&base36).unwrap().serial.value(),
            100_000
        );
    }

    #[test]
    fn gemmi_pdb_atom_field_decoder_uses_strict_occupancy_and_b_factor_gates() {
        let line = pdb_atom_field_line();
        // Pinned Gemmi 5cc1c23 oracle: lengths 58/64 retain Atom defaults;
        // 59/65 activate the source fields, and B is stored through float.
        assert_eq!(
            decode_pdb_atom_field_line(&line[..58]).unwrap().occupancy,
            1.0
        );
        assert_eq!(
            decode_pdb_atom_field_line(&line[..59]).unwrap().occupancy,
            0.5
        );
        assert_eq!(decode_pdb_atom_field_line(&line[..64]).unwrap().b_iso, 20.0);
        assert_eq!(
            decode_pdb_atom_field_line(&line[..65]).unwrap().b_iso,
            f64::from(12.2_f32)
        );
    }

    #[test]
    fn gemmi_pdb_atom_field_decoder_preserves_element_gate_precedence_and_isotope() {
        let mut line = pdb_atom_field_line();

        // At len 76, the explicit element columns are outside the source gate;
        // the padded atom-name inference still recognizes DXYZ as deuterium.
        line[12..16].copy_from_slice(b"DXYZ");
        line[76..78].copy_from_slice(b" C");
        let inferred = decode_pdb_atom_field_line(&line[..76]).unwrap();
        assert_eq!(
            (inferred.element, inferred.isotope_mass_number),
            (Element::H, Some(2))
        );

        // At len 77, an alphabetic first explicit column is sufficient: the
        // source buffer supplies NUL as the second byte of Element("D\0").
        line[76] = b'D';
        let explicit_one_byte = decode_pdb_atom_field_line(&line[..77]).unwrap();
        assert_eq!(
            (
                explicit_one_byte.element,
                explicit_one_byte.isotope_mass_number
            ),
            (Element::H, Some(2))
        );

        // A selected unknown explicit symbol remains X; it does not fall back
        // to the otherwise recognized D atom-name inference.
        line[76..78].copy_from_slice(b" Q");
        let unknown = decode_pdb_atom_field_line(&line).unwrap();
        assert_eq!(
            (unknown.element, unknown.isotope_mass_number),
            (Element::DUMMY, None)
        );

        // A recognized explicit C overrides the D name, as in the source.
        line[76..78].copy_from_slice(b" C");
        let overridden = decode_pdb_atom_field_line(&line).unwrap();
        assert_eq!(
            (overridden.element, overridden.isotope_mass_number),
            (Element::C, None)
        );
    }

    #[test]
    fn gemmi_pdb_atom_field_decoder_applies_charge_gate_and_propagates_source_errors() {
        let mut line = pdb_atom_field_line();
        assert_eq!(
            decode_pdb_atom_field_line(&line[..78])
                .unwrap()
                .formal_charge,
            0
        );

        // The source gate is len > 78, so its fixed buffer's NUL second byte
        // still permits the one-byte "2" charge at len 79.
        line[78..80].copy_from_slice(b"2+");
        assert_eq!(
            decode_pdb_atom_field_line(&line[..79])
                .unwrap()
                .formal_charge,
            2
        );
        assert_eq!(decode_pdb_atom_field_line(&line).unwrap().formal_charge, 2);

        line[78..80].copy_from_slice(b"+2");
        assert_eq!(decode_pdb_atom_field_line(&line).unwrap().formal_charge, 2);
        line[78..80].copy_from_slice(b"2-");
        assert_eq!(decode_pdb_atom_field_line(&line).unwrap().formal_charge, -2);

        line[78..80].copy_from_slice(b"2x");
        assert!(matches!(
            decode_pdb_atom_field_line(&line),
            Err(PdbAtomFieldError::Charge(PdbChargeError {
                digit: b'2',
                sign: b'x'
            }))
        ));
    }

    #[test]
    fn gemmi_pdb_name_inference_matches_pinned_element_and_isotope_branches() {
        // Pinned Gemmi source probe compiled against third_party/gemmi/src/pdb.cpp:
        //   " H  " -> El::H (enum 1), "1HB " -> El::H (1),
        //   " D  ", "DXYZ", "1D  ", "D1  ", "dXYZ" -> El::D (119),
        //   "CL20"/"cl20" -> Cl (17), "CG11" -> X (0), "DY  " -> Dy (66),
        //   " Q  " -> X (0). El::D's atomic number is 1 while its enum/name stay D.
        let cases = [
            (b" H  ", Element::H, None),
            (b"1HB ", Element::H, None),
            (b" D  ", Element::H, Some(2)),
            (b"DXYZ", Element::H, Some(2)),
            (b"1D  ", Element::H, Some(2)),
            (b"D1  ", Element::H, Some(2)),
            (b"dXYZ", Element::H, Some(2)),
            (b"CL20", Element::CL, None),
            (b"cl20", Element::CL, None),
            (b"CG11", Element::DUMMY, None),
            (b"DY  ", Element::DY, None),
            (b" Q  ", Element::DUMMY, None),
        ];

        for (name, expected_element, expected_isotope) in cases {
            assert_eq!(
                infer_element_from_padded_name(&name),
                (expected_element, expected_isotope),
                "name={name:?}"
            );
        }
    }

    #[test]
    fn gemmi_pdb_explicit_element_precedes_name_only_under_source_gate() {
        // Pinned Gemmi source probe: Element(" D") -> El::D (119),
        // Element(" H") -> El::H (1), Element(" C") -> El::C (6),
        // Element("CL")/Element("Cl") -> El::Cl (17), and Element(" Q") -> El::X.
        let deuterium_name = *b"DXYZ";
        let carbon_name = *b" C  ";
        let mut line = [b' '; 80];

        // With no source explicit-element columns, DXYZ uses Gemmi's D branch.
        assert_eq!(
            resolve_pdb_atom_element(&deuterium_name, &line[..76]),
            (Element::H, Some(2))
        );
        line[76] = b' ';
        assert_eq!(
            resolve_pdb_atom_element(&deuterium_name, &line[..77]),
            (Element::H, Some(2))
        );

        // The caller's len > 76 + ASCII-alpha predicate selects the explicit
        // columns; the explicit value overrides a conflicting atom name.
        line[76..78].copy_from_slice(b" C");
        assert_eq!(
            resolve_pdb_atom_element(&deuterium_name, &line),
            (Element::C, None)
        );
        line[76..78].copy_from_slice(b" D");
        assert_eq!(
            resolve_pdb_atom_element(&carbon_name, &line),
            (Element::H, Some(2))
        );
        line[76..78].copy_from_slice(b"cl");
        assert_eq!(
            resolve_pdb_atom_element(&deuterium_name, &line),
            (Element::CL, None)
        );

        // A selected but unknown explicit symbol remains source X; it does not
        // fall back to the otherwise recognized D atom name.
        line[76..78].copy_from_slice(b" Q");
        assert_eq!(
            resolve_pdb_atom_element(&deuterium_name, &line),
            (Element::DUMMY, None)
        );

        // If column 77 is alphabetic and column 78 is absent, the fixed source
        // buffer contributes its NUL terminator; a one-letter element wins.
        line[76] = b'C';
        assert_eq!(
            resolve_pdb_atom_element(&deuterium_name, &line[..77]),
            (Element::C, None)
        );
    }

    fn pdb_matrix_line(row: u8, values: [&str; 4]) -> [u8; 55] {
        let mut line = [b' '; 55];
        line[5] = row;
        for (start, value) in [10, 20, 30, 45].into_iter().zip(values) {
            let bytes = value.as_bytes();
            assert!(bytes.len() <= 10);
            line[start..start + bytes.len()].copy_from_slice(bytes);
        }
        line
    }

    fn pdb_origx_line(row: u8, values: [&str; 4]) -> [u8; 55] {
        let mut line = pdb_matrix_line(row, values);
        line[..5].copy_from_slice(b"ORIGX");
        line
    }

    fn pdb_scale_line(row: u8, values: [&str; 4]) -> [u8; 55] {
        let mut line = pdb_matrix_line(row, values);
        line[..5].copy_from_slice(b"SCALE");
        line
    }

    fn pdb_mtrix_line(
        row: u8,
        id: &str,
        values: [&str; 4],
        length: usize,
        given: Option<u8>,
    ) -> Vec<u8> {
        assert!(length >= 55);
        assert!(id.len() <= 3);
        let mut line = vec![b' '; length];
        line[..5].copy_from_slice(b"MTRIX");
        line[5] = row;
        line[7..7 + id.len()].copy_from_slice(id.as_bytes());
        for (start, value) in [10, 20, 30, 45].into_iter().zip(values) {
            let bytes = value.as_bytes();
            assert!(bytes.len() <= 10);
            line[start..start + bytes.len()].copy_from_slice(bytes);
        }
        if let Some(given) = given {
            assert!(length > 59);
            line[59] = given;
        }
        line
    }

    fn pdb_cryst1_line(
        a: f64,
        b: f64,
        c: f64,
        alpha: f64,
        beta: f64,
        gamma: f64,
        space_group: &str,
        z_pdb: &str,
    ) -> [u8; 68] {
        let mut line = [b' '; 68];
        line[..6].copy_from_slice(b"CRYST1");
        for (offset, width, value) in [
            (6, 9, a),
            (15, 9, b),
            (24, 9, c),
            (33, 7, alpha),
            (40, 7, beta),
            (47, 7, gamma),
        ] {
            let field = format!("{value:>width$.2}");
            assert_eq!(field.len(), width);
            line[offset..offset + width].copy_from_slice(field.as_bytes());
        }
        assert!(space_group.len() <= 11);
        line[55..55 + space_group.len()].copy_from_slice(space_group.as_bytes());
        assert!(z_pdb.len() <= 4);
        line[66..66 + z_pdb.len()].copy_from_slice(z_pdb.as_bytes());
        line
    }

    #[test]
    fn gemmi_pdb_read_matrix_matches_short_rows_vector_fields_and_completion() {
        // Pinned Gemmi read_matrix probe: short length returns 0 without
        // mutation; suffixes 1, 2, and 3 update only their matching matrix
        // row and vector component; suffix 4 returns 4 without mutation.
        // The source Transform default is identity matrix plus zero vector.
        let identity = BioTransform::identity();
        let mut transform = identity;
        let row1 = pdb_matrix_line(b'1', ["1.25", "2.5", "-3.75", "4.125"]);

        assert_eq!(read_matrix(&mut transform, &[]), 0);
        assert_eq!(transform, identity);
        assert_eq!(read_matrix(&mut transform, &row1[..45]), 0);
        assert_eq!(transform, identity);

        assert_eq!(read_matrix(&mut transform, &row1), 1);
        assert_eq!(
            transform,
            BioTransform::new(
                [[1.25, 2.5, -3.75], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                [4.125, 0.0, 0.0]
            )
        );

        let row4 = pdb_matrix_line(b'4', ["91", "92", "93", "94"]);
        let after_row1 = transform;
        assert_eq!(read_matrix(&mut transform, &row4), 4);
        assert_eq!(transform, after_row1);

        let row2 = pdb_matrix_line(b'2', ["5.25", "-6.5", "7.75", "-8.125"]);
        assert_eq!(read_matrix(&mut transform, &row2), 2);
        assert_eq!(
            transform,
            BioTransform::new(
                [[1.25, 2.5, -3.75], [5.25, -6.5, 7.75], [0.0, 0.0, 1.0]],
                [4.125, -8.125, 0.0]
            )
        );

        let row3 = pdb_matrix_line(b'3', ["9.25", "10.5", "11.75", "12.125"]);
        assert_eq!(read_matrix(&mut transform, &row3), 3);
        assert_eq!(
            transform,
            BioTransform::new(
                [[1.25, 2.5, -3.75], [5.25, -6.5, 7.75], [9.25, 10.5, 11.75]],
                [4.125, -8.125, 12.125]
            )
        );

        // The source's len >= 46 guard permits a one-byte final translation
        // prefix. Its zero-terminated line buffer makes this read as 6.0.
        let mut truncated = [b' '; 46];
        truncated[5] = b'3';
        truncated[10..12].copy_from_slice(b"13");
        truncated[20..22].copy_from_slice(b"14");
        truncated[30..32].copy_from_slice(b"15");
        truncated[45] = b'6';
        let mut partial = identity;
        assert_eq!(read_matrix(&mut partial, &truncated), 3);
        assert_eq!(
            partial,
            BioTransform::new(
                [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [13.0, 14.0, 15.0]],
                [0.0, 0.0, 6.0]
            )
        );
    }

    #[test]
    fn gemmi_pdb_origx_distinguishes_absent_identity_and_present_unread_rows() {
        let identity = BioTransform::identity();
        let absent = PdbReaderState::new("origx.pdb", false, BioStructureSourceState::default());
        assert!(!absent.source_state.has_origx);
        assert_eq!(absent.source_state.origx, identity);

        let mut explicit_identity =
            PdbReaderState::new("origx.pdb", false, BioStructureSourceState::default());
        for row in [
            pdb_origx_line(b'1', ["1", "0", "0", "0"]),
            pdb_origx_line(b'2', ["0", "1", "0", "0"]),
            pdb_origx_line(b'3', ["0", "0", "1", "0"]),
        ] {
            explicit_identity.origx_record(&row);
        }
        assert!(explicit_identity.source_state.has_origx);
        assert_eq!(explicit_identity.source_state.origx, identity);

        // Gemmi sets has_origx before read_matrix's short-row early return.
        let mut short = PdbReaderState::new("origx.pdb", false, BioStructureSourceState::default());
        short.origx_record(b"ORIGX1");
        assert!(short.source_state.has_origx);
        assert_eq!(short.source_state.origx, identity);

        // A complete but out-of-range row label is also present but does not
        // change the Transform, because read_matrix returns that suffix alone.
        let mut invalid_row =
            PdbReaderState::new("origx.pdb", false, BioStructureSourceState::default());
        let row4 = pdb_origx_line(b'4', ["91", "92", "93", "94"]);
        invalid_row.origx_record(&row4);
        assert!(invalid_row.source_state.has_origx);
        assert_eq!(invalid_row.source_state.origx, identity);
    }

    #[test]
    fn gemmi_pdb_origx_updates_partial_rows_and_completes_in_record_order() {
        let mut partial =
            PdbReaderState::new("origx.pdb", false, BioStructureSourceState::default());
        let row2 = pdb_origx_line(b'2', ["2.5", "-3.5", "4.5", "5.5"]);
        partial.origx_record(&row2);
        assert!(partial.source_state.has_origx);
        assert_eq!(
            partial.source_state.origx,
            BioTransform::new(
                [[1.0, 0.0, 0.0], [2.5, -3.5, 4.5], [0.0, 0.0, 1.0]],
                [0.0, 5.5, 0.0]
            )
        );

        let mut complete =
            PdbReaderState::new("origx.pdb", false, BioStructureSourceState::default());
        for row in [
            pdb_origx_line(b'3', ["9.25", "10.5", "11.75", "12.125"]),
            pdb_origx_line(b'1', ["1.25", "2.5", "-3.75", "4.125"]),
            pdb_origx_line(b'2', ["5.25", "-6.5", "7.75", "-8.125"]),
        ] {
            complete.origx_record(&row);
        }
        assert!(complete.source_state.has_origx);
        assert_eq!(
            complete.source_state.origx,
            BioTransform::new(
                [[1.25, 2.5, -3.75], [5.25, -6.5, 7.75], [9.25, 10.5, 11.75]],
                [4.125, -8.125, 12.125]
            )
        );
    }

    #[test]
    fn gemmi_pdb_mtrix_identity_stores_only_the_completed_row_id_in_source_info() {
        let mut state = PdbReaderState::new("mtrix.pdb", false, BioStructureSourceState::default());
        state
            .source_state
            .info
            .insert("_struct_ncs_oper.id".to_owned(), "old".to_owned());

        for row in [
            pdb_mtrix_line(b'1', "A1", ["1", "0", "0", "0"], 80, None),
            pdb_mtrix_line(b'2', "B2", ["0", "1", "0", "0"], 80, None),
            // Identity operations store the row-3 ID before consulting `given`.
            pdb_mtrix_line(b'3', "I3", ["0", "0", "1", "0"], 80, Some(b'1')),
        ] {
            state.mtrix_record(&row).unwrap();
        }

        assert_eq!(
            state
                .source_state
                .info
                .get("_struct_ncs_oper.id")
                .map(String::as_str),
            Some("I3")
        );
        assert!(state.ncs_operators.is_empty());
        assert_eq!(state.matrix, BioTransform::identity());
    }

    #[test]
    fn gemmi_pdb_mtrix_nonidentity_preserves_row3_id_given_order_and_reset() {
        let mut state = PdbReaderState::new("mtrix.pdb", false, BioStructureSourceState::default());
        for row in [
            pdb_mtrix_line(b'1', "A1", ["1", "0", "0", "4"], 80, None),
            pdb_mtrix_line(b'2', "B2", ["0", "2", "0", "5"], 80, None),
            pdb_mtrix_line(b'3', "N7", ["0", "0", "3", "6"], 80, Some(b'1')),
        ] {
            state.mtrix_record(&row).unwrap();
        }
        assert_eq!(state.source_state.info.get("_struct_ncs_oper.id"), None);
        assert_eq!(state.ncs_operators.len(), 1);
        assert_eq!(state.ncs_operators[0].id, "N7");
        assert!(state.ncs_operators[0].given);
        assert_eq!(
            state.ncs_operators[0].transform,
            BioTransform::new(
                [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]],
                [4.0, 5.0, 6.0]
            )
        );
        assert_eq!(state.matrix, BioTransform::identity());

        // A row after completion starts from the source-reset identity scratch.
        let partial = pdb_mtrix_line(b'1', "P1", ["1", "0", "0", "7"], 80, None);
        state.mtrix_record(&partial).unwrap();
        assert_eq!(
            state.matrix,
            BioTransform::new(
                [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                [7.0, 0.0, 0.0]
            )
        );

        for row in [
            pdb_mtrix_line(b'1', "C1", ["1", "0", "0", "4"], 80, None),
            pdb_mtrix_line(b'2', "C2", ["0", "2", "0", "5"], 80, None),
            // With source len == 59, column 60 is absent and `given` is false.
            pdb_mtrix_line(b'3', "N8", ["0", "0", "3", "6"], 59, None),
        ] {
            state.mtrix_record(&row).unwrap();
        }

        assert_eq!(state.ncs_operators.len(), 2);
        assert_eq!(state.ncs_operators[1].id, "N8");
        assert!(!state.ncs_operators[1].given);
        assert_eq!(
            state.ncs_operators[1].transform,
            BioTransform::new(
                [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]],
                [4.0, 5.0, 6.0]
            )
        );
        assert_eq!(state.matrix, BioTransform::identity());
    }

    #[test]
    fn gemmi_pdb_cryst1_record_obeys_source_length_gates_and_trims_fields() {
        // Pinned Gemmi populate_structure_from_pdb_stream uses strict
        // `len > 54`, `len > 56`, and `len > 67` gates. The fixed-width
        // read_string helper left/right trims and stops at EOL/NUL. A pinned
        // Gemmi 5cc1c23 probe produced default cell/no group/no Z at length
        // 54; cell only at 55 and 56; a two-byte group prefix `P` at 57;
        // complete group but no Z at 67; and Z `5` at 68.
        let line = pdb_cryst1_line(10.0, 20.0, 30.0, 90.0, 90.0, 90.0, "P 1", "5");

        let mut state =
            PdbReaderState::new("length.pdb", false, BioStructureSourceState::default());
        let initial = state.crystal_info.clone();
        state.cryst1_record(&line[..54]).unwrap();
        assert_eq!(state.crystal_info, initial);

        for length in [55, 56] {
            let mut state =
                PdbReaderState::new("length.pdb", false, BioStructureSourceState::default());
            state.cryst1_record(&line[..length]).unwrap();
            assert_eq!(
                state.crystal_info.cell(),
                BioCrystalCell {
                    a: 10.0,
                    b: 20.0,
                    c: 30.0,
                    alpha: 90.0,
                    beta: 90.0,
                    gamma: 90.0,
                }
            );
            assert_eq!(state.crystal_info.space_group_hm(), None);
            assert_eq!(state.crystal_info.z_pdb(), None);
        }

        let mut group_prefix =
            PdbReaderState::new("length.pdb", false, BioStructureSourceState::default());
        group_prefix.cryst1_record(&line[..57]).unwrap();
        assert_eq!(group_prefix.crystal_info.space_group_hm(), Some("P"));
        assert_eq!(group_prefix.crystal_info.z_pdb(), None);

        let mut no_z = PdbReaderState::new("length.pdb", false, BioStructureSourceState::default());
        no_z.cryst1_record(&line[..67]).unwrap();
        assert_eq!(no_z.crystal_info.space_group_hm(), Some("P 1"));
        assert_eq!(no_z.crystal_info.z_pdb(), None);

        let mut complete =
            PdbReaderState::new("length.pdb", false, BioStructureSourceState::default());
        complete.cryst1_record(&line).unwrap();
        assert_eq!(complete.crystal_info.space_group_hm(), Some("P 1"));
        assert_eq!(complete.crystal_info.z_pdb(), Some("5"));
    }

    #[test]
    fn gemmi_pdb_cryst1_record_preserves_gamma_zero_and_field_order() {
        let mut state = PdbReaderState::new("order.pdb", false, BioStructureSourceState::default());
        let initial = pdb_cryst1_line(10.0, 20.0, 30.0, 90.0, 90.0, 90.0, "P 1", "8");
        state.cryst1_record(&initial).unwrap();
        let old_cell = state.crystal_info.cell();
        let old_fractional = *state.crystal_info.fractional();
        let old_orthogonal = *state.crystal_info.orthogonal();
        let old_volume = state.crystal_info.volume();

        // Gemmi's UnitCell::set returns without changing any cell-derived
        // state when gamma is zero, but the caller continues with group and Z.
        let partial = pdb_cryst1_line(99.0, 98.0, 97.0, 180.0, 179.0, 0.0, "", "");
        state.cryst1_record(&partial).unwrap();
        assert_eq!(state.crystal_info.cell(), old_cell);
        assert_eq!(*state.crystal_info.fractional(), old_fractional);
        assert_eq!(*state.crystal_info.orthogonal(), old_orthogonal);
        assert_eq!(state.crystal_info.volume(), old_volume);
        assert_eq!(state.crystal_info.space_group_hm(), Some(""));
        assert_eq!(state.crystal_info.z_pdb(), Some("8"));

        let replacement = pdb_cryst1_line(11.0, 21.0, 31.0, 90.0, 90.0, 90.0, "P 2", "12");
        state.cryst1_record(&replacement).unwrap();
        assert_eq!(state.crystal_info.space_group_hm(), Some("P 2"));
        assert_eq!(state.crystal_info.z_pdb(), Some("12"));
    }

    #[test]
    fn gemmi_pdb_cryst1_error_keeps_source_ordered_cell_mutation_only() {
        let mut state = PdbReaderState::new("error.pdb", false, BioStructureSourceState::default());
        let initial = pdb_cryst1_line(10.0, 20.0, 30.0, 90.0, 90.0, 90.0, "P 1", "8");
        state.cryst1_record(&initial).unwrap();
        let old_volume = state.crystal_info.volume();
        let invalid = pdb_cryst1_line(40.0, 50.0, 60.0, 0.0, 90.0, 90.0, "P 2", "9");

        assert_eq!(
            state.cryst1_record(&invalid),
            Err(PdbCrystalError::BioStructureData(
                BioStructureError::ImpossibleCrystalAngle
            ))
        );
        assert_eq!(
            state.crystal_info.cell(),
            BioCrystalCell {
                a: 40.0,
                b: 50.0,
                c: 60.0,
                alpha: 0.0,
                beta: 90.0,
                gamma: 90.0,
            }
        );
        assert_eq!(state.crystal_info.volume(), old_volume);
        assert_eq!(state.crystal_info.space_group_hm(), Some("P 1"));
        assert_eq!(state.crystal_info.z_pdb(), Some("8"));
    }

    #[test]
    fn gemmi_pdb_scale_record_commits_after_row_three_and_resets_staging() {
        let mut state = PdbReaderState::new("scale.pdb", false, BioStructureSourceState::default());
        let initial = state.crystal_info.clone();

        let row1 = pdb_matrix_line(b'1', ["0.5", "0", "0", "0.25"]);
        state.scale_record(&row1);
        assert_eq!(state.crystal_info, initial);
        assert_eq!(
            state.matrix,
            BioTransform::new(
                [[0.5, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                [0.25, 0.0, 0.0]
            )
        );

        let row2 = pdb_matrix_line(b'2', ["0", "0.25", "0", "-0.5"]);
        state.scale_record(&row2);
        assert_eq!(state.crystal_info, initial);

        let row3 = pdb_matrix_line(b'3', ["0", "0", "0.2", "0.75"]);
        state.scale_record(&row3);
        assert_eq!(
            *state.crystal_info.fractional(),
            BioTransform::new(
                [[0.5, 0.0, 0.0], [0.0, 0.25, 0.0], [0.0, 0.0, 0.2]],
                [0.25, -0.5, 0.75]
            )
        );
        assert_eq!(
            *state.crystal_info.orthogonal(),
            BioTransform::new(
                [[2.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 5.0]],
                [-0.5, 2.0, -3.75]
            )
        );
        assert!(state.crystal_info.explicit_matrices());
        assert_eq!(state.matrix, BioTransform::identity());
    }

    #[test]
    fn gemmi_pdb_cryst1_preserves_explicit_scale_matrices_in_source_order() {
        let mut state =
            PdbReaderState::new("precedence.pdb", false, BioStructureSourceState::default());
        let first_cell = pdb_cryst1_line(10.0, 20.0, 30.0, 90.0, 90.0, 90.0, "P 1", "5");
        state.cryst1_record(&first_cell).unwrap();

        state.scale_record(&pdb_matrix_line(b'1', ["0.125", "0", "0", "0.25"]));
        state.scale_record(&pdb_matrix_line(b'2', ["0", "0.0625", "0", "0.5"]));
        state.scale_record(&pdb_matrix_line(b'3', ["0", "0", "0.04166667", "0.75"]));
        let explicit_fractional = *state.crystal_info.fractional();
        let explicit_orthogonal = *state.crystal_info.orthogonal();
        assert!(state.crystal_info.explicit_matrices());

        let later_cell = pdb_cryst1_line(5.0, 6.0, 7.0, 90.0, 90.0, 90.0, "P 2", "9");
        state.cryst1_record(&later_cell).unwrap();
        assert_eq!(
            state.crystal_info.cell(),
            BioCrystalCell {
                a: 5.0,
                b: 6.0,
                c: 7.0,
                alpha: 90.0,
                beta: 90.0,
                gamma: 90.0,
            }
        );
        assert_eq!(*state.crystal_info.fractional(), explicit_fractional);
        assert_eq!(*state.crystal_info.orthogonal(), explicit_orthogonal);
        assert!(state.crystal_info.explicit_matrices());
        assert_eq!(state.crystal_info.space_group_hm(), Some("P 2"));
        assert_eq!(state.crystal_info.z_pdb(), Some("9"));
        assert_eq!(state.crystal_info.volume(), 210.0);
    }

    #[test]
    fn gemmi_pdb_scale_rejection_still_resets_staging_transform() {
        let mut state = PdbReaderState::new(
            "rejected-scale.pdb",
            false,
            BioStructureSourceState::default(),
        );
        let initial = state.crystal_info.clone();

        // UnitCell::set_matrices_from_fract rejects f.mat[0][0] == 0 when the
        // current source cell is the default 1x1x1 non-crystal cell. The PDB
        // caller resets its shared staging matrix after the setter returns,
        // whether or not the setter installs the matrix.
        state.scale_record(&pdb_matrix_line(b'1', ["0", "0", "0", "0.25"]));
        state.scale_record(&pdb_matrix_line(b'2', ["0", "0.5", "0", "0.5"]));
        state.scale_record(&pdb_matrix_line(b'3', ["0", "0", "0.25", "0.75"]));

        assert_eq!(state.crystal_info, initial);
        assert_eq!(state.matrix, BioTransform::identity());
    }

    #[test]
    fn gemmi_pdb_date_format_to_iso_matches_source_positions_and_non_calendar_dates() {
        let cases: &[(&[u8], &str)] = &[
            (b"", ""),
            (b"28-MAR-0", ""),
            (b"x8-MAR-07", ""),
            (b"28-MAR-x7", ""),
            (b"28-MAR-07", "2007-03-28"),
            (b"28-Mar-2007", "2007-03-28"),
            (b"39-FEB-00", "2000-02-39"),
            (b"00-JAN-99", "1999-01-00"),
            (b"28-FOO-07", "2007-xx-28"),
            (b"28-MAR-07x", "2007-03-28"),
            (b"28-MAR-207x", "2020-03-28"),
            (b"28-MAR-20x7", "2020-03-28"),
            (b"28-MAR-2007x", "2007-03-28"),
            (b"28-MAR-07tail", "2007-03-28"),
            (b"28xMAR?07tail", "2007-03-28"),
            (b"28-N01-07", "2007-xx-28"),
            (b"28-n01-07", "2007-xx-28"),
            (b"28-222-07", "2007-xx-28"),
            (b"28-MAR-69", "2069-03-28"),
            (b"28-MAR-70", "1970-03-28"),
            (b"28-\0\0\0-07", "2007-01-28"),
        ];

        for (input, expected) in cases {
            assert_eq!(
                pdb_date_format_to_iso(input),
                *expected,
                "source input={input:?}"
            );
        }
    }

    #[test]
    fn gemmi_pdb_author_name_reorders_initials_with_source_continuation_rules() {
        let cases = [
            ("A.-B.DOE", "DOE, A.-B."),
            ("A.B.DOE", "DOE, A.B."),
            ("JU.DOE", "DOE, JU."),
            ("PON.DOE", "DOE, PON."),
            ("A. DOE", "A. DOE"),
            ("  A.B.DOE", "DOE, A.B."),
            ("\t A.B.DOE", "DOE, \t A.B."),
            ("A.B.C.DOE", "DOE, A.B.C."),
            ("DOE", "DOE"),
            ("", ""),
            ("   ", ""),
        ];

        for (input, expected) in cases {
            let mut actual = input.to_owned();
            change_author_name_format_to_mmcif(&mut actual);
            assert_eq!(actual, expected, "source input={input:?}");
        }
    }

    #[test]
    fn gemmi_pdb_read_int_matches_blank_sign_and_decimal_prefix_rules() {
        for blank in [b"".as_slice(), b"\t\n\x0b\x0c\r "] {
            assert_eq!(read_int(blank), Some(0), "field={blank:?}");
        }

        for no_conversion in [b"+".as_slice(), b"-", b"--5", b"+-5", b"x9"] {
            assert_eq!(read_int(no_conversion), Some(0), "field={no_conversion:?}");
        }

        assert_eq!(read_int(b"+17tail"), Some(17));
        assert_eq!(read_int(b" \t-17rest"), Some(-17));
        assert_eq!(read_int(b"12+3"), Some(12));
        assert_eq!(read_int(b"-0suffix"), Some(0));
    }

    #[test]
    fn gemmi_pdb_read_int_respects_field_width_and_defined_i32_boundaries() {
        let field = b"123456tail";
        assert_eq!(read_int(&field[..5]), Some(12345));
        assert_eq!(read_int(&field[..6]), Some(123456));
        assert_eq!(read_int(b"+"), Some(0));
        assert_eq!(read_int(b"+7"), Some(7));

        assert_eq!(read_int(b"2147483647"), Some(i32::MAX));
        assert_eq!(read_int(b"-2147483648"), Some(i32::MIN));
        assert_eq!(read_int(b"+2147483647"), Some(i32::MAX));
    }

    #[test]
    fn gemmi_pdb_read_int_does_not_invent_values_for_cpp_signed_overflow() {
        // Gemmi's source performs unchecked signed-int arithmetic here, whose
        // result is undefined in C++; `None` is this port's explicit safety
        // boundary, not a claimed Gemmi output.
        assert_eq!(read_int(b"2147483648"), None);
        assert_eq!(read_int(b"-2147483649"), None);
        assert_eq!(read_int(b"+2147483648"), None);
    }

    #[test]
    fn gemmi_pdb_read_double_preserves_initial_zero_and_source_prefix_rules() {
        let cases: &[(&str, &[u8], u64)] = &[
            ("empty", b"", 0x0000_0000_0000_0000),
            ("blank_c_space", b" \t\x0b\x0c\r", 0x0000_0000_0000_0000),
            ("invalid_no_conversion", b"x123", 0x0000_0000_0000_0000),
            ("single_plus_no_conversion", b"+", 0x0000_0000_0000_0000),
            ("second_plus_is_not_removed", b"++1", 0x0000_0000_0000_0000),
            ("one_plus_then_negative", b"+-1", 0xbff0_0000_0000_0000),
            (
                "all_c_space_then_plus_and_prefix",
                b" \t\n\x0b\x0c\r+1.5tail",
                0x3ff8_0000_0000_0000,
            ),
            (
                "positive_numeric_prefix",
                b"12.5tail",
                0x4029_0000_0000_0000,
            ),
            ("negative_numeric_prefix", b"-12.3z", 0xc028_9999_9999_999a),
            ("incomplete_exponent", b"1e", 0x3ff0_0000_0000_0000),
            ("incomplete_exponent_sign", b"1e+", 0x3ff0_0000_0000_0000),
            (
                "complete_exponent_prefix",
                b"1e+2tail",
                0x4059_0000_0000_0000,
            ),
            ("negative_zero_prefix", b"-0tail", 0x8000_0000_0000_0000),
            (
                "negative_sign_without_conversion",
                b"-",
                0x0000_0000_0000_0000,
            ),
            ("hex_is_only_zero_prefix", b"0x1p3", 0x0000_0000_0000_0000),
        ];

        for (case, field, expected_bits) in cases {
            assert_eq!(read_double(field).to_bits(), *expected_bits, "case={case}");
        }
    }

    #[test]
    fn gemmi_pdb_read_double_keeps_fast_float_assignment_when_status_is_ignored() {
        // In the pinned header these overflow/nonzero-underflow conversions
        // set ec=result_out_of_range (34) only after assigning the result.
        // Gemmi read_double discards ec and returns the assigned bits.
        let cases: &[(&str, &[u8], u64)] = &[
            ("positive_overflow", b"1e309", 0x7ff0_0000_0000_0000),
            ("negative_overflow", b"-1e309", 0xfff0_0000_0000_0000),
            ("positive_underflow", b"1e-999", 0x0000_0000_0000_0000),
            (
                "negative_underflow_signed_zero",
                b"-1e-999",
                0x8000_0000_0000_0000,
            ),
            ("representable_subnormal", b"1e-323", 0x0000_0000_0000_0002),
            ("nan_payload", b"nan(payload)", 0x7ff8_0000_0000_0000),
            ("negative_nan", b"-NaN", 0xfff8_0000_0000_0000),
            (
                "nan_invalid_payload_still_matches_prefix",
                b"nan(x$)",
                0x7ff8_0000_0000_0000,
            ),
            (
                "case_insensitive_infinity",
                b"InFiNiTy!",
                0x7ff0_0000_0000_0000,
            ),
            ("inf_matches_prefix", b"infoo", 0x7ff0_0000_0000_0000),
            (
                "ordinary_rounding",
                b"1.234567890123456789",
                0x3ff3_c0ca_428c_59fb,
            ),
            (
                "numeric_prefix_before_nul",
                b"7.25\0tail",
                0x401d_0000_0000_0000,
            ),
        ];

        for (case, field, expected_bits) in cases {
            assert_eq!(read_double(field).to_bits(), *expected_bits, "case={case}");
        }
    }

    #[test]
    fn gemmi_pdb_read_charge_matches_source_sign_order_and_zero_fallback() {
        let cases: &[(&str, u8, u8, i8)] = &[
            ("blank_pair", b' ', b' ', 0),
            ("digit_then_plus", b'2', b'+', 2),
            ("plus_then_digit", b'+', b'2', 2),
            ("digit_then_minus", b'2', b'-', -2),
            ("minus_then_digit", b'-', b'2', -2),
            ("digit_without_sign", b'7', b'\0', 7),
            ("tab_sign", b'4', b'\t', 4),
            ("newline_sign", b'4', b'\n', 4),
            ("vertical_tab_sign", b'4', b'\x0b', 4),
            ("form_feed_sign", b'4', b'\x0c', 4),
            ("carriage_return_sign", b'4', b'\r', 4),
            ("space_sign", b'4', b' ', 4),
            ("zero_with_minus", b'0', b'-', 0),
            // Pinned read_charge accepts each decimal digit, including 9,
            // even though the Atom model comment documents the usual [-8,+8]
            // PDB range.
            ("source_signed_char_negative_nine", b'9', b'-', -9),
            ("nondigit_pair_falls_back", b'x', b'?', 0),
            ("nondigit_with_bad_sign_falls_back", b' ', b'x', 0),
        ];

        for (case, digit, sign, expected) in cases {
            let actual: i8 = read_charge(*digit, *sign).unwrap();
            assert_eq!(actual, *expected, "case={case}");
        }
    }

    #[test]
    fn gemmi_pdb_read_charge_reports_invalid_sign_after_source_digit_swap() {
        let cases: &[(&str, u8, u8, u8, u8)] = &[
            ("bad_sign_after_digit", b'5', b'x', b'5', b'x'),
            ("bad_sign_after_swap", b'x', b'5', b'5', b'x'),
            ("two_digits_swap_then_fail", b'2', b'3', b'3', b'2'),
            ("high_byte_is_not_c_space", b'8', 0x80, b'8', 0x80),
        ];

        for (case, digit, sign, error_digit, error_sign) in cases {
            assert_eq!(
                read_charge(*digit, *sign),
                Err(PdbChargeError {
                    digit: *error_digit,
                    sign: *error_sign,
                }),
                "case={case}"
            );
        }
    }

    #[test]
    fn gemmi_pdb_read_serial_matches_decimal_and_case_insensitive_base36() {
        assert_eq!(read_serial(b"    0"), Some(0));
        assert_eq!(read_serial(b"99999"), Some(99_999));
        assert_eq!(read_serial(b"+1234"), Some(1_234));
        assert_eq!(read_serial(b"-1234"), Some(-1_234));

        assert_eq!(read_serial(b"A0000"), Some(100_000));
        assert_eq!(read_serial(b"A0001"), Some(100_001));
        assert_eq!(read_serial(b"Z0000"), Some(42_090_400));
        assert_eq!(read_serial(b"ZZZZZ"), Some(43_770_015));
        assert_eq!(read_serial(b"a0000"), Some(100_000));
        assert_eq!(read_serial(b"z0000"), Some(42_090_400));
        assert_eq!(read_serial(b"zzzzz"), Some(43_770_015));
    }

    #[test]
    fn gemmi_pdb_read_serial_preserves_source_dispatch_and_strtol_prefix_edges() {
        assert_eq!(read_serial(b"     "), Some(0));
        assert_eq!(read_serial(b"+A001"), Some(0));
        assert_eq!(read_serial(b"A000?"), Some(-16_229_600));
        assert_eq!(read_serial(b"A????"), Some(-16_696_150));
        assert_eq!(read_serial(b"[0000"), Some(-16_696_160));
        assert_eq!(read_serial(b"[????"), Some(-16_696_160));
        assert_eq!(read_serial(b"?0000"), Some(0));

        // The pinned source compiler's plain char is signed: high-byte first
        // bytes take the decimal no-conversion path, not the base36 offset.
        assert_eq!(read_serial(&[0x80, b'0', b'0', b'0', b'0']), Some(0));
        assert_eq!(read_serial(&[0xc3, 0xa9, b'0', b'0', b'0']), Some(0));

        // read_base36 copies all five bytes, then strtol stops at the first NUL.
        assert_eq!(read_serial(&[b'A', b'0', 0, b'9', b'9']), Some(-16_695_800));
    }

    #[test]
    fn gemmi_pdb_read_seq_id_matches_blank_decimal_and_hybrid36_fields() {
        assert_eq!(read_seq_id(b"     "), Some(PdbSeqId::new(i32::MIN, None)));
        assert_eq!(read_seq_id(b"   7A"), Some(PdbSeqId::new(7, Some(b'A'))));
        assert_eq!(read_seq_id(b"-123B"), Some(PdbSeqId::new(-123, Some(b'B'))));
        assert_eq!(read_seq_id(b"+123C"), Some(PdbSeqId::new(123, Some(b'C'))));

        assert_eq!(read_seq_id(b"A000 "), Some(PdbSeqId::new(10_000, None)));
        assert_eq!(
            read_seq_id(b"ZZZZD"),
            Some(PdbSeqId::new(1_223_055, Some(b'D')))
        );
        assert_eq!(
            read_seq_id(b"a000E"),
            Some(PdbSeqId::new(10_000, Some(b'E')))
        );
    }

    #[test]
    fn gemmi_pdb_read_seq_id_preserves_source_dispatch_and_insertion_bytes() {
        // '[' is >= 'A' and selects base36: strtol makes no conversion, then
        // the source's fixed offset is still applied. '?' selects decimal,
        // whose unchecked source parser returns zero for no conversion.
        assert_eq!(
            read_seq_id(b"[000F"),
            Some(PdbSeqId::new(-456_560, Some(b'F')))
        );
        assert_eq!(read_seq_id(b"?000G"), Some(PdbSeqId::new(0, Some(b'G'))));

        assert_eq!(read_seq_id(b"1234\r"), Some(PdbSeqId::new(1234, None)));
        assert_eq!(read_seq_id(b"1234\n"), Some(PdbSeqId::new(1234, None)));
        assert_eq!(
            read_seq_id(&[b'1', b'2', b'3', b'4', 0]),
            Some(PdbSeqId::new(1234, Some(0)))
        );
        assert_eq!(
            read_seq_id(&[0x80, b'1', b'2', b'3', b'H']),
            Some(PdbSeqId::new(0, Some(b'H')))
        );
    }

    #[test]
    fn gemmi_pdb_header_info_and_author_continuations_match_pinned_source() {
        let mut state = PdbReaderState::new("p24.pdb", false, BioStructureSourceState::default());

        let mut header = [b' '; 122];
        header[..6].copy_from_slice(b"HEADER");
        header[10..21].copy_from_slice(b"  CLASS    ");
        header[50..59].copy_from_slice(b"01-JAN-00");
        header[62..66].copy_from_slice(b"AB  ");
        header[80] = b'\n';
        header[81] = 0;
        assert!(state.header_record(&header, 81).unwrap());

        let record = |kind: [u8; 6], payload: &[u8]| {
            let mut line = [b' '; 122];
            line[..6].copy_from_slice(&kind);
            let end = 10 + payload.len();
            line[10..end].copy_from_slice(payload);
            line[end] = b'\n';
            line[end + 1] = 0;
            (line, end + 1)
        };
        for (kind, payload) in [
            (*b"TITLE ", &b"FIRST  "[..]),
            (*b"title ", &b"SECOND\t "[..]),
            (*b"KEYWDS", &b"ALPHA  "[..]),
            (*b"KEYWDS", &b" BETA"[..]),
            (*b"EXPDTA", &b"  X-RAY "[..]),
            (*b"EXPDTA", &b" EM \t"[..]),
            (*b"AUTHOR", &b"Doe, A.-"[..]),
            (*b"AUTHOR", &b"B.Smith, C.Jones"[..]),
        ] {
            let (line, line_len) = record(kind, payload);
            assert!(state.header_record(&line, line_len).unwrap(), "{kind:?}");
        }

        assert_eq!(
            state
                .source_state
                .info
                .get("_struct_keywords.pdbx_keywords")
                .map(String::as_str),
            Some("  CLASS")
        );
        assert_eq!(
            state
                .source_state
                .info
                .get("_pdbx_database_status.recvd_initial_deposition_date")
                .map(String::as_str),
            Some("2000-01-01")
        );
        assert_eq!(
            state.source_state.info.get("_entry.id").map(String::as_str),
            Some("AB")
        );
        assert_eq!(
            state
                .source_state
                .info
                .get("_struct.title")
                .map(String::as_str),
            Some("FIRSTSECOND")
        );
        assert_eq!(
            state
                .source_state
                .info
                .get("_struct_keywords.text")
                .map(String::as_str),
            Some("ALPHA BETA")
        );
        assert_eq!(
            state
                .source_state
                .info
                .get("_exptl.method")
                .map(String::as_str),
            Some("X-RAYEM")
        );

        state.finalize_author_names();
        assert_eq!(
            state.metadata.authors,
            vec![
                "Doe".to_owned(),
                "Smith, A.-B.".to_owned(),
                "Jones, C.".to_owned()
            ]
        );
    }

    #[test]
    fn gemmi_pdb_header_length_gates_and_invalid_utf8_are_explicit() {
        let mut state = PdbReaderState::new("p24.pdb", false, BioStructureSourceState::default());
        let mut full_header = [b' '; 122];
        full_header[..6].copy_from_slice(b"HEADER");
        full_header[10..15].copy_from_slice(b"CLASS");
        full_header[50..59].copy_from_slice(b"01-JAN-00");
        full_header[62..66].copy_from_slice(b"AB  ");
        full_header[80] = b'\n';
        full_header[81] = 0;
        assert!(state.header_record(&full_header, 81).unwrap());

        // Gemmi writes the keyword field at len > 50, while date and entry
        // remain unchanged when their independent physical-length gates close.
        let mut short_header = [b' '; 122];
        short_header[..6].copy_from_slice(b"HEADER");
        short_header[50] = b'\n';
        short_header[51] = 0;
        assert!(state.header_record(&short_header, 51).unwrap());
        assert_eq!(
            state
                .source_state
                .info
                .get("_struct_keywords.pdbx_keywords")
                .map(String::as_str),
            Some("")
        );
        assert_eq!(
            state
                .source_state
                .info
                .get("_pdbx_database_status.recvd_initial_deposition_date")
                .map(String::as_str),
            Some("2000-01-01")
        );
        assert_eq!(
            state.source_state.info.get("_entry.id").map(String::as_str),
            Some("AB")
        );

        let mut invalid_title = [b' '; 122];
        invalid_title[..6].copy_from_slice(b"TITLE ");
        invalid_title[10] = 0x80;
        invalid_title[11] = b'\n';
        invalid_title[12] = 0;
        let before = state.source_state.info.clone();
        assert_eq!(
            state.header_record(&invalid_title, 12),
            Err(PdbHeaderError::TextFieldNotUtf8 {
                field_offset: 10,
                valid_up_to: 0
            })
        );
        assert_eq!(state.source_state.info, before);

        let (line, line_len) = {
            let mut line = [b' '; 122];
            line[..6].copy_from_slice(b"REMARK");
            line[6] = b'\n';
            line[7] = 0;
            (line, 7)
        };
        assert!(!state.header_record(&line, line_len).unwrap());
    }

    #[test]
    fn gemmi_pdb_raw_remarks_preserve_order_and_remove_only_one_source_terminator() {
        // Pinned Gemmi Step306 oracle retained the same raw rows with
        // skip_remarks=false and true. The Rust capture helper has no metadata
        // interpretation flag, matching its earlier source-dispatch stage.
        let input = concat!(
            "remark lower lf\n",
            "REMARK crlf\r\n",
            "REMARK double-crlf\r\r\n",
            "REMARK repeated\n",
            "REMARK repeated\n",
            "OTHER not-a-remark\n",
        );
        let mut cursor = PdbLineCursor::new(input, PdbInputOptions::default());
        let mut state = PdbReaderState::new("p25.pdb", false, BioStructureSourceState::default());
        let mut recognized = Vec::new();
        while let Some(line) = cursor.copy_line() {
            let line_len = line.len();
            let mut source_line_buffer = [0; 122];
            source_line_buffer[..line_len].copy_from_slice(line);
            recognized.push(state.remark_record(&source_line_buffer, line_len).unwrap());
        }
        assert_eq!(recognized, [true, true, true, true, true, false]);
        assert_eq!(
            state.source_state.raw_remarks,
            [
                "remark lower lf",
                "REMARK crlf",
                "REMARK double-crlf\r",
                "REMARK repeated",
                "REMARK repeated",
            ]
        );

        let capture_single_eof_record = |input: &str| {
            let mut cursor = PdbLineCursor::new(input, PdbInputOptions::default());
            let mut state =
                PdbReaderState::new("p25-eof.pdb", false, BioStructureSourceState::default());
            while let Some(line) = cursor.copy_line() {
                let line_len = line.len();
                let mut source_line_buffer = [0; 122];
                source_line_buffer[..line_len].copy_from_slice(line);
                state.remark_record(&source_line_buffer, line_len).unwrap();
            }
            state.source_state.raw_remarks
        };
        assert_eq!(
            capture_single_eof_record("REMARK no-terminator"),
            ["REMARK no-terminator"]
        );
        assert_eq!(
            capture_single_eof_record("REMARK cr-only\r"),
            ["REMARK cr-only"]
        );
    }

    #[test]
    fn gemmi_pdb_raw_remark_representation_errors_are_typed_and_atomic() {
        let mut state = PdbReaderState::new("p25.pdb", false, BioStructureSourceState::default());
        let mut invalid_utf8 = [0; 122];
        invalid_utf8[..8].copy_from_slice(b"REMARK \xff");
        assert_eq!(
            state.remark_record(&invalid_utf8, 8),
            Err(PdbRemarkError::TextFieldNotUtf8 { valid_up_to: 7 })
        );
        assert!(state.source_state.raw_remarks.is_empty());

        assert_eq!(
            state.remark_record(&invalid_utf8, 123),
            Err(PdbRemarkError::LineLengthOutsideBuffer { line_len: 123 })
        );
        assert!(state.source_state.raw_remarks.is_empty());
    }

    #[test]
    fn gemmi_pdb_remark3_refinement_scalar_fields_match_pinned_source() {
        let mut state = PdbReaderState::new("p26.pdb", false, BioStructureSourceState::default());
        feed_remark3(&mut state, 1, "PROGRAM: PRE-REFINEMENT 1.0");
        assert!(state.metadata.refinement.is_empty());
        assert_eq!(state.metadata.software[0].name, "PRE-REFINEMENT");
        assert_eq!(state.metadata.software[0].version, "1.0");
        feed_remark3(&mut state, 1, "DATA USED IN REFINEMENT.");

        // A per-bin scalar is ignored until the source has created a bin.
        feed_remark3(&mut state, 1, "BIN RESOLUTION RANGE HIGH       (A): 1.6");
        assert!(state.metadata.refinement[0].bins.is_empty());

        feed_remark3(
            &mut state,
            1,
            "PROGRAM: autoPROC (Version 1.3.0), AIMLESS, STARANISO (VERSION Jan 26, 2018)",
        );
        feed_remark3(&mut state, 1, "RESOLUTION RANGE HIGH (ANGSTROMS): 2.1");
        feed_remark3(&mut state, 1, "RESOLUTION RANGE LOW  (ANGSTROMS): 0.8");
        feed_remark3(&mut state, 1, "COMPLETENESS FOR RANGE        (%): 99.1");
        feed_remark3(
            &mut state,
            1,
            "NUMBER OF REFLECTIONS: \u{000b}+123 trailing",
        );
        feed_remark3(&mut state, 1, "CROSS-VALIDATION METHOD: Leave-One-Out");
        feed_remark3(&mut state, 1, "FREE R VALUE TEST SET SELECTION: FREERFLAG");
        feed_remark3(&mut state, 1, "R VALUE     (WORKING + TEST SET): 0.2");
        feed_remark3(&mut state, 1, "R VALUE            (WORKING SET): 0.19");
        feed_remark3(&mut state, 1, "FREE R VALUE: 0.24");
        feed_remark3(&mut state, 1, "FREE R VALUE TEST SET COUNT: -55 trailing");
        feed_remark3(&mut state, 1, "TOTAL NUMBER OF BINS USED: no-digits");
        feed_remark3(&mut state, 1, "FIT IN THE HIGHEST RESOLUTION BIN.");
        feed_remark3(&mut state, 1, "BIN RESOLUTION RANGE HIGH       (A): 2.0");
        feed_remark3(&mut state, 1, "BIN RESOLUTION RANGE LOW        (A): 0.8");
        feed_remark3(&mut state, 1, "BIN COMPLETENESS (WORKING+TEST) (%): 100.0");
        feed_remark3(&mut state, 1, "REFLECTIONS IN BIN   (WORKING+TEST): 456");
        feed_remark3(&mut state, 1, "BIN R VALUE          (WORKING+TEST): 0.21");
        feed_remark3(&mut state, 1, "REFLECTIONS IN BIN    (WORKING SET): 400");
        feed_remark3(&mut state, 1, "BIN R VALUE           (WORKING SET): 0.2");
        feed_remark3(&mut state, 1, "BIN FREE R VALUE: 0.25");
        feed_remark3(&mut state, 1, "BIN FREE R VALUE TEST SET COUNT: 56");
        feed_remark3(&mut state, 1, "FROM WILSON PLOT           (A**2): 77.7");
        feed_remark3(&mut state, 1, "MEAN B VALUE      (OVERALL, A**2): 12.5");
        feed_remark3(&mut state, 1, "B11 (A**2): 11");
        feed_remark3(&mut state, 1, "B22 (A**2): 22");
        feed_remark3(&mut state, 1, "B33 (A**2): 33");
        feed_remark3(&mut state, 1, "B12 (A**2): 12");
        feed_remark3(&mut state, 1, "B13 (A**2): 13");
        feed_remark3(&mut state, 1, "B23 (A**2): 23");
        feed_remark3(
            &mut state,
            1,
            "ESD FROM LUZZATI PLOT                    (A): 0.22",
        );
        feed_remark3(
            &mut state,
            1,
            "DPI (BLOW EQ-10) BASED ON R VALUE        (A): 0.31",
        );
        feed_remark3(
            &mut state,
            1,
            "DPI (BLOW EQ-9) BASED ON FREE R VALUE    (A): 0.32",
        );
        feed_remark3(
            &mut state,
            1,
            "DPI (CRUICKSHANK) BASED ON R VALUE       (A): 0.33",
        );
        feed_remark3(
            &mut state,
            1,
            "DPI (CRUICKSHANK) BASED ON FREE R VALUE  (A): 0.34",
        );
        feed_remark3(&mut state, 1, "CORRELATION COEFFICIENT FO-FC: 0.91");
        feed_remark3(&mut state, 1, "CORRELATION COEFFICIENT FO-FC FREE: 0.92");

        assert_eq!(state.metadata.refinement.len(), 1);
        assert_eq!(state.metadata.refinement[0].id, "1");
        let refinement = &state.metadata.refinement[0];
        assert_eq!(refinement.basic.resolution_high, 2.1);
        assert_eq!(refinement.basic.resolution_low, 0.8);
        assert_eq!(refinement.basic.completeness, 99.1);
        assert_eq!(refinement.basic.reflection_count, 123);
        assert_eq!(refinement.cross_validation_method, "Leave-One-Out");
        assert_eq!(refinement.rfree_selection_method, "FREERFLAG");
        assert_eq!(refinement.basic.r_all, 0.2);
        assert_eq!(refinement.basic.r_work, 0.19);
        assert_eq!(refinement.basic.r_free, 0.24);
        assert_eq!(refinement.basic.rfree_set_count, -55);
        assert_eq!(refinement.bin_count, 0);
        assert_eq!(refinement.mean_b, 12.5);
        assert_eq!(refinement.aniso_b, [11.0, 22.0, 33.0, 12.0, 13.0, 23.0]);
        assert_eq!(refinement.luzzati_error, 0.22);
        assert_eq!(refinement.dpi_blow_r, 0.31);
        assert_eq!(refinement.dpi_blow_rfree, 0.32);
        assert_eq!(refinement.dpi_cruickshank_r, 0.33);
        assert_eq!(refinement.dpi_cruickshank_rfree, 0.34);
        assert_eq!(refinement.basic.cc_fo_fc_work, 0.91);
        assert_eq!(refinement.basic.cc_fo_fc_free, 0.92);
        assert_eq!(refinement.bins.len(), 1);
        assert_eq!(refinement.bins[0].resolution_high, 2.0);
        assert_eq!(refinement.bins[0].resolution_low, 0.8);
        assert_eq!(refinement.bins[0].completeness, 100.0);
        assert_eq!(refinement.bins[0].reflection_count, 456);
        assert_eq!(refinement.bins[0].work_set_count, 400);
        assert_eq!(refinement.bins[0].rfree_set_count, 56);
        assert_eq!(refinement.bins[0].r_all, 0.21);
        assert_eq!(refinement.bins[0].r_work, 0.2);
        assert_eq!(refinement.bins[0].r_free, 0.25);
        assert_eq!(state.source_state.resolution, 2.1);
        assert!(state.metadata.experiments.is_empty());

        assert_eq!(state.metadata.software.len(), 4);
        assert_eq!(state.metadata.software[0].name, "PRE-REFINEMENT");
        assert_eq!(state.metadata.software[0].version, "1.0");
        assert_eq!(
            state.metadata.software[0].classification,
            BioSoftwareClassification::Refinement
        );
        assert_eq!(state.metadata.software[1].name, "autoPROC");
        assert_eq!(state.metadata.software[1].version, "1.3.0");
        assert_eq!(
            state.metadata.software[1].classification,
            BioSoftwareClassification::Refinement
        );
        assert_eq!(state.metadata.software[2].name, "AIMLESS");
        assert_eq!(state.metadata.software[2].version, "");
        assert_eq!(state.metadata.software[3].name, "STARANISO");
        assert_eq!(state.metadata.software[3].version, "Jan 26, 2018");
    }

    #[test]
    fn gemmi_pdb_remark3_software_segments_dates_null_and_order_match_pinned_source() {
        let mut state = PdbReaderState::new(
            "p29-software.pdb",
            false,
            BioStructureSourceState::default(),
        );
        feed_remark3(&mut state, 0, "DATA USED IN REFINEMENT.");
        feed_remark3(&mut state, 0, "PROGRAM: ALPHA, 1.25, BETA (Version 2.0)");
        feed_remark3(&mut state, 0, "PROGRAM: BETA 3, DELTA 1.0 (28-MAR-07)");
        feed_remark3(&mut state, 0, "PROGRAM: EPSILON VERSION 1.0 (28-Mar-2007)");
        feed_remark3(&mut state, 0, "PROGRAM: ZETA 1.0 (28-FOO-07), ETA 1.5)");
        feed_remark3(&mut state, 0, "PROGRAM: NO_VERSION (");
        feed_remark3(&mut state, 0, "PROGRAM: DUPLICATE 1");
        feed_remark3(&mut state, 0, "PROGRAM: DUPLICATE 1");
        feed_remark3(&mut state, 0, "PROGRAM: NULL");
        feed_remark3(&mut state, 0, "PROGRAM:");

        let software = &state.metadata.software;
        assert_eq!(software.len(), 11);
        assert_eq!(
            software
                .iter()
                .map(|item| item.name.as_str())
                .collect::<Vec<_>>(),
            [
                "ALPHA,",
                "BETA",
                "BETA",
                "DELTA",
                "EPSILON",
                "ZETA",
                "ETA",
                "NO_VERSION",
                "DUPLICATE",
                "DUPLICATE",
                "",
            ]
        );
        assert_eq!(
            software
                .iter()
                .map(|item| item.version.as_str())
                .collect::<Vec<_>>(),
            [
                "1.25",
                "2.0",
                "3",
                "1.0",
                "1.0",
                "1.0 (28-FOO-07)",
                "1.5",
                "",
                "1",
                "1",
                "",
            ]
        );
        assert_eq!(
            software
                .iter()
                .map(|item| item.date.as_str())
                .collect::<Vec<_>>(),
            [
                "",
                "",
                "",
                "2007-03-28",
                "2007-03-28",
                "",
                "",
                "",
                "",
                "",
                "",
            ]
        );
        assert!(software.iter().all(|item| {
            item.classification == BioSoftwareClassification::Refinement
                && item.description.is_empty()
                && item.contact_author.is_empty()
                && item.contact_author_email.is_empty()
        }));
    }

    #[test]
    fn gemmi_pdb_remark3_program_obeys_tls_continuation_before_dispatch() {
        let mut state = PdbReaderState::new(
            "p29-continuation.pdb",
            false,
            BioStructureSourceState::default(),
        );
        feed_remark3(&mut state, 0, "DATA USED IN REFINEMENT.");
        feed_remark3(&mut state, 0, "TLS GROUP: 1");
        feed_remark3(&mut state, 2, "SET: CHAIN A");
        assert_eq!(
            state.metadata.refinement[0].tls_groups[0].selections[0].details,
            "CHAIN A"
        );

        feed_remark3(&mut state, 7, "PROGRAM: SHOULD_NOT_BE_A_SOFTWARE_ROW");
        assert_eq!(
            state.metadata.refinement[0].tls_groups[0].selections[0].details,
            "CHAIN A PROGRAM"
        );
        assert!(state.metadata.software.is_empty());
        assert!(state.remark3_continuation.is_some());

        feed_remark3(&mut state, 0, "PROGRAM: RELEASED 4.2");
        assert_eq!(
            state.metadata.refinement[0].tls_groups[0].selections[0].details,
            "CHAIN A PROGRAM"
        );
        assert!(state.remark3_continuation.is_none());
        assert_eq!(state.metadata.software.len(), 1);
        assert_eq!(state.metadata.software[0].name, "RELEASED");
        assert_eq!(state.metadata.software[0].version, "4.2");
        assert_eq!(
            state.metadata.software[0].classification,
            BioSoftwareClassification::Refinement
        );
    }

    #[test]
    fn gemmi_pdb_remark3_restraint_key_mapping_and_order_match_pinned_source() {
        let mut state = PdbReaderState::new("p28.pdb", false, BioStructureSourceState::default());
        feed_remark3(&mut state, 1, "DATA USED IN REFINEMENT.");

        let source_rows = [
            (
                "BOND LENGTHS",
                "17 ; 0.25 ; HARMONIC",
                "t_bond_d",
                17,
                0.25,
                "HARMONIC",
            ),
            (
                "BOND ANGLES",
                "18 ; -0.5 ; ANGLE",
                "t_angle_deg",
                18,
                -0.5,
                "ANGLE",
            ),
            (
                "TORSION ANGLES",
                "19 ; 1.25 ; TORSION",
                "t_dihedral_angle_d",
                19,
                1.25,
                "TORSION",
            ),
            (
                "TRIGONAL CARBON PLANES",
                "20 ; 2 ; TRIGONAL",
                "t_trig_c_planes",
                20,
                2.0,
                "TRIGONAL",
            ),
            (
                "GENERAL PLANES",
                "21 ; 3 ; GENERAL",
                "t_gen_planes",
                21,
                3.0,
                "GENERAL",
            ),
            (
                "ISOTROPIC THERMAL FACTORS",
                "22 ; 4 ; ISOTROPIC",
                "t_it",
                22,
                4.0,
                "ISOTROPIC",
            ),
            (
                "BAD NON-BONDED CONTACTS",
                "23 ; 5 ; CONTACT",
                "t_nbd",
                23,
                5.0,
                "CONTACT",
            ),
            (
                "IMPROPER TORSIONS",
                "24 ; 6 ; IMPROPER",
                "t_improper_torsion",
                24,
                6.0,
                "IMPROPER",
            ),
            (
                "CHIRAL IMPROPER TORSION",
                "25 ; 7 ; CHIRAL",
                "t_chiral_improper_torsion",
                25,
                7.0,
                "CHIRAL",
            ),
            (
                "SUM OF OCCUPANCIES",
                "26 ; 8 ; OCCUPANCY",
                "t_sum_occupancies",
                26,
                8.0,
                "OCCUPANCY",
            ),
            (
                "UTILITY DISTANCES",
                "27 ; 9 ; DISTANCE",
                "t_utility_distance",
                27,
                9.0,
                "DISTANCE",
            ),
            (
                "UTILITY ANGLES",
                "28 ; 10 ; ANGLE",
                "t_utility_angle",
                28,
                10.0,
                "ANGLE",
            ),
            (
                "UTILITY TORSION",
                "29 ; 11 ; TORSION",
                "t_utility_torsion",
                29,
                11.0,
                "TORSION",
            ),
            (
                "IDEAL-DIST CONTACT TERM",
                "30 ; 12 ; IDEAL",
                "t_ideal_dist_contact",
                30,
                12.0,
                "IDEAL",
            ),
        ];

        for &(key, value, _, _, _, _) in &source_rows {
            feed_remark3(&mut state, 1, &format!("{key}: {value}"));
        }

        let restraints = &state.metadata.refinement[0].restr_stats;
        assert_eq!(restraints.len(), source_rows.len());
        for (actual, (_, _, name, count, weight, function)) in restraints.iter().zip(source_rows) {
            assert_eq!(actual.name, name);
            assert_eq!(actual.count, count);
            assert_eq!(actual.weight, weight);
            assert_eq!(actual.function, function);
            assert!(actual.dev_ideal.is_nan());
        }
    }

    #[test]
    fn gemmi_pdb_remark3_restraint_cursor_defaults_and_repeats_match_pinned_source() {
        let mut state =
            PdbReaderState::new("p28-cursor.pdb", false, BioStructureSourceState::default());
        feed_remark3(&mut state, 1, "DATA USED IN REFINEMENT.");
        feed_remark3(&mut state, 1, "BOND LENGTHS: 31 ;");
        feed_remark3(&mut state, 1, "BOND ANGLES: 32");
        feed_remark3(&mut state, 1, "TORSION ANGLES: 33; bad ; BAD-WEIGHT");
        feed_remark3(
            &mut state,
            1,
            "TRIGONAL CARBON PLANES: +34 ; 0.75 ; +SPRING",
        );
        feed_remark3(&mut state, 1, "GENERAL PLANES: NONSENSE");
        feed_remark3(&mut state, 1, "ISOTROPIC THERMAL FACTORS: null");
        feed_remark3(&mut state, 1, "BAD NON-BONDED CONTACTS: NULL");
        feed_remark3(&mut state, 1, "BOND LENGTHS: 35 ; 0.7 ; REPEAT-ONE");
        feed_remark3(&mut state, 1, "BOND LENGTHS: 36 ; 0.8 ; REPEAT-TWO");

        let restraints = &state.metadata.refinement[0].restr_stats;
        assert_eq!(restraints.len(), 7);

        assert_eq!(restraints[0].name, "t_bond_d");
        assert_eq!(restraints[0].count, 31);
        assert_eq!(restraints[0].weight, 0.0);
        assert!(restraints[0].function.is_empty());

        assert_eq!(restraints[1].name, "t_angle_deg");
        assert_eq!(restraints[1].count, 32);
        assert!(restraints[1].weight.is_nan());
        assert!(restraints[1].function.is_empty());

        assert_eq!(restraints[2].name, "t_dihedral_angle_d");
        assert_eq!(restraints[2].count, 33);
        assert_eq!(restraints[2].weight, 0.0);
        assert_eq!(restraints[2].function, "BAD-WEIGHT");

        assert_eq!(restraints[3].name, "t_trig_c_planes");
        assert_eq!(restraints[3].count, 0);
        assert_eq!(restraints[3].weight, 0.75);
        assert_eq!(restraints[3].function, "+SPRING");

        assert_eq!(restraints[4].name, "t_it");
        assert_eq!(restraints[4].count, 0);
        assert!(restraints[4].weight.is_nan());
        assert!(restraints[4].function.is_empty());

        assert_eq!(restraints[5].name, "t_bond_d");
        assert_eq!(restraints[5].count, 35);
        assert_eq!(restraints[5].weight, 0.7);
        assert_eq!(restraints[5].function, "REPEAT-ONE");
        assert_eq!(restraints[6].name, "t_bond_d");
        assert_eq!(restraints[6].count, 36);
        assert_eq!(restraints[6].weight, 0.8);
        assert_eq!(restraints[6].function, "REPEAT-TWO");
    }

    #[test]
    fn gemmi_pdb_remark3_restraint_function_uses_source_50_byte_field() {
        let mut state = PdbReaderState::new(
            "p28-function.pdb",
            false,
            BioStructureSourceState::default(),
        );
        feed_remark3(&mut state, 1, "DATA USED IN REFINEMENT.");
        let function = "X".repeat(60);
        feed_remark3(
            &mut state,
            1,
            &format!("GENERAL PLANES: 7 ; 0.5 ; {function}"),
        );

        let restraints = &state.metadata.refinement[0].restr_stats;
        assert_eq!(restraints.len(), 1);
        assert_eq!(restraints[0].name, "t_gen_planes");
        assert_eq!(restraints[0].count, 7);
        assert_eq!(restraints[0].weight, 0.5);
        assert_eq!(restraints[0].function, "X".repeat(49));
    }

    #[test]
    fn gemmi_pdb_remark3_null_defaults_continuation_and_resolution_order_match_source() {
        let mut state =
            PdbReaderState::new("p26-order.pdb", false, BioStructureSourceState::default());
        feed_remark3(&mut state, 1, "DATA USED IN REFINEMENT.");
        feed_remark3(&mut state, 1, "RESOLUTION RANGE HIGH (ANGSTROMS): 0");
        feed_remark3(&mut state, 1, "FREE R VALUE: 0.6");
        feed_remark3(&mut state, 1, "FREE R VALUE: NULL");
        feed_remark3(&mut state, 1, "FREE R VALUE: null");
        feed_remark3(&mut state, 1, "data used in refinement.");
        assert_eq!(state.metadata.refinement.len(), 1);
        feed_remark3(&mut state, 1, "DATA USED IN REFINEMENT.");
        feed_remark3(&mut state, 1, "RESOLUTION RANGE HIGH (ANGSTROMS): nan");
        assert!(state.metadata.refinement[1].basic.resolution_high.is_nan());
        assert_eq!(state.metadata.refinement[1].basic.reflection_count, -1);
        assert_eq!(state.metadata.refinement[1].bin_count, -1);
        feed_remark3(&mut state, 1, "DATA USED IN REFINEMENT.");
        feed_remark3(&mut state, 1, "RESOLUTION RANGE HIGH (ANGSTROMS): 3.3");
        feed_remark3(&mut state, 1, "RESOLUTION RANGE HIGH (ANGSTROMS): 4.4");

        let mut group = BioTlsGroup {
            id: "1".to_owned(),
            num_id: 1,
            ..BioTlsGroup::default()
        };
        group.selections.push(BioTlsSelection {
            details: "CHAIN A".to_owned(),
            ..BioTlsSelection::default()
        });
        state.metadata.refinement[2].tls_groups.push(group);
        state.remark3_continuation = Some(PdbRemark3Continuation {
            refinement_index: 2,
            tls_group_index: 0,
            selection_index: 0,
        });
        feed_remark3(&mut state, 8, "AND CHAIN B");
        assert_eq!(
            state.metadata.refinement[2].tls_groups[0].selections[0].details,
            "CHAIN A AND CHAIN B"
        );
        assert!(state.remark3_continuation.is_some());
        feed_remark3(&mut state, 1, "FREE R VALUE: 0.45");

        assert_eq!(state.metadata.refinement.len(), 3);
        assert_eq!(state.metadata.refinement[0].id, "1");
        assert_eq!(state.metadata.refinement[1].id, "2");
        assert_eq!(state.metadata.refinement[2].id, "3");
        assert_eq!(state.metadata.refinement[0].basic.r_free, 0.0);
        assert!(state.metadata.refinement[1].basic.r_free.is_nan());
        assert_eq!(state.metadata.refinement[2].basic.r_free, 0.45);
        assert_eq!(state.source_state.resolution, 3.3);
        assert!(state.remark3_continuation.is_none());

        let mut prior = BioStructureSourceState::default();
        prior.resolution = 1.25;
        let mut prepopulated = PdbReaderState::new("p26-prior.pdb", false, prior);
        feed_remark3(&mut prepopulated, 1, "DATA USED IN REFINEMENT.");
        feed_remark3(
            &mut prepopulated,
            1,
            "RESOLUTION RANGE HIGH (ANGSTROMS): 2.5",
        );
        assert_eq!(prepopulated.source_state.resolution, 1.25);
    }

    #[test]
    fn gemmi_pdb_remark3_tls_groups_selections_and_tensors_match_pinned_source() {
        let mut state = PdbReaderState::new("p27.pdb", false, BioStructureSourceState::default());
        feed_remark3(&mut state, 0, "DATA USED IN REFINEMENT.");
        feed_remark3(&mut state, 0, "TLS GROUP: 42");
        feed_remark3(&mut state, 3, "SET: CHAIN A");
        feed_remark3(&mut state, 7, "AND CHAIN B");
        feed_remark3(&mut state, 3, "SELECTION: CHAIN B");
        feed_remark3(
            &mut state,
            4,
            "SELECTION          : NCS should not enter TLS",
        );
        feed_remark3(&mut state, 0, "ORIGIN FOR THE GROUP (A): 1.25 -2.5 3.75");
        feed_remark3(
            &mut state,
            0,
            "T11: 1.1 T12: 1.2 T13: 1.3 T21: 2.1 T22: 2.2 T23: 2.3 T31: 3.1 T32: 3.2 T33: 3.3",
        );
        feed_remark3(
            &mut state,
            0,
            "L11: 4.1 L12: 4.2 L13: 4.3 L21: 5.1 L22: 5.2 L23: 5.3 L31: 6.1 L32: 6.2 L33: 6.3",
        );
        feed_remark3(
            &mut state,
            0,
            "S11: 7.1 S12: 7.2 S13: 7.3 S21: 8.1 S22: 8.2 S23: 8.3 S31: 9.1 S32: 9.2 S33: 9.3",
        );
        feed_remark3(&mut state, 0, "RESIDUE RANGE:A        1     A       10 ");
        feed_remark3(&mut state, 0, "TLS GROUP: -7");
        feed_remark3(&mut state, 0, "TLS GROUP: +8");
        feed_remark3(&mut state, 0, "TLS GROUP: 32768");

        assert_eq!(state.metadata.refinement.len(), 1);
        let groups = &state.metadata.refinement[0].tls_groups;
        assert_eq!(groups.len(), 4);

        let group = &groups[0];
        assert_eq!(group.id, "42");
        assert_eq!(group.num_id, 42);
        assert_eq!(group.selections.len(), 3);
        assert_eq!(group.selections[0].details, "CHAIN A AND CHAIN B");
        assert_eq!(group.selections[1].details, "CHAIN B");
        assert_eq!(group.selections[2].chain.as_str(), "A");
        assert_eq!(group.selections[2].res_begin, PdbSeqId::new(1, None));
        assert_eq!(group.selections[2].res_end, PdbSeqId::new(10, None));
        assert_eq!(group.origin, [1.25, -2.5, 3.75]);
        assert_eq!(group.t, [1.1, 2.2, 3.3, 2.1, 3.1, 3.2]);
        assert_eq!(group.l, [4.1, 5.2, 6.3, 5.1, 6.1, 6.2]);
        assert_eq!(
            group.s,
            [[7.1, 7.2, 7.3], [8.1, 8.2, 8.3], [9.1, 9.2, 9.3],]
        );

        assert_eq!(groups[1].id, "-7");
        assert_eq!(groups[1].num_id, 0);
        assert_eq!(groups[2].id, "+8");
        assert_eq!(groups[2].num_id, 0);
        assert_eq!(groups[3].id, "32768");
        assert_eq!(groups[3].num_id, i16::MIN);
        for default_group in &groups[1..] {
            assert!(default_group.selections.is_empty());
            assert_eq!(default_group.origin, [0.0; 3]);
            assert!(default_group.t.iter().all(|component| component.is_nan()));
            assert!(default_group.l.iter().all(|component| component.is_nan()));
            assert!(
                default_group
                    .s
                    .iter()
                    .flatten()
                    .all(|component| component.is_nan())
            );
        }
    }

    #[test]
    fn gemmi_pdb_remark3_tls_gates_ranges_short_fields_and_latest_group_match_source() {
        let mut state = PdbReaderState::new(
            "p27-branches.pdb",
            false,
            BioStructureSourceState::default(),
        );
        feed_remark3(&mut state, 0, "DATA USED IN REFINEMENT.");

        // Source TLS assignments with no group do not synthesize one.
        feed_remark3(&mut state, 3, "SET: orphan selection");
        feed_remark3(&mut state, 3, "SELECTION: orphan selection");
        feed_remark3(&mut state, 0, "ORIGIN FOR THE GROUP (A): 1 2 3");
        feed_remark3(&mut state, 0, "T11: 4.0");
        assert!(state.metadata.refinement[0].tls_groups.is_empty());

        feed_remark3(&mut state, 0, "TLS GROUP: 10");
        // The source returns without assignment unless the origin has exactly
        // three tokens.
        feed_remark3(&mut state, 0, "ORIGIN FOR THE GROUP (A): 1 2");
        feed_remark3(&mut state, 3, "SET: CHAIN A");
        feed_remark3(&mut state, 7, "AND CHAIN B");
        feed_remark3(&mut state, 3, "SELECTION: CHAIN B");
        feed_remark3(
            &mut state,
            4,
            "SELECTION          : NCS selection is not TLS",
        );
        feed_remark3(&mut state, 0, "RESIDUE RANGE:A        1A    A       10 ");
        feed_remark3(&mut state, 0, "RESIDUE RANGE:A       BAD    A       11 ");
        feed_remark3(&mut state, 0, "RESIDUE RANGE:A        1     B       12 ");
        // The last fixed-width sequence field has only one source byte before
        // the zero-filled source line buffer supplies its terminator.
        feed_remark3(&mut state, 0, "RESIDUE RANGE:A        1     A       1");
        feed_remark3(&mut state, 0, "T11: 1.0 T12: 2.0 S33: 3.0 odd");
        feed_remark3(&mut state, 0, "T11: 4.0 T99: 9.0");

        let first = &state.metadata.refinement[0].tls_groups[0];
        assert_eq!(first.id, "10");
        assert_eq!(first.num_id, 10);
        assert_eq!(first.origin, [0.0; 3]);
        assert_eq!(first.selections.len(), 4);
        assert_eq!(first.selections[0].details, "CHAIN A AND CHAIN B");
        assert_eq!(first.selections[1].details, "CHAIN B");
        assert_eq!(first.selections[2].res_begin, PdbSeqId::new(1, Some(b'a')));
        assert_eq!(first.selections[2].res_end, PdbSeqId::new(10, None));
        assert_eq!(first.selections[3].res_begin, PdbSeqId::new(1, None));
        assert_eq!(first.selections[3].res_end, PdbSeqId::new(1, None));
        assert_eq!(first.t[0], 4.0);
        assert_eq!(first.t[3], 2.0);
        assert!(first.t[1].is_nan());
        assert!(first.t[2].is_nan());
        assert!(first.t[4].is_nan());
        assert!(first.t[5].is_nan());
        assert_eq!(first.s[2][2], 3.0);
        assert!(first.s[0][0].is_nan());

        feed_remark3(&mut state, 0, "TLS GROUP: 11");
        feed_remark3(&mut state, 0, "ORIGIN FOR THE GROUP (A): 4 5 6");
        assert_eq!(state.metadata.refinement[0].tls_groups[0].origin, [0.0; 3]);
        assert_eq!(state.metadata.refinement[0].tls_groups[1].id, "11");
        assert_eq!(
            state.metadata.refinement[0].tls_groups[1].origin,
            [4.0, 5.0, 6.0]
        );

        feed_remark3(&mut state, 0, "DATA USED IN REFINEMENT.");
        feed_remark3(&mut state, 0, "TLS GROUP: 99");
        feed_remark3(&mut state, 0, "ORIGIN FOR THE GROUP (A): 7 8 9");
        assert_eq!(state.metadata.refinement.len(), 2);
        assert_eq!(state.metadata.refinement[1].tls_groups[0].id, "99");
        assert_eq!(
            state.metadata.refinement[1].tls_groups[0].origin,
            [7.0, 8.0, 9.0]
        );
        assert_eq!(
            state.metadata.refinement[0].tls_groups[1].origin,
            [4.0, 5.0, 6.0]
        );
    }

    #[test]
    fn gemmi_pdb_remark3_tls_chain_width_boundary_is_typed_and_non_mutating() {
        let mut state = PdbReaderState::new(
            "p27-chain-width.pdb",
            false,
            BioStructureSourceState::default(),
        );
        feed_remark3(&mut state, 0, "DATA USED IN REFINEMENT.");
        feed_remark3(&mut state, 0, "TLS GROUP: 1");

        let row = format!(
            "REMARK   3 RESIDUE RANGE:{}{}{}{}{}",
            "ABCDE", "    1 ", "    ", "ABCDE", "   10 "
        );
        assert_eq!(
            state.remark3_record(&row),
            Err(PdbRemark3Error::TlsChainIdNotRepresentable { width: 5 })
        );
        assert!(
            state.metadata.refinement[0].tls_groups[0]
                .selections
                .is_empty()
        );
    }

    #[test]
    fn gemmi_pdb_remark200_experimental_defaults_scattering_types_and_ids_match_source() {
        let mut state = PdbReaderState::new(
            "p30-defaults.pdb",
            false,
            BioStructureSourceState::default(),
        );
        feed_remark200(&mut state, 200, "EXPERIMENTAL DETAILS");
        feed_remark200(&mut state, 230, "EXPERIMENTAL DETAILS");
        feed_remark200(&mut state, 240, "EXPERIMENTAL DETAILS");
        feed_remark200(&mut state, 240, "IN THE HIGHEST RESOLUTION SHELL.");

        assert_eq!(state.metadata.crystals.len(), 3);
        assert_eq!(state.metadata.experiments.len(), 3);
        let expected_scattering = ["x-ray", "neutron", "electron"];
        for (index, (crystal, experiment)) in state
            .metadata
            .crystals
            .iter()
            .zip(&state.metadata.experiments)
            .enumerate()
        {
            let id = (index + 1).to_string();
            assert_eq!(crystal.id, id);
            assert_eq!(crystal.diffractions.len(), 1);
            assert_eq!(crystal.description, "");
            let diffraction = &crystal.diffractions[0];
            assert_eq!(diffraction.id, id);
            assert_eq!(diffraction.scattering_type, expected_scattering[index]);
            assert!(diffraction.temperature.is_nan());
            assert_eq!(diffraction.source, "");
            assert_eq!(diffraction.source_type, "");
            assert_eq!(diffraction.synchrotron, "");
            assert_eq!(diffraction.beamline, "");
            assert_eq!(diffraction.wavelengths, "");
            assert_eq!(diffraction.mono_or_laue, 0);
            assert_eq!(diffraction.monochromator, "");
            assert_eq!(diffraction.collection_date, "");
            assert_eq!(diffraction.optics, "");
            assert_eq!(diffraction.detector, "");
            assert_eq!(diffraction.detector_make, "");
            assert_eq!(experiment.diffraction_ids, [id]);
            assert_eq!(experiment.method, "");
            assert_eq!(experiment.number_of_crystals, -1);
            assert_eq!(experiment.unique_reflections, -1);
            assert!(experiment.b_wilson.is_nan());
            assert!(experiment.reflections.resolution_high.is_nan());
            assert!(experiment.reflections.resolution_low.is_nan());
            assert!(experiment.reflections.completeness.is_nan());
            assert!(experiment.reflections.redundancy.is_nan());
            assert!(experiment.reflections.r_merge.is_nan());
            assert!(experiment.reflections.r_sym.is_nan());
            assert!(experiment.reflections.mean_i_over_sigma.is_nan());
            assert!(crystal.ph.is_nan());
            assert_eq!(crystal.ph_range, "");
            assert!(experiment.shells.is_empty() || index == 2);
        }
        assert_eq!(state.metadata.experiments[2].shells.len(), 1);
        assert!(
            state.metadata.experiments[2].shells[0]
                .resolution_high
                .is_nan()
        );
    }

    #[test]
    fn gemmi_pdb_remark200_fields_software_continuation_and_shells_match_pinned_oracle() {
        let mut state =
            PdbReaderState::new("p30-fields.pdb", false, BioStructureSourceState::default());
        // The source ignores experiment fields until EXPERIMENTAL DETAILS has
        // established the linked experiment/crystal/diffraction rows.
        feed_remark200(&mut state, 200, "EXPERIMENT TYPE: SHOULD BE IGNORED");
        feed_remark200(
            &mut state,
            200,
            "INTENSITY-INTEGRATION SOFTWARE: MOSFLM 1.0",
        );
        feed_remark200(&mut state, 200, "DATA SCALING SOFTWARE: AIMLESS 0.2");
        feed_remark200(&mut state, 200, "SOFTWARE USED: PHASER 2.8");
        feed_remark200(
            &mut state,
            200,
            "METHOD USED TO DETERMINE THE STRUCTURE: Molecular replacement",
        );
        feed_remark200(&mut state, 200, "STARTING MODEL: 1ABC");
        feed_remark200(&mut state, 200, "EXPERIMENTAL DETAILS");
        feed_remark200(&mut state, 200, "EXPERIMENT TYPE: X-RAY DIFFRACTION");
        feed_remark200(&mut state, 200, "NUMBER OF CRYSTALS USED: +2 trailing");
        feed_remark200(&mut state, 200, "PH: 6.5");
        feed_remark200(&mut state, 200, "DATE OF DATA COLLECTION: 28-MAR-07");
        feed_remark200(&mut state, 200, "TEMPERATURE           (KELVIN): 100.25");
        feed_remark200(&mut state, 200, "SYNCHROTRON              (Y/N): Y");
        feed_remark200(&mut state, 200, "RADIATION SOURCE: SSRL");
        feed_remark200(&mut state, 200, "BEAMLINE: BL 9-2");
        feed_remark200(&mut state, 200, "X-RAY GENERATOR MODEL: GENERATOR X");
        feed_remark200(&mut state, 200, "MONOCHROMATIC OR LAUE    (M/L): M");
        feed_remark200(&mut state, 200, "WAVELENGTH OR RANGE        (A): 0.9793");
        feed_remark200(&mut state, 200, "MONOCHROMATOR: Si(111)");
        feed_remark200(&mut state, 200, "OPTICS: mirrors");
        feed_remark200(&mut state, 200, "DETECTOR TYPE: PILATUS");
        feed_remark200(&mut state, 200, "DETECTOR MANUFACTURER: Dectris");
        feed_remark200(&mut state, 200, "NUMBER OF UNIQUE REFLECTIONS: 2345tail");
        feed_remark200(&mut state, 200, "RESOLUTION RANGE HIGH      (A): 1.25");
        feed_remark200(&mut state, 200, "RESOLUTION RANGE LOW       (A): 30.0");
        feed_remark200(&mut state, 200, "COMPLETENESS FOR RANGE     (%): 98.5");
        feed_remark200(&mut state, 200, "DATA REDUNDANCY: 3.2");
        feed_remark200(&mut state, 200, "R MERGE                    (I): 0.08");
        feed_remark200(&mut state, 200, "R SYM                      (I): 0.09");
        feed_remark200(&mut state, 200, "<I/SIGMA(I)> FOR THE DATA SET: 12.4");
        feed_remark200(&mut state, 200, "IN THE HIGHEST RESOLUTION SHELL.");
        feed_remark200(
            &mut state,
            200,
            "HIGHEST RESOLUTION SHELL, RANGE HIGH (A): 1.25",
        );
        feed_remark200(
            &mut state,
            200,
            "HIGHEST RESOLUTION SHELL, RANGE LOW  (A): 1.32",
        );
        feed_remark200(&mut state, 200, "COMPLETENESS FOR SHELL     (%): 90.1");
        feed_remark200(&mut state, 200, "DATA REDUNDANCY IN SHELL: 2.4");
        feed_remark200(&mut state, 200, "R MERGE FOR SHELL          (I): 0.3");
        feed_remark200(&mut state, 200, "R SYM FOR SHELL            (I): 0.4");
        feed_remark200(&mut state, 200, "<I/SIGMA(I)> FOR SHELL: 5.6");
        feed_remark200(&mut state, 200, "REMARK: first line");
        feed_remark200(&mut state, 200, " second line   ");
        feed_remark200(&mut state, 230, "STARTING MODEL: 2DEF");
        feed_remark200(&mut state, 200, "EXPERIMENT TYPE: NULL");

        assert_eq!(state.metadata.solved_by, "Molecular replacement");
        assert_eq!(state.metadata.starting_model, "2DEF");
        let software = &state.metadata.software;
        assert_eq!(software.len(), 3);
        assert_eq!(
            software
                .iter()
                .map(|item| (
                    item.name.as_str(),
                    item.version.as_str(),
                    item.classification
                ))
                .collect::<Vec<_>>(),
            [
                ("MOSFLM", "1.0", BioSoftwareClassification::DataReduction),
                ("AIMLESS", "0.2", BioSoftwareClassification::DataScaling),
                ("PHASER", "2.8", BioSoftwareClassification::Phasing),
            ]
        );

        let crystal = &state.metadata.crystals[0];
        assert_eq!(crystal.id, "1");
        assert_eq!(crystal.description, "first line second line");
        assert_eq!(crystal.ph, 6.5);
        assert_eq!(crystal.ph_range, "");
        let diffraction = &crystal.diffractions[0];
        assert_eq!(diffraction.id, "1");
        assert_eq!(diffraction.scattering_type, "x-ray");
        assert_eq!(diffraction.collection_date, "2007-03-28");
        assert_eq!(diffraction.temperature, 100.25);
        assert_eq!(diffraction.source, "SYNCHROTRON");
        assert_eq!(diffraction.synchrotron, "SSRL");
        assert_eq!(diffraction.beamline, "BL 9-2");
        assert_eq!(diffraction.source_type, "GENERATOR X");
        assert_eq!(diffraction.mono_or_laue, b'M');
        assert_eq!(diffraction.wavelengths, "0.9793");
        assert_eq!(diffraction.monochromator, "Si(111)");
        assert_eq!(diffraction.optics, "mirrors");
        assert_eq!(diffraction.detector, "PILATUS");
        assert_eq!(diffraction.detector_make, "Dectris");

        let experiment = &state.metadata.experiments[0];
        assert_eq!(experiment.method, "X-RAY DIFFRACTION");
        assert_eq!(experiment.number_of_crystals, 2);
        assert_eq!(experiment.unique_reflections, 2345);
        assert_eq!(experiment.diffraction_ids, ["1"]);
        let reflections = &experiment.reflections;
        assert_eq!(reflections.resolution_high, 1.25);
        assert_eq!(reflections.resolution_low, 30.0);
        assert_eq!(reflections.completeness, 98.5);
        assert_eq!(reflections.redundancy, 3.2);
        assert_eq!(reflections.r_merge, 0.08);
        assert_eq!(reflections.r_sym, 0.09);
        assert_eq!(reflections.mean_i_over_sigma, 12.4);
        assert_eq!(experiment.shells.len(), 1);
        let shell = &experiment.shells[0];
        assert_eq!(shell.resolution_high, 1.25);
        assert_eq!(shell.resolution_low, 1.32);
        assert_eq!(shell.completeness, 90.1);
        assert_eq!(shell.redundancy, 2.4);
        assert_eq!(shell.r_merge, 0.3);
        assert_eq!(shell.r_sym, 0.4);
        assert_eq!(shell.mean_i_over_sigma, 5.6);
    }

    #[test]
    fn gemmi_pdb_remark200_pH_range_numeric_prefix_and_undefined_screening_are_explicit() {
        let mut state =
            PdbReaderState::new("p30-ph.pdb", false, BioStructureSourceState::default());
        feed_remark200(&mut state, 200, "EXPERIMENTAL DETAILS");
        feed_remark200(&mut state, 200, "PH: not a number");
        assert!(state.metadata.crystals[0].ph.is_nan());
        assert_eq!(state.metadata.crystals[0].ph_range, "not a number");

        // Gemmi's `is_space` and `fast_from_chars` both accept C-locale
        // whitespace; `skip_blank` does not consume vertical tab, so it must
        // remain in the numeric parser's input.
        feed_remark200(&mut state, 200, "PH:\u{000b}+6.5");
        assert_eq!(state.metadata.crystals[0].ph, 6.5);
        assert_eq!(state.metadata.crystals[0].ph_range, "not a number");

        feed_remark200(&mut state, 200, "NUMBER OF CRYSTALS USED: no digits");
        feed_remark200(&mut state, 200, "NUMBER OF UNIQUE REFLECTIONS:");
        assert_eq!(state.metadata.experiments[0].number_of_crystals, 0);
        assert_eq!(state.metadata.experiments[0].unique_reflections, 0);

        // The source's pre-increment plus loop-body increment reads beyond
        // NUL for this two-digit fractional spelling. Rust keeps that source-
        // undefined case outside parity and returns the documented boundary.
        assert_eq!(
            state.remark_200_230_240_record("REMARK 200 PH: 7.00"),
            Err(super::PdbRemark200Error::UndefinedPhNumericScreening)
        );
        assert_eq!(state.metadata.crystals[0].ph, 6.5);
    }

    #[test]
    fn gemmi_pdb_remark2_keeps_first_resolution_in_source_order() {
        let resolution_remark =
            |value: &str| format!("REMARK   2 RESOLUTION.    {value} ANGSTROMS.");

        let mut state = PdbReaderState::new(
            "p31-resolution.pdb",
            false,
            BioStructureSourceState::default(),
        );
        state
            .remark_2_300_record(&resolution_remark("9.99").replace("ANGSTROMS", "angstroms"))
            .unwrap();
        state
            .remark_2_300_record("REMARK   2 RESOLUTION.    8.88 OTHER")
            .unwrap();
        assert_eq!(state.source_state.resolution, 0.0);

        state
            .remark_2_300_record(&resolution_remark("2.50"))
            .unwrap();
        state
            .remark_2_300_record(&resolution_remark("1.25"))
            .unwrap();
        assert_eq!(state.source_state.resolution, 2.5);

        // Gemmi applies the REMARK 3 high-resolution fallback in raw source
        // order. Once it has supplied a value, a later REMARK 2 cannot replace
        // it because the source condition requires resolution == 0.0.
        let mut remark3_first = PdbReaderState::new(
            "p31-remark3-first.pdb",
            false,
            BioStructureSourceState::default(),
        );
        feed_remark3(&mut remark3_first, 1, "DATA USED IN REFINEMENT.");
        feed_remark3(
            &mut remark3_first,
            1,
            "RESOLUTION RANGE HIGH (ANGSTROMS): 1.75",
        );
        remark3_first
            .remark_2_300_record(&resolution_remark("2.50"))
            .unwrap();
        assert_eq!(remark3_first.source_state.resolution, 1.75);

        let mut prior = BioStructureSourceState::default();
        prior.resolution = 1.25;
        let mut prepopulated = PdbReaderState::new("p31-prior.pdb", false, prior);
        prepopulated
            .remark_2_300_record(&resolution_remark("2.50"))
            .unwrap();
        assert_eq!(prepopulated.source_state.resolution, 1.25);

        let mut nul_terminated =
            PdbReaderState::new("p31-nul.pdb", false, BioStructureSourceState::default());
        nul_terminated
            .remark_2_300_record("REMARK   2 RESOLUTION.    2.50 \0ANGSTROMS.")
            .unwrap();
        assert_eq!(nul_terminated.source_state.resolution, 0.0);
    }

    #[test]
    fn gemmi_pdb_remark300_initial_and_continuation_trimming_match_source() {
        let mut state =
            PdbReaderState::new("p31-detail.pdb", false, BioStructureSourceState::default());

        state
            .remark_2_300_record("REMARK 300 ignored before detail")
            .unwrap();
        state
            .remark_2_300_record("REMARK 300 remark: wrong case")
            .unwrap();
        assert!(state.metadata.remark_300_detail.is_empty());

        state
            .remark_2_300_record("REMARK 300 REMARK:  Biological details  \t")
            .unwrap();
        assert_eq!(state.metadata.remark_300_detail, "Biological details");

        state
            .remark_2_300_record("REMARK 300   continuation detail \t  ")
            .unwrap();
        assert_eq!(
            state.metadata.remark_300_detail,
            "Biological details\n  continuation detail"
        );

        // Once detail is nonempty the source appends each sufficiently long
        // REMARK 300 payload, even if it itself begins with the initial marker.
        state.remark_2_300_record("REMARK 300   \t").unwrap();
        state
            .remark_2_300_record("REMARK 300 REMARK: later marker")
            .unwrap();
        assert_eq!(
            state.metadata.remark_300_detail,
            "Biological details\n  continuation detail\n\nREMARK: later marker"
        );

        // The outer source loop skips all rows of 11 bytes or fewer, including
        // a blank REMARK 300 row after detail has already started.
        let before = state.metadata.remark_300_detail.clone();
        state.remark_2_300_record("REMARK 300 ").unwrap();
        assert_eq!(state.metadata.remark_300_detail, before);
    }

    #[test]
    fn gemmi_pdb_remark2_short_triggered_field_is_source_undefined_and_atomic() {
        let mut state =
            PdbReaderState::new("p31-short.pdb", false, BioStructureSourceState::default());
        assert_eq!(
            state.remark_2_300_record("REMARK   2 ANGSTROM"),
            Err(PdbRemarkMetadataError::UndefinedResolutionField {
                remark_length: "REMARK   2 ANGSTROM".len()
            })
        );
        assert_eq!(state.source_state.resolution, 0.0);
        assert!(state.metadata.remark_300_detail.is_empty());

        // Source short-circuit order avoids touching the fixed field when a
        // prior resolution is already present.
        state.source_state.resolution = 1.25;
        state.remark_2_300_record("REMARK   2 ANGSTROM").unwrap();
        assert_eq!(state.source_state.resolution, 1.25);
    }

    #[test]
    fn gemmi_pdb_remark350_assemblies_metadata_chains_and_operators_match_pinned_oracle() {
        let mut state = PdbReaderState::new(
            "p32-assemblies.pdb",
            false,
            BioStructureSourceState::default(),
        );

        assert_eq!(
            state.remark_350_record(&remark350_biomolecule("1")),
            Ok(true)
        );
        for row in [
            remark350_biomt(b'1', "  1", [1.0, 0.0, 0.0, 10.0]),
            remark350_biomt(b'2', "  1", [0.0, 1.0, 0.0, 20.0]),
            remark350_biomt(b'3', "  1", [0.0, 0.0, 1.0, 30.0]),
        ] {
            assert_eq!(state.remark_350_record(&row), Ok(false));
        }
        assert_eq!(
            state.remark_350_record(&remark350_apply("A,  B,,LONG", true)),
            Ok(false)
        );
        assert_eq!(
            state.remark_350_record(&remark350_apply("B,UNRESOLVED A", false)),
            Ok(false)
        );
        assert_eq!(
            state.remark_350_record(&remark350_biomt(b'3', "  1", [0.0, 0.0, 1.0, 30.0])),
            Ok(false)
        );
        assert_eq!(
            state.remark_350_record(&remark350_biomt(b'3', "  2", [9.0, 8.0, 7.0, 70.0])),
            Ok(false)
        );

        assert_eq!(
            state.remark_350_record(&remark350_biomolecule("2")),
            Ok(true)
        );
        assert_eq!(
            state.remark_350_record(&remark350_metadata_row(
                "AUTHOR DETERMINED",
                44,
                45,
                "dimeric assembly",
            )),
            Ok(false)
        );
        assert_eq!(
            state.remark_350_record(&remark350_metadata_row(
                "SOFTWARE DETERMINED",
                51,
                52,
                "software oligomer",
            )),
            Ok(false)
        );
        assert_eq!(
            state.remark_350_record(&remark350_metadata_row(
                "SOFTWARE USED",
                24,
                25,
                "Gemmi 0.7.5"
            )),
            Ok(false)
        );
        assert_eq!(
            state.remark_350_record(&remark350_metadata_row(
                "TOTAL BURIED SURFACE AREA",
                36,
                37,
                "123.5",
            )),
            Ok(false)
        );
        assert_eq!(
            state.remark_350_record(&remark350_metadata_row(
                "SURFACE AREA OF THE COMPLEX",
                38,
                39,
                "456.75",
            )),
            Ok(false)
        );
        assert_eq!(
            state.remark_350_record(&remark350_metadata_row(
                "CHANGE IN SOLVENT FREE ENERGY",
                40,
                41,
                "-7.25",
            )),
            Ok(false)
        );

        assert_eq!(state.assemblies.len(), 2);
        let first = &state.assemblies[0];
        assert_eq!(first.name, "1");
        assert!(!first.author_determined);
        assert!(!first.software_determined);
        assert_eq!(first.special_kind, BioAssemblySpecialKind::NotApplicable);
        assert_eq!(first.oligomeric_count, 0);
        assert!(first.oligomeric_details.is_empty());
        assert!(first.software_name.is_empty());
        assert!(first.buried_surface_area.is_nan());
        assert!(first.surface_area.is_nan());
        assert!(first.solvent_free_energy_change.is_nan());
        assert_eq!(first.generators.len(), 1);
        assert_eq!(
            first.generators[0].chains,
            ["A", "B", "LONG", "B", "UNRESOLVED", "A"]
        );
        let operators = &first.generators[0].operators;
        assert_eq!(operators.len(), 2);
        assert_eq!(operators[0].name.as_deref(), Some("1"));
        assert_eq!(operators[0].operator_type, None);
        assert_eq!(
            *operators[0].transform.matrix(),
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        );
        assert_eq!(*operators[0].transform.translation(), [10.0, 20.0, 30.0]);
        assert_eq!(operators[1].name.as_deref(), Some("2"));
        assert_eq!(
            *operators[1].transform.matrix(),
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [9.0, 8.0, 7.0]]
        );
        assert_eq!(*operators[1].transform.translation(), [0.0, 0.0, 70.0]);

        let second = &state.assemblies[1];
        assert_eq!(second.name, "2");
        assert!(second.author_determined);
        assert!(second.software_determined);
        assert_eq!(second.oligomeric_details, "software oligomer");
        assert_eq!(second.software_name, "Gemmi 0.7.5");
        assert_eq!(second.buried_surface_area, 123.5);
        assert_eq!(second.surface_area, 456.75);
        assert_eq!(second.solvent_free_energy_change, -7.25);
        assert!(second.generators.is_empty());
    }

    #[test]
    fn gemmi_pdb_remark350_continue_and_partial_matrix_state_match_source() {
        let mut state =
            PdbReaderState::new("p32-partial.pdb", false, BioStructureSourceState::default());

        // Before a BIOMOLECULE row, Gemmi's outer-loop continue ignores a
        // valid-looking APPLY row and retains no synthetic assembly.
        assert_eq!(
            state.remark_350_record(&remark350_apply("MISSING", true)),
            Ok(true)
        );
        assert!(state.assemblies.is_empty());

        assert_eq!(
            state.remark_350_record(&remark350_biomolecule("3")),
            Ok(true)
        );
        // Continuation without a preceding generator is another explicit
        // source `continue`; no generator is inferred from this row.
        assert_eq!(
            state.remark_350_record(&remark350_apply("ORPHAN", false)),
            Ok(true)
        );
        // Rows 1 and 2 arrive before a generator. Row 3 therefore must not
        // emit/reset an operator, but its transform row remains in scratch.
        assert_eq!(
            state.remark_350_record(&remark350_biomt(b'1', "  1", [1.0, 0.0, 0.0, 10.0])),
            Ok(false)
        );
        assert_eq!(
            state.remark_350_record(&remark350_biomt(b'2', "  1", [0.0, 1.0, 0.0, 20.0])),
            Ok(false)
        );
        assert_eq!(
            state.remark_350_record(&remark350_biomt(b'3', "  1", [0.0, 0.0, 1.0, 30.0])),
            Ok(false)
        );
        assert!(state.assemblies[0].generators.is_empty());
        assert_eq!(state.remark350_matrix.matrix()[2], [0.0, 0.0, 1.0]);
        assert_eq!(state.remark350_matrix.translation()[2], 30.0);

        assert_eq!(
            state.remark_350_record(&remark350_apply("KNOWN,UNKNOWN", true)),
            Ok(false)
        );
        assert_eq!(
            state.remark_350_record(&remark350_biomt(b'3', "  1", [0.0, 0.0, 1.0, 30.0])),
            Ok(false)
        );
        assert_eq!(state.assemblies[0].generators.len(), 1);
        assert_eq!(
            state.assemblies[0].generators[0].chains,
            ["KNOWN", "UNKNOWN"]
        );
        assert_eq!(state.assemblies[0].generators[0].operators.len(), 1);
        assert_eq!(
            *state.assemblies[0].generators[0].operators[0]
                .transform
                .matrix(),
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        );
        assert_eq!(
            *state.assemblies[0].generators[0].operators[0]
                .transform
                .translation(),
            [10.0, 20.0, 30.0]
        );
        assert_eq!(
            *state.remark350_matrix.matrix(),
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        );
        assert_eq!(*state.remark350_matrix.translation(), [0.0, 0.0, 0.0]);
    }

    #[test]
    fn gemmi_pdb_remark350_source_undefined_matrix_field_is_typed_and_atomic() {
        let mut state = PdbReaderState::new(
            "p32-short-matrix.pdb",
            false,
            BioStructureSourceState::default(),
        );
        state
            .remark_350_record(&remark350_biomolecule("4"))
            .unwrap();
        state
            .remark_350_record(&remark350_apply("A", true))
            .unwrap();
        let before_assemblies = state.assemblies.clone();
        let before_matrix = state.remark350_matrix;
        let mut short_row = remark350_biomt(b'3', "  1", [0.0, 0.0, 1.0, 30.0]);
        short_row.truncate(65);

        assert_eq!(
            state.remark_350_record(&short_row),
            Err(PdbRemark350Error::SourceUndefinedAccess {
                offset: 58,
                width: 10,
                record_length: 65,
            })
        );
        assert_eq!(state.assemblies.len(), before_assemblies.len());
        for (actual, before) in state.assemblies.iter().zip(&before_assemblies) {
            assert_eq!(actual.name, before.name);
            assert_eq!(actual.author_determined, before.author_determined);
            assert_eq!(actual.software_determined, before.software_determined);
            assert_eq!(actual.special_kind, before.special_kind);
            assert_eq!(actual.oligomeric_count, before.oligomeric_count);
            assert_eq!(actual.oligomeric_details, before.oligomeric_details);
            assert_eq!(actual.software_name, before.software_name);
            assert_eq!(
                actual.buried_surface_area.to_bits(),
                before.buried_surface_area.to_bits()
            );
            assert_eq!(actual.surface_area.to_bits(), before.surface_area.to_bits());
            assert_eq!(
                actual.solvent_free_energy_change.to_bits(),
                before.solvent_free_energy_change.to_bits()
            );
            assert_eq!(actual.generators, before.generators);
        }
        assert_eq!(state.remark350_matrix, before_matrix);
    }

    fn append_dispatch_line(text: &mut String, bytes: &[u8]) {
        text.push_str(std::str::from_utf8(bytes).expect("fixed PDB fixture is ASCII"));
        if !bytes.ends_with(b"\n") {
            text.push('\n');
        }
    }

    fn dispatch_test_reader<'text>(
        text: &'text str,
        options: PdbInputOptions,
    ) -> PdbBioReader<'text, 'static> {
        PdbBioReader::new(text, "dispatch.pdb", options, false, false)
    }

    fn dispatch_hierarchy_lines(reader: &mut PdbBioReader<'_, '_>) -> Result<(), PdbReaderError> {
        loop {
            let Some(copied_len) = reader.cursor.copy_line().map(|line| line.len()) else {
                break;
            };
            reader
                .state
                .record_line_checked(&reader.cursor.line_buffer[..copied_len])
                .map_err(|cause| reader.stream_error(cause))?;
            reader.dispatch_hierarchy_record(copied_len)?;
            if reader.transition.stopped_at_end_record() {
                break;
            }
        }
        Ok(())
    }

    fn dispatch_metadata_lines(reader: &mut PdbBioReader<'_, '_>) -> Result<(), PdbReaderError> {
        loop {
            let Some(copied_len) = reader.cursor.copy_line().map(|line| line.len()) else {
                break;
            };
            reader
                .state
                .record_line_checked(&reader.cursor.line_buffer[..copied_len])
                .map_err(|cause| reader.stream_error(cause))?;

            if reader.dispatch_hierarchy_record(copied_len)? {
                if reader.transition.stopped_at_end_record() {
                    break;
                }
                continue;
            }

            if reader.dispatch_metadata_record(copied_len)? == Some(PdbDbRefAction::Stop) {
                break;
            }
        }
        Ok(())
    }

    fn dispatch_crystal_lines(
        reader: &mut PdbBioReader<'_, '_>,
    ) -> Result<Vec<Vec<u8>>, PdbReaderError> {
        let mut unhandled_records = Vec::new();
        loop {
            let Some(copied_len) = reader.cursor.copy_line().map(|line| line.len()) else {
                break;
            };
            reader
                .state
                .record_line_checked(&reader.cursor.line_buffer[..copied_len])
                .map_err(|cause| reader.stream_error(cause))?;

            if reader.dispatch_hierarchy_record(copied_len)? {
                if reader.transition.stopped_at_end_record() {
                    break;
                }
                continue;
            }

            if reader.dispatch_metadata_record(copied_len)? == Some(PdbDbRefAction::Stop) {
                break;
            }
            if reader.dispatch_crystal_format_record(copied_len)? {
                continue;
            }

            unhandled_records.push(reader.cursor.line_buffer[..copied_len].to_vec());
        }
        Ok(unhandled_records)
    }

    fn dispatch_metadata_row(record: [u8; 6], payload: &[u8]) -> Vec<u8> {
        let mut row = Vec::with_capacity(11 + payload.len());
        row.extend_from_slice(&record);
        row.extend_from_slice(b"    ");
        row.extend_from_slice(payload);
        row.push(b'\n');
        row
    }

    fn append_dispatch_metadata_row(text: &mut String, record: [u8; 6], payload: &[u8]) {
        let row = dispatch_metadata_row(record, payload);
        append_dispatch_line(text, &row);
    }

    #[test]
    fn pdb_metadata_dispatch_orphan_dbref2_stops_before_later_atom() {
        // Gemmi creates/promotes the DBREF2 entity, then breaks the source
        // record loop if no DBREF1/DBREF row exists. The later ATOM must not
        // be dispatched, while finalization remains the caller's next stage.
        let mut text = String::new();
        let mut orphan = [b' '; 80];
        orphan[..5].copy_from_slice(b"DBREF");
        orphan[5] = b'2';
        orphan[11..13].copy_from_slice(b"C ");
        append_dispatch_line(&mut text, &orphan);
        append_dispatch_atom(&mut text, *b"ATOM  ", b'C', 1);

        let mut reader = dispatch_test_reader(&text, PdbInputOptions::default());
        dispatch_metadata_lines(&mut reader).unwrap();

        assert_eq!(reader.state.line_number, 1);
        assert!(reader.cursor.position < reader.cursor.input.len());
        assert_eq!(reader.state.entities.entities.len(), 1);
        assert_eq!(reader.state.entities.entities[0].source_name, b"C");
        assert_eq!(
            reader.state.entities.entities[0].entity_kind,
            EntityKind::Polymer
        );
        assert!(reader.state.entities.entities[0].dbrefs.is_empty());
        assert!(reader.transition.grouping.models.is_empty());
    }

    #[test]
    fn pdb_metadata_dispatch_stages_connections_in_order_and_routes_conect_seqres() {
        use std::collections::BTreeMap;

        let mut text = String::new();
        let connection_rows = [
            dispatch_metadata_row(*b"SSBOND", b"first"),
            dispatch_metadata_row(*b"LINK  ", b"second"),
            dispatch_metadata_row(*b"CISPEP", b"third"),
        ];
        for row in &connection_rows {
            append_dispatch_line(&mut text, row);
        }
        let conect = conect_source_line(9, &[20, 7], 21);
        append_dispatch_line(&mut text, &conect[..21]);
        let sequence = seqres_source_line(*b"SEQRES", *b"A ", &[(19, *b"ALA"), (27, *b"GLY")], 30);
        append_dispatch_line(&mut text, &sequence);

        let mut reader = dispatch_test_reader(&text, PdbInputOptions::default());
        dispatch_metadata_lines(&mut reader).unwrap();

        assert_eq!(reader.connection_records, connection_rows);
        assert_eq!(
            reader.state.source_state.conect_map,
            BTreeMap::from([(9, vec![20, 7])])
        );
        assert_eq!(reader.state.entities.entities.len(), 1);
        let entity = &reader.state.entities.entities[0];
        assert_eq!(entity.source_name, b"A");
        assert_eq!(entity.entity_kind, EntityKind::Polymer);
        assert_eq!(entity.full_sequence, [b"ALA".to_vec(), b"GLY".to_vec()]);
    }

    #[test]
    fn pdb_metadata_dispatch_author_continuation_and_raw_remarks_ignore_skip_gate() {
        let mut text = String::new();
        append_dispatch_metadata_row(&mut text, *b"AUTHOR", b"Doe, A.-");
        append_dispatch_metadata_row(&mut text, *b"AUTHOR", b"B.Smith, C.Jones");
        append_dispatch_metadata_row(&mut text, *b"REMARK", b"  2 first\r");
        append_dispatch_metadata_row(&mut text, *b"REMARK", b"  3 second");

        let expected_remarks = [
            "REMARK      2 first".to_owned(),
            "REMARK      3 second".to_owned(),
        ];

        for skip_remarks in [false, true] {
            let mut reader = PdbBioReader::new(
                &text,
                "metadata.pdb",
                PdbInputOptions::default(),
                false,
                skip_remarks,
            );
            dispatch_metadata_lines(&mut reader).unwrap();
            reader.state.finalize_author_names();

            assert_eq!(
                reader.state.metadata.authors,
                [
                    "Doe".to_owned(),
                    "Smith, A.-B.".to_owned(),
                    "Jones, C.".to_owned()
                ]
            );
            assert_eq!(reader.state.source_state.raw_remarks, expected_remarks);
        }
    }

    #[test]
    fn pdb_metadata_dispatch_short_helix_and_sheet_rows_are_consumed_without_values() {
        let helix = helix_source_line(
            *b"HELIX ", *b" A", *b"  12B", *b"ALA", *b" B", *b"  34C", *b"GLY", *b" 5", *b"   17",
            38,
        );
        let sheet = sheet_source_line(*b"SHEET ", 38, &[]);
        let mut text = String::new();
        append_dispatch_line(&mut text, &helix[..38]);
        append_dispatch_line(&mut text, &sheet[..38]);

        let mut reader = dispatch_test_reader(&text, PdbInputOptions::default());
        dispatch_metadata_lines(&mut reader).unwrap();

        assert_eq!(reader.state.line_number, 2);
        assert!(reader.state.helices.is_empty());
        assert!(reader.state.sheets.is_empty());
    }

    #[test]
    fn pdb_crystal_dispatch_interleaves_scale_and_mtrix_staging_and_keeps_identity_ncs() {
        let mut text = String::new();

        for row in [
            pdb_mtrix_line(b'1', "I1", ["1", "0", "0", "0"], 80, None),
            pdb_mtrix_line(b'2', "I2", ["0", "1", "0", "0"], 80, None),
            pdb_mtrix_line(b'3', "I3", ["0", "0", "1", "0"], 80, Some(b'1')),
        ] {
            append_dispatch_line(&mut text, &row);
        }

        append_dispatch_line(&mut text, &pdb_scale_line(b'1', ["0.5", "0", "0", "0.25"]));
        append_dispatch_line(&mut text, b"NOPE ignored row");
        append_dispatch_line(
            &mut text,
            &pdb_mtrix_line(b'2', "M2", ["0", "0.25", "0", "0.5"], 80, None),
        );
        append_dispatch_line(
            &mut text,
            &pdb_mtrix_line(b'3', "N7", ["0", "0", "0.2", "0.75"], 80, Some(b'1')),
        );

        append_dispatch_line(
            &mut text,
            &pdb_mtrix_line(b'1', "M1", ["0.75", "0", "0", "1"], 80, None),
        );
        append_dispatch_line(&mut text, &pdb_scale_line(b'2', ["0", "3", "0", "2"]));
        append_dispatch_line(&mut text, &pdb_scale_line(b'3', ["0", "0", "4", "3"]));

        let mut reader = dispatch_test_reader(&text, PdbInputOptions::default());
        let unhandled = dispatch_crystal_lines(&mut reader).unwrap();

        assert_eq!(unhandled, [b"NOPE ignored row\n".to_vec()]);
        assert_eq!(reader.state.line_number, 10);
        assert_eq!(reader.state.ncs_operators.len(), 1);
        assert_eq!(reader.state.ncs_operators[0].id, "N7");
        assert!(reader.state.ncs_operators[0].given);
        assert_eq!(
            reader.state.ncs_operators[0].transform,
            BioTransform::new(
                [[0.5, 0.0, 0.0], [0.0, 0.25, 0.0], [0.0, 0.0, 0.2]],
                [0.25, 0.5, 0.75]
            )
        );
        assert_eq!(
            reader
                .state
                .source_state
                .info
                .get("_struct_ncs_oper.id")
                .map(String::as_str),
            Some("I3")
        );
        assert_eq!(
            *reader.state.crystal_info.fractional(),
            BioTransform::new(
                [[0.75, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0]],
                [1.0, 2.0, 3.0]
            )
        );
        assert!(reader.state.crystal_info.explicit_matrices());
        assert_eq!(reader.state.matrix, BioTransform::identity());
    }

    #[test]
    fn pdb_crystal_dispatch_preserves_origx_cryst1_scale_order_and_unknown_noop() {
        let mut text = String::new();
        append_dispatch_line(
            &mut text,
            &pdb_cryst1_line(10.0, 20.0, 30.0, 90.0, 90.0, 90.0, "P 1", "5"),
        );
        append_dispatch_line(
            &mut text,
            &pdb_origx_line(b'2', ["5.25", "-6.5", "7.75", "-8.125"]),
        );
        append_dispatch_line(
            &mut text,
            &pdb_origx_line(b'1', ["1.25", "2.5", "-3.75", "4.125"]),
        );
        append_dispatch_line(&mut text, b"UNKN must not alter matrix staging");
        for row in [
            pdb_scale_line(b'1', ["0.5", "0", "0", "0.25"]),
            pdb_scale_line(b'2', ["0", "0.25", "0", "-0.5"]),
            pdb_scale_line(b'3', ["0", "0", "0.2", "0.75"]),
        ] {
            append_dispatch_line(&mut text, &row);
        }
        append_dispatch_line(
            &mut text,
            &pdb_cryst1_line(11.0, 21.0, 31.0, 90.0, 90.0, 90.0, "P 2", "9"),
        );
        append_dispatch_line(
            &mut text,
            &pdb_origx_line(b'3', ["9.25", "10.5", "11.75", "12.125"]),
        );

        let mut reader = dispatch_test_reader(&text, PdbInputOptions::default());
        let unhandled = dispatch_crystal_lines(&mut reader).unwrap();

        assert_eq!(
            unhandled,
            [b"UNKN must not alter matrix staging\n".to_vec()]
        );
        assert_eq!(reader.state.line_number, 9);
        assert!(reader.state.source_state.has_origx);
        assert_eq!(
            reader.state.source_state.origx,
            BioTransform::new(
                [[1.25, 2.5, -3.75], [5.25, -6.5, 7.75], [9.25, 10.5, 11.75]],
                [4.125, -8.125, 12.125]
            )
        );
        assert_eq!(
            reader.state.crystal_info.cell(),
            BioCrystalCell {
                a: 11.0,
                b: 21.0,
                c: 31.0,
                alpha: 90.0,
                beta: 90.0,
                gamma: 90.0,
            }
        );
        assert_eq!(reader.state.crystal_info.space_group_hm(), Some("P 2"));
        assert_eq!(reader.state.crystal_info.z_pdb(), Some("9"));
        assert_eq!(
            *reader.state.crystal_info.fractional(),
            BioTransform::new(
                [[0.5, 0.0, 0.0], [0.0, 0.25, 0.0], [0.0, 0.0, 0.2]],
                [0.25, -0.5, 0.75]
            )
        );
        assert!(reader.state.crystal_info.explicit_matrices());
    }

    #[test]
    fn pdb_crystal_dispatch_checks_format_before_and_after_active_model() {
        let mut before_model = dispatch_test_reader("DaTa_model\n", PdbInputOptions::default());
        let error = dispatch_crystal_lines(&mut before_model).unwrap_err();
        assert_eq!(error.line_number, Some(1));
        assert_eq!(error.stage, PdbReaderStage::Record(*b"DaTa"));
        assert!(matches!(
            error.cause,
            PdbReaderCause::InputFormat(PdbInputFormatError::Cif { source })
                if source == "dispatch.pdb"
        ));

        let mut active_then_inactive = String::new();
        let model = pdb_model_source_line(1);
        append_dispatch_line(&mut active_then_inactive, &model[..15]);
        append_dispatch_line(&mut active_then_inactive, b"data_model");
        append_dispatch_line(&mut active_then_inactive, b"ENDMDL");
        append_dispatch_line(&mut active_then_inactive, b"data_model");
        let mut reader = dispatch_test_reader(&active_then_inactive, PdbInputOptions::default());
        let error = dispatch_crystal_lines(&mut reader).unwrap_err();
        assert_eq!(reader.state.line_number, 4);
        assert_eq!(error.line_number, Some(4));
        assert_eq!(error.stage, PdbReaderStage::Record(*b"data"));
        assert!(matches!(
            error.cause,
            PdbReaderCause::InputFormat(PdbInputFormatError::Cif { source })
                if source == "dispatch.pdb"
        ));

        let mut mmjson = dispatch_test_reader("{\"DaTa_\n", PdbInputOptions::default());
        let error = dispatch_crystal_lines(&mut mmjson).unwrap_err();
        assert_eq!(error.line_number, Some(1));
        assert!(matches!(
            error.cause,
            PdbReaderCause::InputFormat(PdbInputFormatError::Mmjson { source })
                if source == "dispatch.pdb"
        ));
    }

    fn append_dispatch_atom(text: &mut String, record: [u8; 6], chain: u8, sequence: i32) {
        let line = grouped_pdb_atom_line(
            record, sequence, *b" CA ", *b"GLY", chain, sequence, b' ', *b"    ",
        );
        append_dispatch_line(text, &line);
    }

    #[test]
    fn pdb_dispatch_implicit_explicit_models_and_endmdl_are_distinct() {
        // Pinned Gemmi populate_structure_from_pdb_stream: an atom without
        // MODEL creates model-count+1; ENDMDL clears the active model; a later
        // atom then creates the next implicit model. This drives the actual
        // reader dispatcher over its fixed-buffer cursor, not helper calls.
        let mut implicit_text = String::new();
        append_dispatch_atom(&mut implicit_text, *b"ATOM  ", b'A', 1);
        let mut implicit = dispatch_test_reader(&implicit_text, PdbInputOptions::default());
        dispatch_hierarchy_lines(&mut implicit).unwrap();
        assert_eq!(implicit.transition.grouping.models.len(), 1);
        assert_eq!(
            implicit.transition.grouping.models[0].source_model_number,
            Some(1)
        );

        let mut text = String::new();
        let model = pdb_model_source_line(5);
        append_dispatch_line(&mut text, &model[..15]);
        append_dispatch_atom(&mut text, *b"ATOM  ", b'A', 1);
        append_dispatch_line(&mut text, b"ENDMDL\n");
        append_dispatch_atom(&mut text, *b"HETATM", b'B', 2);
        let mut reader = dispatch_test_reader(&text, PdbInputOptions::default());
        dispatch_hierarchy_lines(&mut reader).unwrap();

        assert_eq!(reader.transition.grouping.models.len(), 2);
        assert_eq!(
            reader
                .transition
                .grouping
                .models
                .iter()
                .map(|model| model.source_model_number)
                .collect::<Vec<_>>(),
            [Some(5), Some(2)]
        );
        assert_eq!(
            reader.transition.grouping.models[0].chains[0]
                .source_id
                .as_str(),
            "A"
        );
        assert_eq!(
            reader.transition.grouping.models[1].chains[0]
                .source_id
                .as_str(),
            "B"
        );
    }

    #[test]
    fn pdb_dispatch_anisou_order_and_duplicate_errors_keep_record_context() {
        // Gemmi associates ANISOU to the last atom and treats nonzero u11 as
        // the duplicate sentinel. Both source errors retain their stream line
        // and the four-byte dispatched record tag.
        let first_anisou = pdb_anisou_source_line(1, *b" CA ", [1, 2, 3, 4, 5, 6]);
        let mut no_atom_text = String::new();
        append_dispatch_line(&mut no_atom_text, &first_anisou[..81]);
        let mut no_atom = dispatch_test_reader(&no_atom_text, PdbInputOptions::default());
        let error = dispatch_hierarchy_lines(&mut no_atom).unwrap_err();
        assert_eq!(error.line_number, Some(1));
        assert_eq!(error.stage, PdbReaderStage::Record(*b"ANIS"));
        assert!(matches!(
            error.cause,
            PdbReaderCause::ModelTransition(PdbModelTransitionError::AnisouWithoutAtom)
        ));

        let atom = grouped_pdb_atom_line(*b"ATOM  ", 1, *b" CA ", *b"GLY", b'A', 1, b' ', *b"    ");
        let first_anisou =
            pdb_anisou_source_line(999, *b" XX ", [10000, 20000, 30000, 40000, 50000, 60000]);
        let duplicate_anisou =
            pdb_anisou_source_line(1, *b" CA ", [12345, 23456, 34567, 45678, 56789, 67890]);
        let mut text = String::new();
        append_dispatch_line(&mut text, &atom);
        append_dispatch_line(&mut text, &first_anisou[..81]);
        append_dispatch_line(&mut text, &duplicate_anisou[..81]);
        let mut reader = dispatch_test_reader(&text, PdbInputOptions::default());
        let error = dispatch_hierarchy_lines(&mut reader).unwrap_err();
        assert_eq!(error.line_number, Some(3));
        assert_eq!(error.stage, PdbReaderStage::Record(*b"ANIS"));
        assert!(matches!(
            error.cause,
            PdbReaderCause::ModelTransition(PdbModelTransitionError::DuplicateAnisou)
        ));
        assert_ne!(
            reader.transition.grouping.models[0].chains[0].residues[0].atoms[0].anisou[0],
            0.0
        );
    }

    #[test]
    fn pdb_dispatch_ter_ignore_split_and_error_rollback_follow_source() {
        let mut ignored_text = String::new();
        append_dispatch_atom(&mut ignored_text, *b"ATOM  ", b'A', 1);
        append_dispatch_line(&mut ignored_text, b"TER \n");
        append_dispatch_atom(&mut ignored_text, *b"ATOM  ", b'A', 2);
        let mut ignored = dispatch_test_reader(
            &ignored_text,
            PdbInputOptions {
                ignore_ter: true,
                ..PdbInputOptions::default()
            },
        );
        dispatch_hierarchy_lines(&mut ignored).unwrap();
        assert_eq!(ignored.state.source_state.ter_status, 0);
        assert_eq!(ignored.transition.grouping.models[0].chains.len(), 1);
        assert_eq!(
            ignored.transition.grouping.models[0].chains[0]
                .residues
                .iter()
                .map(|residue| residue.entity_kind)
                .collect::<Vec<_>>(),
            [EntityKind::Unknown, EntityKind::Unknown]
        );

        let mut split_text = String::new();
        append_dispatch_atom(&mut split_text, *b"ATOM  ", b'A', 1);
        append_dispatch_line(&mut split_text, b"TER \n");
        append_dispatch_atom(&mut split_text, *b"ATOM  ", b'A', 2);
        let mut split = dispatch_test_reader(
            &split_text,
            PdbInputOptions {
                split_chain_on_ter: true,
                ..PdbInputOptions::default()
            },
        );
        dispatch_hierarchy_lines(&mut split).unwrap();
        assert_eq!(split.state.source_state.ter_status, b'y');
        let split_chains = &split.transition.grouping.models[0].chains;
        assert_eq!(split_chains.len(), 2);
        assert!(split_chains.iter().all(|chain| {
            chain.residues.len() == 1 && chain.residues[0].entity_kind == EntityKind::Unknown
        }));

        let mut repeated_ter_text = String::new();
        append_dispatch_atom(&mut repeated_ter_text, *b"ATOM  ", b'A', 1);
        append_dispatch_line(&mut repeated_ter_text, b"TER \n");
        append_dispatch_line(&mut repeated_ter_text, b"TER \n");
        let mut repeated_ter = dispatch_test_reader(&repeated_ter_text, PdbInputOptions::default());
        dispatch_hierarchy_lines(&mut repeated_ter).unwrap();
        assert_eq!(repeated_ter.state.source_state.ter_status, b'e');
        assert_eq!(
            repeated_ter.transition.grouping.models[0].chains[0].residues[0].entity_kind,
            EntityKind::Polymer
        );
        repeated_ter
            .transition
            .finalize_ter_entity_types(&repeated_ter.state.source_state);
        assert_eq!(
            repeated_ter.transition.grouping.models[0].chains[0].residues[0].entity_kind,
            EntityKind::Unknown
        );
    }

    #[test]
    fn pdb_dispatch_endmdl_continues_but_end_stops_before_later_records() {
        let mut text = String::new();
        let model = pdb_model_source_line(5);
        append_dispatch_line(&mut text, &model[..15]);
        append_dispatch_atom(&mut text, *b"ATOM  ", b'A', 1);
        append_dispatch_line(&mut text, b"ENDMDL\n");
        append_dispatch_atom(&mut text, *b"ATOM  ", b'B', 2);
        append_dispatch_line(&mut text, b"END\n");
        append_dispatch_atom(&mut text, *b"ATOM  ", b'C', 3);
        let mut reader = dispatch_test_reader(&text, PdbInputOptions::default());
        dispatch_hierarchy_lines(&mut reader).unwrap();

        assert!(reader.transition.stopped_at_end_record());
        assert_eq!(reader.state.line_number, 5);
        assert_eq!(reader.transition.grouping.models.len(), 2);
        assert_eq!(
            reader
                .transition
                .grouping
                .models
                .iter()
                .map(|model| model.source_model_number)
                .collect::<Vec<_>>(),
            [Some(5), Some(2)]
        );
        assert!(
            reader.transition.grouping.models[0]
                .chains
                .iter()
                .all(|chain| chain.source_id.as_str() != "C")
        );

        let mut invalid_text = String::new();
        let first_model = pdb_model_source_line(1);
        append_dispatch_line(&mut invalid_text, &first_model[..15]);
        append_dispatch_atom(&mut invalid_text, *b"ATOM  ", b'A', 1);
        let second_model = pdb_model_source_line(2);
        append_dispatch_line(&mut invalid_text, &second_model[..15]);
        let mut invalid = dispatch_test_reader(&invalid_text, PdbInputOptions::default());
        let error = dispatch_hierarchy_lines(&mut invalid).unwrap_err();
        assert_eq!(error.line_number, Some(3));
        assert_eq!(error.stage, PdbReaderStage::Record(*b"MODE"));
        assert!(matches!(
            error.cause,
            PdbReaderCause::ModelTransition(PdbModelTransitionError::ModelWithoutEnd)
        ));
    }

    #[test]
    fn pdb_complete_reader_finalizes_empty_header_and_unterminated_eof_inputs() {
        let empty = PdbBioReader::new(
            "",
            "/fixtures/empty.pdb.gz",
            PdbInputOptions::default(),
            false,
            false,
        )
        .read_structure()
        .unwrap();
        assert_eq!(empty.input_format(), BioCoordinateFormat::Pdb);
        assert_eq!(empty.source_state().name, "empty");
        assert_eq!(empty.models().len(), 1);
        assert_eq!(empty.atoms().len(), 0);

        let header = PdbBioReader::new(
            "HEADER\n",
            "header.pdb",
            PdbInputOptions::default(),
            false,
            false,
        )
        .read_structure()
        .unwrap();
        assert_eq!(header.models().len(), 1);
        assert!(header.models()[0].chain_span().is_empty());
        assert!(header.atoms().is_empty());

        let final_atom = pdb_atom_field_line();
        let final_line = std::str::from_utf8(&final_atom).unwrap();
        let unterminated = PdbBioReader::new(
            final_line,
            "unterminated.pdb",
            PdbInputOptions::default(),
            false,
            false,
        )
        .read_structure()
        .unwrap();
        assert_eq!(unterminated.models().len(), 1);
        assert_eq!(unterminated.atoms().len(), 1);
    }

    #[test]
    fn pdb_complete_reader_end_and_orphan_dbref2_breaks_still_finalize() {
        let mut end_input = String::new();
        append_dispatch_atom(&mut end_input, *b"ATOM  ", b'A', 1);
        append_dispatch_line(&mut end_input, b"END\n");
        append_dispatch_atom(&mut end_input, *b"ATOM  ", b'B', 2);
        let ended = PdbBioReader::new(
            &end_input,
            "end.pdb",
            PdbInputOptions::default(),
            false,
            false,
        )
        .read_structure()
        .unwrap();
        assert_eq!(ended.models().len(), 1);
        assert_eq!(ended.atoms().len(), 1);
        assert_eq!(ended.residues().len(), 1);

        let mut dbref2_input = String::new();
        let mut orphan = [b' '; 80];
        orphan[..5].copy_from_slice(b"DBREF");
        orphan[5] = b'2';
        orphan[11..13].copy_from_slice(b"C ");
        append_dispatch_line(&mut dbref2_input, &orphan);
        append_dispatch_atom(&mut dbref2_input, *b"ATOM  ", b'C', 1);
        let after_orphan = PdbBioReader::new(
            &dbref2_input,
            "orphan.pdb",
            PdbInputOptions::default(),
            false,
            false,
        )
        .read_structure()
        .unwrap();
        assert_eq!(after_orphan.models().len(), 1);
        assert!(after_orphan.atoms().is_empty());
        assert_eq!(after_orphan.entities().len(), 1);
        assert_eq!(after_orphan.entities()[0].kind(), EntityKind::Polymer);
        assert!(after_orphan.entities()[0].dbrefs().is_empty());
    }

    #[test]
    fn pdb_complete_reader_applies_bounded_line_and_non_ascii_stream_options() {
        let mut input = "\u{00e9}UNKNOWN\n".to_owned();
        let first =
            grouped_pdb_atom_line(*b"ATOM  ", 1, *b" CA ", *b"GLY", b'A', 1, b' ', *b"    ");
        input.push_str(std::str::from_utf8(&first).unwrap());
        input.push_str("discarded overlong suffix\n");
        let final_atom =
            grouped_pdb_atom_line(*b"ATOM  ", 2, *b" CA ", *b"GLY", b'A', 2, b' ', *b"    ");
        input.push_str(std::str::from_utf8(&final_atom).unwrap());

        let parsed = PdbBioReader::new(
            &input,
            "/fixtures/stream.pdb.gz",
            PdbInputOptions {
                max_line_length: 80,
                ..PdbInputOptions::default()
            },
            true,
            false,
        )
        .read_structure()
        .unwrap();
        assert_eq!(parsed.source_state().name, "stream");
        assert_eq!(parsed.source_state().non_ascii_line, 1);
        assert_eq!(parsed.atoms().len(), 2);
        assert_eq!(parsed.residues().len(), 2);
    }

    #[test]
    fn pdb_complete_reader_uses_reusable_buffer_after_embedded_nul() {
        let complete_discriminator = PdbBioReader::new(
            "{\"data_\0\nNEXT\n",
            "embedded.json",
            PdbInputOptions::default(),
            false,
            false,
        )
        .read_structure()
        .unwrap_err();
        assert_eq!(complete_discriminator.line_number, Some(1));
        assert_eq!(
            complete_discriminator.stage,
            PdbReaderStage::Record(*b"{\"da")
        );
        assert!(matches!(
            complete_discriminator.cause,
            PdbReaderCause::InputFormat(PdbInputFormatError::Mmjson { source })
                if source == "embedded.json"
        ));

        let incomplete_discriminator = PdbBioReader::new(
            "{\"data\0_\nNEXT\n",
            "not-mmjson.pdb",
            PdbInputOptions::default(),
            false,
            false,
        )
        .read_structure()
        .unwrap();
        assert_eq!(incomplete_discriminator.models().len(), 1);
        assert!(incomplete_discriminator.atoms().is_empty());
    }

    #[test]
    fn pdb_complete_reader_finalizes_ssbond_after_stream_dispatch() {
        let mut input = String::new();
        let mut ssbond = [b' '; 80];
        ssbond[..6].copy_from_slice(b"SSBOND");
        ssbond[7..10].copy_from_slice(b"  1");
        ssbond[11..14].copy_from_slice(b"CYS");
        ssbond[14..16].copy_from_slice(b"A ");
        ssbond[17..22].copy_from_slice(b"   1 ");
        ssbond[25..28].copy_from_slice(b"CYS");
        ssbond[28..30].copy_from_slice(b"A ");
        ssbond[31..36].copy_from_slice(b"   2 ");
        append_dispatch_line(&mut input, &ssbond);

        for sequence in [1, 2] {
            let mut atom = grouped_pdb_atom_line(
                *b"ATOM  ", sequence, *b" SG ", *b"CYS", b'A', sequence, b' ', *b"    ",
            );
            atom[16] = b' ';
            atom[76..78].copy_from_slice(b" S");
            append_dispatch_line(&mut input, &atom);
        }

        let parsed = PdbBioReader::new(
            &input,
            "disulfide.pdb",
            PdbInputOptions::default(),
            false,
            false,
        )
        .read_structure()
        .unwrap();
        assert_eq!(parsed.atoms().len(), 2);
        assert_eq!(parsed.connections().len(), 1);
        let connection = &parsed.connections()[0];
        assert_eq!(connection.name, "disulf1");
        assert_eq!(connection.kind, BioConnectionKind::Disulf);
        assert_eq!(connection.partner1.logical_atom_name(), "SG");
        assert_eq!(connection.partner2.logical_atom_name(), "SG");
    }

    #[test]
    fn pdb_complete_reader_finalizes_ccd_entity_and_remark_options() {
        let mut input = String::new();
        let sequence = seqres_source_line(*b"SEQRES", *b"A ", &[(19, *b"AAA")], 23);
        append_dispatch_line(&mut input, &sequence);
        let alias = hetnam_source_line(*b"HETNAM", *b"AAA", b"BBB", 80, b' ');
        append_dispatch_line(&mut input, &alias[..80]);
        let atom = grouped_pdb_atom_line(*b"HETATM", 1, *b" CA ", *b"AAA", b'A', 1, b' ', *b"    ");
        append_dispatch_line(&mut input, &atom);
        append_dispatch_line(&mut input, b"TER   \n");

        let parsed = PdbBioReader::new(
            &input,
            "aliases.pdb",
            PdbInputOptions::default(),
            false,
            false,
        )
        .read_structure()
        .unwrap();
        assert_eq!(
            parsed.residues()[0].name(),
            ResidueName::from_ascii(b"BBB").unwrap()
        );
        assert_eq!(parsed.entities().len(), 1);
        assert_eq!(parsed.entities()[0].kind(), EntityKind::Polymer);
        assert_eq!(parsed.entities()[0].full_sequence(), &["BBB".to_owned()]);
        assert!(!parsed.entities()[0].subchains().is_empty());

        let remark = "REMARK   2 RESOLUTION.    2.50 ANGSTROMS.\n";
        let interpreted = PdbBioReader::new(
            remark,
            "remarks.pdb",
            PdbInputOptions::default(),
            false,
            false,
        )
        .read_structure()
        .unwrap();
        let skipped = PdbBioReader::new(
            remark,
            "remarks.pdb",
            PdbInputOptions::default(),
            false,
            true,
        )
        .read_structure()
        .unwrap();
        assert_eq!(interpreted.source_state().raw_remarks.len(), 1);
        assert_eq!(skipped.source_state().raw_remarks.len(), 1);
        assert_eq!(interpreted.source_state().resolution, 2.5);
        assert_eq!(skipped.source_state().resolution, 0.0);
    }
}
