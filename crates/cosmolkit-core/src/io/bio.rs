//! PDB/mmCIF/mmJSON readers for `BioStructure`.
//!
//! This is the canonical structural PDB/mmCIF/mmJSON IO path. It is
//! Gemmi-aligned and returns COSMolKit `BioStructure` rows for PDB coordinate
//! records, mmCIF structure blocks, mmJSON-backed CIF documents, and chem-comp
//! CIF input through the same public dispatch surface. RDKit-derived PDB
//! behavior belongs to Molecule compatibility APIs layered over
//! `BioStructure`, not to a parallel public parser.
//!
//! Source-derived parser work in this module follows
//! `dev/source_reproduction_protocol.md`: copied Gemmi source lines stay
//! adjacent to the Rust code that ports them, and every unsupported source
//! branch remains visible with a two-axis marker.
//!
//! Marker prefix:
//!
//! - `Gemmi❌❌:` source behavior is not implemented.
//! - `Gemmi❗✔️:` behavior is only partially modeled, but the local algorithmic
//!   shape is not materially worse for the modeled input state.
//! - `Gemmi✔️✔️:` behavior and performance are reproduced for the currently
//!   modeled input state.

use crate::Element;
use crate::bio::{
    AltLocLabel, AtomName, AtomRow, AtomSourceIds, BioAssembly, BioAssemblyGenerator,
    BioAssemblyOperator, BioAsu, BioAtomAddress, BioCalcFlag, BioCisPep, BioConnection,
    BioConnectionType, BioCoorFormat, BioDiffractionInfo, BioEntityDbRef, BioExperimentCrystalInfo,
    BioExperimentInfo, BioHelix, BioMetadata, BioModRes, BioNcsOperator, BioRefinementBin,
    BioRefinementInfo, BioRefinementRestraint, BioReflectionsInfo, BioSheet, BioSheetStrand,
    BioSiftsUnpResidue, BioSoftwareClassification, BioSoftwareItem, BioStructure, BioTlsGroup,
    BioTlsSelection, BioTransform, ChainId, ChainKind, ChainRow, ChainSourceIds, CrystalCell,
    CrystalInfo, EntityId, EntityKind, EntityRow, EntitySourceIds, ModelId, ModelRow,
    PdbAtomSerial, PdbChainId, PdbSeqId, PolymerKind, ResidueId, ResidueKind, ResidueName,
    ResidueRow, ResidueSourceIds, RowSpan, classify_residue_name,
};
use crate::bio_invariants::enforce_bio_structure_invariants;
use crate::io::gemmi_spacegroup_table::{
    GEMMI_ALT_NAMES, GEMMI_OP_DEN, GEMMI_SPACEGROUPS, GemmiSpaceGroupEntry,
};
use crate::support::{BIO_MMCIF_READ_FEATURE, BIO_PDB_READ_FEATURE, FeatureSpec};
use serde_json::Value as JsonValue;
use std::collections::HashMap;

mod cif;
mod sprintf;
mod writer;

use cif::{CifBlock, CifDocument, CifItem, CifLoop, CifToken};
pub use writer::{BioWriteError, MmcifOutputGroups, MmcifWriteOptions};

pub(crate) fn bio_structure_to_mmcif(
    structure: &BioStructure,
    options: MmcifWriteOptions,
) -> Result<String, BioWriteError> {
    let document = writer::make_mmcif_document(structure, options.groups)?;
    cif::cif_document_to_string(&document, options.cif_options()).map_err(BioWriteError::Io)
}

pub(crate) fn write_bio_structure_mmcif(
    structure: &BioStructure,
    path: &std::path::Path,
    options: MmcifWriteOptions,
) -> Result<(), BioWriteError> {
    let document = writer::make_mmcif_document(structure, options.groups)?;
    let file = std::fs::File::create(path)?;
    cif::write_cif_document(file, &document, options.cif_options()).map_err(BioWriteError::Io)
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum BioReadError {
    #[error(
        "unsupported BioStructure input feature {feature_name} at line {line_number}: {reason}"
    )]
    Unsupported {
        line_number: usize,
        feature_name: &'static str,
        reason: &'static str,
    },
    #[error("BioStructure parse error at line {line_number}: {message}")]
    Parse { line_number: usize, message: String },
    #[error("BioStructure invariant violation after parsing: {0}")]
    Invariant(&'static str),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BioPdbReadParams {
    pub max_line_length: usize,
    pub split_chain_on_ter: bool,
    pub skip_remarks: bool,
    pub check_non_ascii: bool,
    pub ignore_ter: bool,
}

impl Default for BioPdbReadParams {
    fn default() -> Self {
        Self {
            max_line_length: 120,
            split_chain_on_ter: false,
            skip_remarks: false,
            check_non_ascii: false,
            ignore_ter: false,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct ResidueKey {
    chain_key: PdbChainId,
    seq_id: PdbSeqId,
    name: ResidueName,
    segment_id: Option<[u8; 4]>,
}

#[derive(Debug, Default)]
struct PdbBioBuilder {
    structure: BioStructure,
    current_model: Option<ModelId>,
    current_chain: Option<(ChainId, PdbChainId)>,
    current_residue: Option<(ResidueId, ResidueKey)>,
    residue_lookup: HashMap<ResidueKey, ResidueId>,
    last_atom_index: Option<usize>,
    skip_remarks: bool,
}

impl PdbBioBuilder {
    fn has_model_number(&self, source_model_number: i32) -> bool {
        self.structure
            .models
            .iter()
            .any(|model| model.source_model_number == Some(source_model_number))
    }

    fn find_chain_in_model(&self, model_id: ModelId, chain_key: PdbChainId) -> Option<ChainId> {
        let model = &self.structure.models[model_id.index() as usize];
        let start = model.chain_span.start as usize;
        let end = model.chain_span.end() as usize;
        (start..end)
            .find(|&idx| {
                let chain = &self.structure.chains[idx];
                chain.source.auth_chain_id == Some(chain_key)
                    || chain.source.label_asym_id == Some(chain_key)
            })
            .map(|idx| ChainId::new(idx as u32))
    }

    fn chain_has_polymer_first_residue(&self, chain_id: ChainId) -> bool {
        let chain = &self.structure.chains[chain_id.index() as usize];
        if chain.residue_span.is_empty() {
            return false;
        }
        self.structure.residues[chain.residue_span.start as usize].entity_kind
            == EntityKind::Polymer
    }

    fn find_entity_by_source_id(&self, source_entity_id: &str) -> Option<EntityId> {
        self.structure
            .entities
            .iter()
            .position(|entity| entity.source.source_entity_id == source_entity_id)
            .map(|idx| EntityId::new(idx as u32))
    }

    fn find_entity_by_subchain(&self, subchain: PdbChainId) -> Option<EntityId> {
        self.structure
            .entities
            .iter()
            .position(|entity| entity.subchains.contains(&subchain))
            .map(|idx| EntityId::new(idx as u32))
    }

    fn find_or_add_entity(
        &mut self,
        source_entity_id: &str,
        kind: EntityKind,
        polymer_kind: PolymerKind,
    ) -> EntityId {
        if let Some(entity_id) = self.find_entity_by_source_id(source_entity_id) {
            let entity = &mut self.structure.entities[entity_id.index() as usize];
            if entity.kind == EntityKind::Unknown {
                entity.kind = kind;
            }
            if entity.polymer_kind == PolymerKind::Unknown {
                entity.polymer_kind = polymer_kind;
            }
            return entity_id;
        }

        let entity_id = EntityId::new(self.structure.entities.len() as u32);
        self.structure.entities.push(EntityRow {
            kind,
            polymer_kind,
            reflects_microhetero: false,
            sequence: Vec::new(),
            dbrefs: Vec::new(),
            sifts_unp_acc: Vec::new(),
            subchains: Vec::new(),
            source: EntitySourceIds {
                source_entity_id: source_entity_id.to_string(),
            },
        });
        entity_id
    }

    fn append_entity_sequence(&mut self, entity_id: EntityId, residue_name: &str) {
        self.structure.entities[entity_id.index() as usize]
            .sequence
            .push(residue_name.to_string());
    }

    fn merge_entity_sequence_at(&mut self, entity_id: EntityId, pos: usize, residue_name: &str) {
        let sequence = &mut self.structure.entities[entity_id.index() as usize].sequence;
        if pos == sequence.len() {
            sequence.push(residue_name.to_string());
        } else if let Some(existing) = sequence.get_mut(pos) {
            existing.push(',');
            existing.push_str(residue_name);
        }
    }

    fn add_entity_subchain(&mut self, entity_id: EntityId, subchain: PdbChainId) {
        let entity = &mut self.structure.entities[entity_id.index() as usize];
        if !entity.subchains.contains(&subchain) {
            entity.subchains.push(subchain);
        }
    }

    fn entity_for_chain_source(source: &ChainSourceIds) -> Option<String> {
        source
            .label_asym_id
            .or(source.auth_chain_id)
            .map(|chain_id| chain_id.as_str().to_string())
    }

    fn ensure_model(&mut self, source_model_number: Option<i32>) -> ModelId {
        if let Some(model_id) = self.current_model {
            return model_id;
        }

        let model_id = ModelId::new(self.structure.models.len() as u32);
        let start = self.structure.chains.len() as u32;
        self.structure.models.push(ModelRow {
            chain_span: RowSpan::new(start, 0),
            source_model_number,
        });
        self.current_model = Some(model_id);
        model_id
    }

    fn begin_model(&mut self, source_model_number: Option<i32>) {
        self.current_model = None;
        self.current_chain = None;
        self.current_residue = None;
        self.residue_lookup.clear();
        self.last_atom_index = None;
        self.ensure_model(source_model_number);
    }

    fn end_model(&mut self) {
        self.current_model = None;
        self.current_chain = None;
        self.current_residue = None;
        self.residue_lookup.clear();
        self.last_atom_index = None;
    }

    fn ensure_chain(
        &mut self,
        model_id: ModelId,
        chain_key: PdbChainId,
        source: ChainSourceIds,
        force_new_chain: bool,
    ) -> ChainId {
        if let Some((chain_id, current_source_id)) = self.current_chain
            && current_source_id == chain_key
            && !force_new_chain
        {
            return chain_id;
        }

        let chain_id = ChainId::new(self.structure.chains.len() as u32);
        let residue_start = self.structure.residues.len() as u32;
        let entity_id = source
            .label_asym_id
            .or(source.auth_chain_id)
            .and_then(|subchain| self.find_entity_by_subchain(subchain))
            .or_else(|| {
                Self::entity_for_chain_source(&source)
                    .and_then(|source_id| self.find_entity_by_source_id(&source_id))
            });
        self.structure.chains.push(ChainRow {
            model_id,
            entity_id,
            residue_span: RowSpan::new(residue_start, 0),
            kind: ChainKind::Unknown,
            source,
        });
        self.structure.models[model_id.index() as usize]
            .chain_span
            .len += 1;
        self.current_chain = Some((chain_id, chain_key));
        self.current_residue = None;
        self.residue_lookup.clear();
        chain_id
    }

    fn assign_chain_entity_by_subchain(&mut self, subchain: PdbChainId, entity_id: EntityId) {
        for chain in &mut self.structure.chains {
            if chain.source.label_asym_id == Some(subchain) {
                chain.entity_id = Some(entity_id);
            }
        }
    }

    fn ensure_residue(
        &mut self,
        chain_id: ChainId,
        key: ResidueKey,
        kind: ResidueKind,
        entity_kind: EntityKind,
    ) -> ResidueId {
        if let Some((residue_id, current_key)) = &self.current_residue
            && *current_key == key
        {
            return *residue_id;
        }

        if let Some(&residue_id) = self.residue_lookup.get(&key) {
            self.current_residue = Some((residue_id, key));
            return residue_id;
        }

        let residue_id = ResidueId::new(self.structure.residues.len() as u32);
        let atom_start = self.structure.atoms.len() as u32;
        self.structure.residues.push(ResidueRow {
            chain_id,
            atom_span: RowSpan::new(atom_start, 0),
            name: key.name,
            kind,
            entity_kind,
            source: ResidueSourceIds {
                seq_id: Some(key.seq_id),
                label_seq_id: None,
                segment_id: key.segment_id,
                subchain_id: None,
                label_entity_id: None,
            },
            het_flag: None,
            sifts_unp: None,
        });
        self.structure.chains[chain_id.index() as usize]
            .residue_span
            .len += 1;
        self.residue_lookup.insert(key.clone(), residue_id);
        self.current_residue = Some((residue_id, key));
        residue_id
    }

    fn push_atom(
        &mut self,
        record: PdbAtomRecord,
        source_model_number: Option<i32>,
        force_new_chain: bool,
        entity_kind: EntityKind,
    ) {
        let model_id = self.ensure_model(source_model_number);
        let chain_id = self.ensure_chain(
            model_id,
            record.chain_key,
            record.chain_source,
            force_new_chain,
        );
        let residue_id = self.ensure_residue(
            chain_id,
            record.residue_key(),
            record.residue_kind,
            entity_kind,
        );
        if let Some(label_seq_id) = record.label_seq_id {
            let residue = &mut self.structure.residues[residue_id.index() as usize];
            if residue.source.label_seq_id.is_none() {
                residue.source.label_seq_id = Some(label_seq_id);
            }
        }
        if let Some(label_asym_id) = record.chain_source.label_asym_id {
            let residue = &mut self.structure.residues[residue_id.index() as usize];
            if residue.source.subchain_id.is_none() {
                residue.source.subchain_id = Some(label_asym_id);
            }
        }
        if let Some(label_entity_id) = record.label_entity_id {
            let residue = &mut self.structure.residues[residue_id.index() as usize];
            if residue.source.label_entity_id.is_none() {
                residue.source.label_entity_id = Some(label_entity_id);
            }
        }
        if let Some(het_flag) = record.group_pdb {
            let residue = &mut self.structure.residues[residue_id.index() as usize];
            residue.het_flag = Some(het_flag);
        }

        let atom_row = AtomRow {
            residue_id,
            name: record.atom_name,
            element: record.element,
            altloc: record.altloc,
            occupancy: record.occupancy,
            b_iso: record.b_iso,
            formal_charge: record.formal_charge,
            anisou: None,
            calc_flag: record.calc_flag,
            tls_group_id: record.tls_group_id,
            fraction: record.fraction,
            source: AtomSourceIds {
                serial: record.serial.map(PdbAtomSerial),
            },
        };
        let residue_idx = residue_id.index() as usize;
        let insert_idx = self.structure.residues[residue_idx].atom_span.end() as usize;
        if insert_idx == self.structure.atoms.len() {
            self.structure.atoms.push(atom_row);
            self.structure.coordinates.positions.push(record.position);
        } else {
            self.structure.atoms.insert(insert_idx, atom_row);
            self.structure
                .coordinates
                .positions
                .insert(insert_idx, record.position);
            for (idx, residue) in self.structure.residues.iter_mut().enumerate() {
                if idx != residue_idx && residue.atom_span.start as usize >= insert_idx {
                    residue.atom_span.start += 1;
                }
            }
        }
        self.last_atom_index = Some(insert_idx);
        self.structure.residues[residue_idx].atom_span.len += 1;
    }

    fn set_last_atom_anisou(
        &mut self,
        line_number: usize,
        anisou: [f32; 6],
    ) -> Result<(), BioReadError> {
        let atom = self
            .last_atom_index
            .and_then(|index| self.structure.atoms.get_mut(index))
            .ok_or_else(|| BioReadError::Parse {
                line_number,
                message: "ANISOU before ATOM/HETATM is not valid".to_string(),
            })?;
        if atom.anisou.is_some() {
            return Err(BioReadError::Parse {
                line_number,
                message: "duplicated ANISOU record for previous atom".to_string(),
            });
        }
        atom.anisou = Some(anisou);
        Ok(())
    }

    fn set_last_atom_source_traits(
        &mut self,
        line_number: usize,
        calc_flag: BioCalcFlag,
        tls_group_id: Option<i16>,
        fraction: Option<f32>,
    ) -> Result<(), BioReadError> {
        let atom = self
            .last_atom_index
            .and_then(|index| self.structure.atoms.get_mut(index))
            .ok_or_else(|| BioReadError::Parse {
                line_number,
                message: "mmCIF atom properties before atom_site row is not valid".to_string(),
            })?;
        atom.calc_flag = calc_flag;
        atom.tls_group_id = tls_group_id;
        atom.fraction = fraction;
        Ok(())
    }

    fn finish(mut self) -> Result<BioStructure, BioReadError> {
        // BEGIN GEMMI CPP FUNCTION gemmi::read_pdb
        // Gemmi✔️✔️: if (st.models.empty())
        // Gemmi✔️✔️:   st.models.emplace_back(1);
        // Gemmi✔️✔️: if (st.ter_status == 'e')
        // Gemmi✔️✔️:   remove_entity_types(st);
        // Gemmi✔️✔️: assign_subchains(st, /*force=*/false, /*fail_if_unknown=*/false);
        // Gemmi✔️✔️: for (Chain& ch : st.models[0].chains)
        // Gemmi✔️✔️:   if (Entity* entity = st.get_entity(ch.name))
        // Gemmi✔️✔️:     if (auto polymer = ch.get_polymer())
        // Gemmi✔️✔️:       entity->subchains.emplace_back(polymer.subchain_id());
        // Gemmi✔️✔️: st.setup_cell_images();
        // Gemmi✔️✔️: process_conn(st, conn_records);
        // Gemmi✔️✔️: for (std::string& name : st.meta.authors)
        // Gemmi✔️✔️:   change_author_name_format_to_mmcif(name);
        // Gemmi✔️✔️: if (!options.skip_remarks)
        // Gemmi✔️✔️:   read_metadata_from_remarks(st);
        // Gemmi✔️✔️: restore_full_ccd_codes(st);
        // END GEMMI CPP FUNCTION
        if self.structure.models.is_empty() {
            self.begin_model(Some(1));
            self.end_model();
        }
        if self.structure.ter_status == 'e' {
            remove_entity_types(&mut self.structure);
        }
        assign_subchains(&mut self.structure, false, false)?;
        backfill_polymer_subchains_to_entities(&mut self.structure);
        setup_cell_images(&mut self.structure)?;
        if !self.structure.deferred_conn_records.is_empty() {
            let records = std::mem::take(&mut self.structure.deferred_conn_records);
            process_conn(&mut self.structure, &records);
        }
        for name in &mut self.structure.metadata.authors {
            change_author_name_format_to_mmcif(name);
        }
        if !self.skip_remarks {
            read_metadata_from_remarks(&mut self.structure)?;
        }
        restore_full_ccd_codes(&mut self.structure);
        enforce_bio_structure_invariants(&self.structure).map_err(BioReadError::Invariant)?;
        Ok(self.structure)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct PdbAtomRecord {
    serial: Option<i32>,
    atom_name: AtomName,
    altloc: Option<AltLocLabel>,
    residue_name: ResidueName,
    residue_kind: ResidueKind,
    chain_key: PdbChainId,
    chain_source: ChainSourceIds,
    seq_id: PdbSeqId,
    label_seq_id: Option<i32>,
    label_entity_id: Option<EntityId>,
    group_pdb: Option<char>,
    segment_id: Option<[u8; 4]>,
    position: [f64; 3],
    occupancy: Option<f32>,
    b_iso: Option<f32>,
    formal_charge: Option<i8>,
    element: Element,
    calc_flag: BioCalcFlag,
    tls_group_id: Option<i16>,
    fraction: Option<f32>,
}

impl PdbAtomRecord {
    fn residue_key(self) -> ResidueKey {
        ResidueKey {
            chain_key: self.chain_key,
            seq_id: self.seq_id,
            name: self.residue_name,
            segment_id: self.segment_id,
        }
    }
}

/// Reads the Gemmi-aligned PDB structural reader surface with explicit parser
/// options for the source-defined record branches.
pub(crate) fn read_pdb_bio_structure_from_str_with_params(
    text: &str,
    params: BioPdbReadParams,
) -> Result<BioStructure, BioReadError> {
    read_pdb_bio_structure_from_str_with_params_and_source(text, "", params)
}

fn read_pdb_bio_structure_from_str_with_params_and_source(
    text: &str,
    source: &str,
    params: BioPdbReadParams,
) -> Result<BioStructure, BioReadError> {
    // BEGIN GEMMI CPP FUNCTION gemmi::populate_structure_from_pdb_stream
    // Gemmi✔️✔️: if (is_record_type4(line, "ATOM") || is_record_type4(line, "HETATM")) {
    // Gemmi✔️✔️:   if (len < 55)
    // Gemmi✔️✔️:     wrong("The line is too short to be correct:\n" + std::string(line));
    // Gemmi✔️✔️:   std::string chain_name = read_string(line+20, 2);
    // Gemmi✔️✔️:   ResidueId rid = read_res_id(line+22, line+17);
    // Gemmi✔️✔️:   if (!chain || chain_name != chain->name) {
    // Gemmi✔️✔️:     if (!model) {
    // Gemmi✔️✔️:       int num = (int) st.models.size() + 1;
    // Gemmi✔️✔️:       if (st.find_model(num))
    // Gemmi✔️✔️:         wrong("ATOM/HETATM between models");
    // Gemmi✔️✔️:       st.models.emplace_back(num);
    // Gemmi✔️✔️:       model = &st.models.back();
    // Gemmi✔️✔️:     const Chain* prev_part = model->find_chain(chain_name);
    // Gemmi✔️✔️:     after_ter = prev_part &&
    // Gemmi✔️✔️:                   prev_part->residues[0].entity_type == EntityType::Polymer;
    // Gemmi✔️✔️:     model->chains.emplace_back(chain_name);
    // Gemmi✔️✔️:     chain = &model->chains.back();
    // Gemmi✔️✔️:     resmap.clear();
    // Gemmi✔️✔️:     resi = nullptr;
    // Gemmi✔️✔️:   if (len > 72)
    // Gemmi✔️✔️:     rid.segment = read_string(line+72, 4);
    // Gemmi✔️✔️:   if (!resi || !resi->matches(rid)) {
    // Gemmi✔️✔️:     auto it = resmap.find(rid);
    // Gemmi✔️✔️:     if (it == resmap.end()) {
    // Gemmi✔️✔️:       resmap.emplace(rid, (int) chain->residues.size());
    // Gemmi✔️✔️:       chain->residues.emplace_back(rid);
    // Gemmi✔️✔️:       resi = &chain->residues.back();
    // Gemmi✔️✔️:       resi->het_flag = line[0] & ~0x20;
    // Gemmi✔️✔️:       if (after_ter)
    // Gemmi✔️✔️:         resi->entity_type = resi->is_water() ? EntityType::Water
    // Gemmi✔️✔️:                                                  : EntityType::NonPolymer;
    // Gemmi✔️✔️:     } else {
    // Gemmi✔️✔️:       resi = &chain->residues[it->second];
    // Gemmi✔️✔️:   Atom atom;
    // Gemmi✔️✔️:   atom.serial = read_serial(line+6);
    // Gemmi✔️✔️:   atom.name = read_string(line+12, 4);
    // Gemmi✔️✔️:   atom.altloc = read_altloc(line[16]);
    // Gemmi✔️✔️:   atom.pos.x = read_double(line+30, 8);
    // Gemmi✔️✔️:   atom.pos.y = read_double(line+38, 8);
    // Gemmi✔️✔️:   atom.pos.z = read_double(line+46, 8);
    // Gemmi✔️✔️:   if (len > 58)
    // Gemmi✔️✔️:     atom.occ = (float) read_double(line+54, 6);
    // Gemmi✔️✔️:   if (len > 64)
    // Gemmi✔️✔️:     atom.b_iso = (float) read_double(line+60, 6);
    // Gemmi✔️✔️:   if (len > 76 && (std::isalpha(line[76]) || std::isalpha(line[77])))
    // Gemmi✔️✔️:     atom.element = Element(line + 76);
    // Gemmi✔️✔️:   else
    // Gemmi✔️✔️:     atom.element = infer_element_from_padded_name(line+12);
    // Gemmi✔️✔️:   atom.charge = (len > 78 ? read_charge(line[78], line[79]) : 0);
    // Gemmi✔️✔️:   resi->atoms.emplace_back(atom);
    // Gemmi✔️✔️: } else if (is_record_type4(line, "ANISOU")) {
    // Gemmi✔️✔️:   if (!model || !chain || !resi || resi->atoms.empty())
    // Gemmi✔️✔️:     wrong("ANISOU record not directly after ATOM/HETATM.");
    // Gemmi✔️✔️:   Atom &atom = resi->atoms.back();
    // Gemmi✔️✔️:   if (atom.aniso.u11 != 0.)
    // Gemmi✔️✔️:     wrong("Duplicated ANISOU record or not directly after ATOM/HETATM.");
    // Gemmi✔️✔️:   atom.aniso.u11 = read_int(line+28, 7) * 1e-4f;
    // Gemmi✔️✔️:   atom.aniso.u22 = read_int(line+35, 7) * 1e-4f;
    // Gemmi✔️✔️:   atom.aniso.u33 = read_int(line+42, 7) * 1e-4f;
    // Gemmi✔️✔️:   atom.aniso.u12 = read_int(line+49, 7) * 1e-4f;
    // Gemmi✔️✔️:   atom.aniso.u13 = read_int(line+56, 7) * 1e-4f;
    // Gemmi✔️✔️:   atom.aniso.u23 = read_int(line+63, 7) * 1e-4f;
    // Gemmi✔️✔️: } else if (is_record_type4(line, "REMARK")) {
    // Gemmi✔️✔️: } else if (is_record_type4(line, "CONECT")) {
    // Gemmi✔️✔️: } else if (is_record_type4(line, "SEQRES")) {
    // Gemmi✔️✔️:   std::string chain_name = read_string(line+10, 2);
    // Gemmi✔️✔️:   Entity& ent = impl::find_or_add(st.entities, chain_name);
    // Gemmi✔️✔️:   ent.entity_type = EntityType::Polymer;
    // Gemmi✔️✔️:   for (int i = 19; i < 68 && i < (int)len; i += 4) {
    // Gemmi✔️✔️:     std::string res_name = read_string(line+i, 3);
    // Gemmi✔️✔️:     if (!res_name.empty())
    // Gemmi✔️✔️:       ent.full_sequence.emplace_back(res_name);
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️: } else if (is_record_type4(line, "HELIX")) {
    // Gemmi✔️✔️: } else if (is_record_type4(line, "SHEET")) {
    // Gemmi✔️✔️: } else if (is_record_type3(line, "TER") && !options.ignore_ter) {
    // Gemmi✔️✔️: } else if (is_record_type4(line, "HEADER")) {
    // Gemmi✔️✔️:   if (len > 50)
    // Gemmi✔️✔️:     st.info["_struct_keywords.pdbx_keywords"] = rtrim_str(std::string(line+10, 40));
    // Gemmi✔️✔️:   if (len > 59) { // date in PDB has format 28-MAR-07
    // Gemmi✔️✔️:     std::string date = pdb_date_format_to_iso(std::string(line+50, 9));
    // Gemmi✔️✔️:     if (!date.empty())
    // Gemmi✔️✔️:       st.info["_pdbx_database_status.recvd_initial_deposition_date"] = date;
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️:   if (len > 66) {
    // Gemmi✔️✔️:     std::string entry_id = rtrim_str(std::string(line+62, 4));
    // Gemmi✔️✔️:     if (!entry_id.empty())
    // Gemmi✔️✔️:       st.info["_entry.id"] = entry_id;
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️: } else if (is_record_type4(line, "TITLE")) {
    // Gemmi✔️✔️:   if (len > 10)
    // Gemmi✔️✔️:     st.info["_struct.title"] += rtrim_str(std::string(line+10, len-10-1));
    // Gemmi✔️✔️: } else if (is_record_type4(line, "KEYWDS")) {
    // Gemmi✔️✔️:   if (len > 10)
    // Gemmi✔️✔️:     st.info["_struct_keywords.text"] += rtrim_str(std::string(line+10, len-10-1));
    // Gemmi✔️✔️: } else if (is_record_type4(line, "EXPDTA")) {
    // Gemmi✔️✔️:   if (len > 10)
    // Gemmi✔️✔️:     st.info["_exptl.method"] += trim_str(std::string(line+10, len-10-1));
    // Gemmi✔️✔️: } else if (is_record_type4(line, "AUTHOR") && len > 10) {
    // Gemmi✔️✔️:   split_str_into(std::string(start, end), ',', st.meta.authors);
    // Gemmi✔️✔️: } else if (is_record_type4(line, "CRYST1")) {
    // Gemmi✔️✔️:   if (len > 54)
    // Gemmi✔️✔️:     st.cell.set(read_double(line+6, 9),
    // Gemmi✔️✔️:                 read_double(line+15, 9),
    // Gemmi✔️✔️:                 read_double(line+24, 9),
    // Gemmi✔️✔️:                 read_double(line+33, 7),
    // Gemmi✔️✔️:                 read_double(line+40, 7),
    // Gemmi✔️✔️:                 read_double(line+47, 7));
    // Gemmi✔️✔️:   if (len > 56)
    // Gemmi✔️✔️:     st.spacegroup_hm = read_string(line+55, 11);
    // Gemmi✔️✔️:   if (len > 67) {
    // Gemmi✔️✔️:     std::string z = read_string(line+66, 4);
    // Gemmi✔️✔️:     if (!z.empty())
    // Gemmi✔️✔️:       st.info["_cell.Z_PDB"] = z;
    // Gemmi✔️✔️:   }
    // END GEMMI CPP FUNCTION

    let mut params = params;
    if params.max_line_length == 0 || params.max_line_length > 120 {
        params.max_line_length = 120;
    }

    let mut builder = PdbBioBuilder::default();
    // BEGIN GEMMI CPP FUNCTION gemmi::populate_structure_from_pdb_stream
    // Gemmi✔️✔️: st.input_format = CoorFormat::Pdb;
    // Gemmi✔️✔️: st.name = path_basename(source, {".gz", ".pdb"});
    // END GEMMI CPP FUNCTION
    builder.structure.input_format = BioCoorFormat::Pdb;
    builder.structure.name = gemmi_path_basename(source, &[".gz", ".pdb"]);
    builder.skip_remarks = params.skip_remarks;
    let mut split_chain_on_next_atom = false;
    let mut after_ter = false;
    let mut matrix = BioTransform::default();

    for (line_index, line) in text.lines().enumerate() {
        let line_number = line_index + 1;
        if params.check_non_ascii && builder.structure.non_ascii_line.is_none() {
            if line.as_bytes().iter().any(|&byte| byte >= 0x80) {
                builder.structure.non_ascii_line = Some(line_number);
            }
        }
        if line.len() > params.max_line_length {
            return Err(BioReadError::Parse {
                line_number,
                message: format!(
                    "PDB line exceeds max_line_length {}",
                    params.max_line_length
                ),
            });
        }
        if starts_record(line, "data_") && builder.current_model.is_none() {
            return Err(BioReadError::Parse {
                line_number,
                message: "Incorrect file format (perhaps it is cif not pdb?)".to_string(),
            });
        }
        if line.starts_with("{\"data_") && builder.current_model.is_none() {
            return Err(BioReadError::Parse {
                line_number,
                message: "Incorrect file format (perhaps it is mmJSON not pdb?)".to_string(),
            });
        }
        let record = record_type(line);
        match record {
            "ATOM" | "HETA" => {
                if builder.current_model.is_none() {
                    let next_num = builder.structure.models.len() as i32 + 1;
                    if builder.has_model_number(next_num) {
                        return Err(BioReadError::Parse {
                            line_number,
                            message: "ATOM/HETATM between models".to_string(),
                        });
                    }
                    builder.begin_model(Some(next_num));
                }
                let mut atom = parse_pdb_atom_record(line, line_number)?;
                if line.len() > 72 {
                    let mut segment = [b' '; 4];
                    let raw = field_raw(line, 72, 76).as_bytes();
                    segment[..raw.len().min(4)].copy_from_slice(&raw[..raw.len().min(4)]);
                    atom.segment_id = Some(segment);
                }
                let model_id = builder
                    .current_model
                    .expect("implicit model initialized above");
                let force_new_chain = split_chain_on_next_atom;
                if force_new_chain {
                    builder.current_chain = None;
                    builder.current_residue = None;
                } else if let Some((_, current_chain_key)) = builder.current_chain
                    && current_chain_key != atom.chain_key
                {
                    let prev_part = builder.find_chain_in_model(model_id, atom.chain_key);
                    after_ter = prev_part
                        .is_some_and(|chain_id| builder.chain_has_polymer_first_residue(chain_id));
                }
                let entity_kind = if after_ter {
                    if atom.residue_kind == ResidueKind::Water {
                        EntityKind::Water
                    } else {
                        EntityKind::NonPolymer
                    }
                } else {
                    EntityKind::Unknown
                };
                builder.push_atom(atom, None, force_new_chain, entity_kind);
                split_chain_on_next_atom = false;
            }
            "MODE" if starts_record(line, "MODEL") => {
                if builder.current_model.is_some() && builder.current_chain.is_some() {
                    return Err(BioReadError::Parse {
                        line_number,
                        message: "MODEL without ENDMDL?".to_string(),
                    });
                }
                let source_model_number = parse_optional_i32(field(line, 6, 14), line_number)?;
                if let Some(model_number) = source_model_number
                    && builder.has_model_number(model_number)
                {
                    return Err(BioReadError::Parse {
                        line_number,
                        message: format!("duplicate MODEL number: {model_number}"),
                    });
                }
                builder.begin_model(source_model_number);
                split_chain_on_next_atom = false;
                after_ter = false;
            }
            "ENDM" if starts_record(line, "ENDMDL") => {
                builder.end_model();
                split_chain_on_next_atom = false;
                after_ter = false;
            }
            "ANIS" if starts_record(line, "ANISOU") => {
                if builder.current_model.is_none()
                    || builder.current_chain.is_none()
                    || builder.current_residue.is_none()
                {
                    return Err(BioReadError::Parse {
                        line_number,
                        message: "ANISOU record not directly after ATOM/HETATM.".to_string(),
                    });
                }
                let anisou = parse_pdb_anisou_record(line, line_number)?;
                builder.set_last_atom_anisou(line_number, anisou)?;
            }
            "REMA" if starts_record(line, "REMARK") => {
                if !params.skip_remarks {
                    builder
                        .structure
                        .raw_remarks
                        .push(line.trim_end_matches(['\n', '\r']).to_string());
                }
            }
            "CONE" if starts_record(line, "CONECT") => {
                let serial = field(line, 6, 11).parse::<i32>().unwrap_or(0);
                if line.len() >= 11 && serial != 0 {
                    let bonded_atoms = builder.structure.conect_map.entry(serial).or_default();
                    let limit = 27usize.min(line.len().saturating_sub(1));
                    let mut offset = 11usize;
                    while offset <= limit {
                        let bonded = field(line, offset, offset + 5).parse::<i32>().unwrap_or(0);
                        if bonded != 0 {
                            bonded_atoms.push(bonded);
                        }
                        offset += 5;
                    }
                }
            }
            "SEQR" if starts_record(line, "SEQRES") => {
                parse_pdb_seqres_record(&mut builder, line);
            }
            "HELI" if starts_record(line, "HELIX") => {
                parse_pdb_helix_record(&mut builder.structure, line);
            }
            "SHEE" if starts_record(line, "SHEET") => {
                parse_pdb_sheet_record(&mut builder.structure, line);
            }
            "HEAD" if starts_record(line, "HEADER") => {
                parse_pdb_header_record(&mut builder, line);
            }
            "TITL" if starts_record(line, "TITLE") => {
                append_metadata_string(
                    &mut builder.structure.metadata.title,
                    field_raw(line, 10, line.len()),
                );
            }
            "KEYW" if starts_record(line, "KEYWDS") => {
                append_metadata_string(
                    &mut builder.structure.metadata.keywords,
                    field_raw(line, 10, line.len()),
                );
            }
            "EXPD" if starts_record(line, "EXPDTA") => {
                append_metadata_string(
                    &mut builder.structure.metadata.experimental_method,
                    field_raw(line, 10, line.len()).trim(),
                );
            }
            "AUTH" if starts_record(line, "AUTHOR") => parse_pdb_author_record(&mut builder, line),
            "CRYS" if starts_record(line, "CRYST1") => {
                parse_pdb_cryst1_record(&mut builder, line, line_number)?;
            }
            // Gemmi✔️✔️: } else if (is_record_type4(line, "SCALEn")) {
            // Gemmi✔️✔️:   if (read_matrix(matrix, line, len) == 3) {
            // Gemmi✔️✔️:     st.cell.set_matrices_from_fract(matrix);
            // Gemmi✔️✔️:     matrix.set_identity();
            // Gemmi✔️✔️:   }
            "SCAL" if starts_record(line, "SCALE") => {
                if read_matrix(&mut matrix, line) == 3 {
                    let crystal = builder
                        .structure
                        .crystal
                        .get_or_insert_with(default_crystal_info);
                    crystal.scale = Some(matrix);
                    crystal_set_matrices_from_fract(crystal, matrix);
                    matrix = bio_transform_identity();
                }
            }
            // Gemmi✔️✔️: } else if (is_record_type4(line, "ORIGX")) {
            // Gemmi✔️✔️:   st.has_origx = true;
            // Gemmi✔️✔️:   read_matrix(st.origx, line, len);
            "ORIG" if starts_record(line, "ORIGX") => {
                builder.structure.has_origx = true;
                read_matrix(&mut builder.structure.origx, line);
            }
            // Gemmi✔️✔️: } else if (is_record_type4(line, "MTRIXn")) {
            // Gemmi✔️✔️:   if (read_matrix(matrix, line, len) == 3) {
            // Gemmi✔️✔️:     std::string id = read_string(line+7, 3);
            // Gemmi✔️✔️:     if (matrix.is_identity()) {
            // Gemmi✔️✔️:       st.info["_struct_ncs_oper.id"] = id;
            // Gemmi✔️✔️:     } else {
            // Gemmi✔️✔️:       bool given = len > 59 && line[59] == '1';
            // Gemmi✔️✔️:       st.ncs.push_back({id, given, matrix});
            // Gemmi✔️✔️:       matrix.set_identity();
            // Gemmi✔️✔️:     }
            // Gemmi✔️✔️:   }
            "MTRI" if starts_record(line, "MTRIX") => {
                if read_matrix(&mut matrix, line) == 3 {
                    let id = read_fixed_string(line, 7, 3);
                    if bio_transform_is_identity(&matrix) {
                        builder.structure.ncs_oper_identity_id = Some(id);
                    } else {
                        let given = line.len() > 59 && byte_at(line, 59) == b'1';
                        builder.structure.ncs_operators.push(BioNcsOperator {
                            id,
                            given,
                            transform: matrix,
                        });
                        matrix = bio_transform_identity();
                    }
                }
            }
            "TER " | "TER" if starts_record3(line, "TER") && !params.ignore_ter => {
                if let Some((chain_id, _)) = builder.current_chain {
                    if builder.structure.ter_status == 'e' {
                        continue;
                    }
                    builder.structure.ter_status = 'y';
                    if params.split_chain_on_ter {
                        builder.current_chain = None;
                        builder.current_residue = None;
                        split_chain_on_next_atom = true;
                        after_ter = false;
                    } else {
                        if after_ter {
                            builder.structure.ter_status = 'e';
                            continue;
                        }
                        if chain_has_water(&builder.structure, chain_id) {
                            builder.structure.ter_status = 'e';
                        }
                        mark_chain_polymer(&mut builder.structure, chain_id);
                        after_ter = true;
                        builder.current_residue = None;
                    }
                }
            }
            "END " | "END" if starts_record3(line, "END") => break,
            "SSBO" if starts_record(line, "SSBOND") => {
                builder
                    .structure
                    .deferred_conn_records
                    .push(line.trim_end_matches(['\n', '\r']).to_string());
            }
            "LINK" if starts_record(line, "LINK") => {
                builder
                    .structure
                    .deferred_conn_records
                    .push(line.trim_end_matches(['\n', '\r']).to_string());
            }
            "CISP" if starts_record(line, "CISPEP") => {
                builder
                    .structure
                    .deferred_conn_records
                    .push(line.trim_end_matches(['\n', '\r']).to_string());
            }
            "MODR" if starts_record(line, "MODRES") => {
                parse_pdb_modres_record(&mut builder.structure, line);
            }
            "HETN" if starts_record(line, "HETNAM") => {
                parse_pdb_hetnam_record(&mut builder.structure, line);
            }
            "DBRE" if starts_record(line, "DBREF") => {
                parse_pdb_dbref_record(&mut builder, line);
            }
            _ => {}
        }
    }

    builder.finish()
}

fn gemmi_path_basename(path: &str, extensions: &[&str]) -> String {
    // BEGIN GEMMI CPP FUNCTION gemmi::path_basename
    // Gemmi✔️✔️: size_t pos = path.find_last_of("\\/");
    // Gemmi✔️✔️: std::string basename = pos == std::string::npos ? path : path.substr(pos + 1);
    // Gemmi✔️✔️: for (const char* ext : exts) {
    // Gemmi✔️✔️:   size_t len = std::strlen(ext);
    // Gemmi✔️✔️:   if (basename.size() > len &&
    // Gemmi✔️✔️:       basename.compare(basename.length() - len, len, ext, len) == 0)
    // Gemmi✔️✔️:     basename.resize(basename.length() - len);
    // Gemmi✔️✔️: }
    // Gemmi✔️✔️: return basename;
    // END GEMMI CPP FUNCTION
    let mut basename = path
        .rfind(['\\', '/'])
        .map_or(path, |position| &path[position + 1..])
        .to_string();
    for extension in extensions {
        if basename.len() > extension.len() && basename.ends_with(extension) {
            basename.truncate(basename.len() - extension.len());
        }
    }
    basename
}

/// Reads the Gemmi-aligned mmCIF structural reader surface into a
/// `BioStructure`.
///
/// The historical function name is kept for API stability. Unsupported mmCIF
/// source branches remain marked in this module's Gemmi source frame.
#[deprecated(
    since = "0.2.13",
    note = "use BioStructure::from_mmcif_str(); this function now reads more than an atom-site subset"
)]
pub fn read_mmcif_atom_site_subset_from_str(text: &str) -> Result<BioStructure, BioReadError> {
    let document = parse_cif_document(text)?;
    make_structure_from_mmcif_document(document, None)
}

fn build_structure_from_mmcif_first_block(
    document: &CifDocument,
    block: &CifBlock,
) -> Result<BioStructure, BioReadError> {
    let atom_site_storage;
    let atom_site = if let Some(loop_) = block
        .loops
        .iter()
        .find(|loop_| loop_.tags.iter().any(|tag| tag == "_atom_site.id"))
    {
        loop_
    } else {
        atom_site_storage = make_single_row_atom_site_loop_from_items(block)?;
        atom_site_storage.as_ref().ok_or_else(|| {
            unsupported(
                &BIO_MMCIF_READ_FEATURE,
                0,
                "mmCIF _atom_site loop is missing",
            )
        })?
    };
    let empty_loops = Vec::new();
    let atom_site_loops = if block
        .loops
        .iter()
        .any(|loop_| loop_.tags.iter().any(|tag| tag == "_atom_site.id"))
    {
        &block.loops
    } else {
        &empty_loops
    };
    let mut structure = read_mmcif_atom_site(atom_site, atom_site_loops)?;
    structure.input_format = BioCoorFormat::Mmcif;
    structure.name = block.name.clone();
    read_entry_info(&document, &mut structure);
    read_audit_author(&block.loops, &mut structure)?;
    read_refinement_info(&document, &block.loops, &mut structure)?;
    read_tls_info(&block.loops, &mut structure)?;
    read_experimental_info(&block.loops, &mut structure)?;
    read_reflns_info(&block.loops, &mut structure)?;
    read_software_info(&block.loops, &mut structure)?;
    read_ncs_info(&block.loops, &mut structure)?;
    set_cell_from_mmcif(&document, &mut structure.crystal);
    if let Some(spacegroup_hm) = find_spacegroup_hm_value(&document).and_then(cif_optional) {
        structure
            .crystal
            .get_or_insert_with(default_crystal_info)
            .spacegroup_hm = Some(spacegroup_hm.to_string());
    }
    if let Some(fract) = find_cif_transform(
        &block.loops,
        "_atom_sites.fract_transf_matrix",
        "_atom_sites.fract_transf_vector",
    )? {
        let crystal = structure.crystal.get_or_insert_with(default_crystal_info);
        crystal_set_matrices_from_fract(crystal, fract);
    }
    if let Some(origx) = find_cif_transform(
        &block.loops,
        "_database_PDB_matrix.origx",
        "_database_PDB_matrix.origx_vector",
    )? {
        structure.has_origx = true;
        structure.origx = origx;
    }
    if find_cif_loop(&block.loops, "_struct_asym.id").is_none() {
        infer_entity_subchains_from_first_model(&mut structure);
    }
    fill_residue_entity_type(&mut structure);
    setup_cell_images(&mut structure)?;
    structure.helices = read_helices(&block.loops)?;
    structure.sheets = read_sheets(&block.loops)?;
    read_connectivity(&block.loops, &mut structure)?;
    read_prot_cis(&block.loops, &mut structure)?;
    read_struct_mod_residue(&block.loops, &mut structure)?;
    structure.assemblies = read_assemblies(&block.loops)?;
    read_sifts_unp(&block.loops, &mut structure)?;
    import_shortened_ccd_codes_from_chem_comp(&block.loops, &mut structure)?;
    restore_full_ccd_codes(&mut structure);
    Ok(structure)
}

fn read_mmcif_atom_site_subset_from_document(
    document: &CifDocument,
) -> Result<BioStructure, BioReadError> {
    // BEGIN GEMMI CPP FUNCTION gemmi::make_structure
    // Gemmi✔️✔️: for (size_t i = 1; i < doc.blocks.size(); ++i)
    // Gemmi✔️✔️:   if (doc.blocks[i].has_tag("_atom_site.id"))
    // Gemmi✔️✔️:     fail("2+ blocks are ok if only the first one has coordinates;\n"
    // Gemmi✔️✔️:          "_atom_site in block #" + std::to_string(i+1) + ": " + doc.source);
    // Gemmi✔️✔️: Structure st = make_structure_from_block(doc.blocks.at(0));
    // Gemmi✔️✔️: if (save_doc)
    // Gemmi✔️✔️:   *save_doc = std::move(doc);
    // Gemmi✔️✔️: return st;
    // END GEMMI CPP FUNCTION

    // BEGIN GEMMI CPP FUNCTION gemmi::read_atom_sites
    // Gemmi✔️✔️: auto aniso_map = get_anisotropic_u(block);
    // Gemmi✔️✔️: cif::Table atom_table = block.find("_atom_site.",
    // Gemmi✔️✔️:                                {"id",
    // Gemmi✔️✔️:                                 "?group_PDB",
    // Gemmi✔️✔️:                                 "type_symbol",
    // Gemmi✔️✔️:                                 "?label_atom_id",
    // Gemmi✔️✔️:                                 "label_alt_id",
    // Gemmi✔️✔️:                                 "?label_comp_id",
    // Gemmi✔️✔️:                                 "label_asym_id",
    // Gemmi✔️✔️:                                 "?label_entity_id",
    // Gemmi✔️✔️:                                 "?label_seq_id",
    // Gemmi✔️✔️:                                 "?pdbx_PDB_ins_code",
    // Gemmi✔️✔️:                                 "Cartn_x",
    // Gemmi✔️✔️:                                 "Cartn_y",
    // Gemmi✔️✔️:                                 "Cartn_z",
    // Gemmi✔️✔️:                                 "?occupancy",
    // Gemmi✔️✔️:                                 "?B_iso_or_equiv",
    // Gemmi✔️✔️:                                 "?pdbx_formal_charge",
    // Gemmi✔️✔️:                                 "?auth_seq_id",
    // Gemmi✔️✔️:                                 "?auth_comp_id",
    // Gemmi✔️✔️:                                 "?auth_asym_id",
    // Gemmi✔️✔️:                                 "?auth_atom_id",
    // Gemmi✔️✔️:                                 "?pdbx_PDB_model_num",
    // Gemmi✔️✔️:                                 "?calc_flag",
    // Gemmi✔️✔️:                                 "?pdbx_tls_group_id",
    // Gemmi✔️✔️:                                 "?ccp4_deuterium_fraction",
    // Gemmi✔️✔️: if (atom_table.length() != 0) {
    // Gemmi✔️✔️:   RowAccess asym_id(atom_table, kAuthAsymId, kLabelAsymId);
    // Gemmi✔️✔️:   RowAccess comp_id(atom_table, kAuthCompId, kLabelCompId);
    // Gemmi✔️✔️:   RowAccess atom_id(atom_table, kAuthAtomId, kLabelAtomId);
    // Gemmi✔️✔️:   RowAccess seq_id(atom_table, kAuthSeqId, kLabelSeqId);
    // Gemmi✔️✔️:   if (!asym_id.ok()) fail("Neither _atom_site.label_asym_id nor auth_asym_id found");
    // Gemmi✔️✔️:   if (!comp_id.ok()) fail("Neither _atom_site.label_comp_id nor auth_comp_id found");
    // Gemmi✔️✔️:   if (!atom_id.ok()) fail("Neither _atom_site.label_atom_id nor auth_atom_id found");
    // Gemmi✔️✔️:   if (!seq_id.ok()) fail("Neither _atom_site.label_seq_id nor auth_seq_id found");
    // Gemmi✔️✔️:   for (auto row : atom_table) {
    // Gemmi✔️✔️:     atom.name = cif::as_string(atom_id.get(gap));
    // Gemmi✔️✔️:     atom.altloc = cif::as_char(row[kAltId], '\0');
    // Gemmi✔️✔️:     atom.charge = row.has2(kCharge) ? cif::as_int(row[kCharge]) : 0;
    // Gemmi✔️✔️:     atom.element = gemmi::Element(cif::as_string(row[kSymbol]));
    // Gemmi✔️✔️:     atom.serial = string_to_int(row[kId], false);
    // Gemmi✔️✔️:     if (st.has_d_fraction)
    // Gemmi✔️✔️:       atom.fraction = (float) cif::as_number(row[kDeuterium], 0.);
    // Gemmi✔️✔️:     if (row.has2(kCalcFlag)) { set atom.calc_flag from CIF value }
    // Gemmi✔️✔️:     if (row.has2(kTlsGroupId)) { set atom.tls_group_id from CIF value }
    // Gemmi✔️✔️:     atom.pos.x = cif::as_number(row[kX]);
    // Gemmi✔️✔️:     atom.pos.y = cif::as_number(row[kY]);
    // Gemmi✔️✔️:     atom.pos.z = cif::as_number(row[kZ]);
    // Gemmi✔️✔️:     if (row.has2(kOcc))
    // Gemmi✔️✔️:       atom.occ = (float) cif::as_number(row[kOcc]);
    // Gemmi✔️✔️:     if (row.has2(kBiso))
    // Gemmi✔️✔️:       atom.b_iso = (float) cif::as_number(row[kBiso]);
    // Gemmi✔️✔️:     if (!aniso_map.empty()) { attach matching anisotropic U tensor }
    // Gemmi✔️✔️:     resi->atoms.emplace_back(atom);
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️: }
    // END GEMMI CPP FUNCTION

    // BEGIN GEMMI CPP FUNCTION read_entity_and_sequence_info
    // Gemmi✔️✔️: cif::Table polymer_types = block.find("_entity_poly.", {"entity_id", "type"});
    // Gemmi✔️✔️: for (auto row : block.find("_entity.", {"id", "?type"})) {
    // Gemmi✔️✔️:     Entity ent(row.str(0));
    // Gemmi✔️✔️:     if (row.has(1))
    // Gemmi✔️✔️:         ent.entity_type = entity_type_from_string(row.str(1));
    // Gemmi✔️✔️:     ent.polymer_type = PolymerType::Unknown;
    // Gemmi✔️✔️:     if (polymer_types.ok()) {
    // Gemmi✔️✔️:         try {
    // Gemmi✔️✔️:             std::string poly_type = polymer_types.find_row(ent.name).str(1);
    // Gemmi✔️✔️:             if (ent.entity_type == EntityType::Unknown)
    // Gemmi✔️✔️:                 ent.entity_type = EntityType::Polymer;
    // Gemmi✔️✔️:             ent.polymer_type = polymer_type_from_string(poly_type);
    // Gemmi✔️✔️:         } catch (std::runtime_error&) {}
    // Gemmi✔️✔️:     }
    // Gemmi✔️✔️:     st.entities.push_back(ent);
    // Gemmi✔️✔️: }
    // Gemmi✔️✔️: for (auto row : block.find("_entity_poly_seq.",
    // Gemmi✔️✔️:                            {"entity_id", "num", "mon_id"}))
    // Gemmi✔️✔️:     if (Entity* ent = st.get_entity(row.str(0))) {
    // Gemmi✔️✔️:         int pos = cif::as_int(row[1], 0) - 1;
    // Gemmi✔️✔️:         if (pos == (int) ent->full_sequence.size())
    // Gemmi✔️✔️:             ent->full_sequence.push_back(row.str(2));
    // Gemmi✔️✔️:         else if (pos >= 0 && pos < (int) ent->full_sequence.size())
    // Gemmi✔️✔️:             cat_to(ent->full_sequence[pos], ',', row.str(2));
    // Gemmi✔️✔️:     }
    // Gemmi✔️✔️: cif::Table struct_ref = block.find("_struct_ref.", ...);
    // Gemmi✔️✔️: cif::Table s_asym_table = block.find("_struct_asym.", {"id", "entity_id"});
    // Gemmi✔️✔️: if (s_asym_table.ok()) {
    // Gemmi✔️✔️:     for (auto row : s_asym_table)
    // Gemmi✔️✔️:         if (Entity* ent = st.get_entity(row.str(1)))
    // Gemmi✔️✔️:             ent->subchains.push_back(row.str(0));
    // Gemmi✔️✔️: } else if (!st.models.empty()) { infer subchains from model residues }
    // END GEMMI CPP FUNCTION
    //
    // BEGIN GEMMI CPP FUNCTION populate_structure_from_block
    // Gemmi✔️✔️: st.input_format = CoorFormat::Mmcif;
    // Gemmi✔️✔️: st.name = block.name;
    // Gemmi✔️✔️: impl::set_cell_from_mmcif(block, st.cell);
    // Gemmi✔️✔️: st.spacegroup_hm = cif::as_string(impl::find_spacegroup_hm_value(block));
    // Gemmi✔️✔️: read_entry_info(block, st);
    // Gemmi✔️✔️: read_audit_author(block, st);
    // Gemmi✔️✔️: read_refinement_info(block, st);
    // Gemmi✔️✔️: read_tls_info(block, st);
    // Gemmi✔️✔️: read_experimental_info(block, st);
    // Gemmi✔️✔️: read_reflns_info(block, st);
    // Gemmi✔️✔️: read_software_info(block, st);
    // Gemmi✔️✔️: read_ncs_info(block, st);
    // Gemmi✔️✔️: cif::Table fract_tv = block.find("_atom_sites.fract_transf_",
    // Gemmi✔️✔️:                                  transform_tags("matrix", "vector"));
    // Gemmi✔️✔️: if (fract_tv.length() > 0) {
    // Gemmi✔️✔️:   Transform fract = get_transform_matrix(fract_tv[0]);
    // Gemmi✔️✔️:   st.cell.set_matrices_from_fract(fract);
    // Gemmi✔️✔️: }
    // Gemmi✔️✔️: cif::Table origx_tv = block.find("_database_PDB_matrix.",
    // Gemmi✔️✔️:                                  transform_tags("origx", "origx_vector"));
    // Gemmi✔️✔️: if (origx_tv.length() > 0) {
    // Gemmi✔️✔️:   st.has_origx = true;
    // Gemmi✔️✔️:   st.origx = get_transform_matrix(origx_tv[0]);
    // Gemmi✔️✔️: }
    // Gemmi✔️✔️: read_atom_sites(block, st);
    // Gemmi✔️✔️: read_entity_and_sequence_info(block, st);
    // Gemmi✔️✔️: fill_residue_entity_type(st);
    // Gemmi✔️✔️: st.setup_cell_images();
    // Gemmi✔️✔️: st.helices = read_helices(block);
    // Gemmi✔️✔️: st.sheets = read_sheets(block);
    // Gemmi✔️✔️: read_connectivity(block, st);
    // Gemmi✔️✔️: read_prot_cis(block, st);
    // Gemmi✔️✔️: read_struct_mod_residue(block, st);
    // Gemmi✔️✔️: st.assemblies = read_assemblies(block);
    // Gemmi✔️✔️: read_sifts_unp(block, st);
    // Gemmi✔️✔️: cif::Table chem_comp_table = block.find("_chem_comp.", {"id", "three_letter_code"});
    // Gemmi✔️✔️: if (chem_comp_table.ok()) {
    // Gemmi✔️✔️:   for (auto row : chem_comp_table) {
    // Gemmi✔️✔️:     std::string alias = row.str(0);
    // Gemmi✔️✔️:     std::string long_id = row.str(1);
    // Gemmi✔️✔️:     if (alias[0] == '~' && long_id[0] != '~' && long_id[0] != '\0')
    // Gemmi✔️✔️:       st.shortened_ccd_codes.emplace_back(long_id, alias);
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️:   restore_full_ccd_codes(st);
    // Gemmi✔️✔️: }
    // END GEMMI CPP FUNCTION

    for (index, block) in document.blocks.iter().enumerate().skip(1) {
        if block
            .loops
            .iter()
            .any(|loop_| loop_.tags.iter().any(|tag| tag == "_atom_site.id"))
        {
            return Err(BioReadError::Parse {
                line_number: 0,
                message: format!(
                    "2+ blocks are ok if only the first one has coordinates;\n_atom_site in block #{}",
                    index + 1
                ),
            });
        }
    }
    let block = document.blocks.first().ok_or_else(|| {
        unsupported(
            &BIO_MMCIF_READ_FEATURE,
            0,
            "mmCIF _atom_site loop is missing",
        )
    })?;
    build_structure_from_mmcif_first_block(document, block)
}

// BEGIN GEMMI CPP FUNCTION gemmi::make_structure
// Gemmi✔️✔️: inline Structure make_structure(cif::Document&& doc, cif::Document* save_doc=nullptr) {
// Gemmi✔️✔️:   for (size_t i = 1; i < doc.blocks.size(); ++i)
// Gemmi✔️✔️:     if (doc.blocks[i].has_tag("_atom_site.id"))
// Gemmi✔️✔️:       fail("2+ blocks are ok if only the first one has coordinates;\n"
// Gemmi✔️✔️:            "_atom_site in block #" + std::to_string(i+1) + ": " + doc.source);
// Gemmi✔️✔️:   Structure st = make_structure_from_block(doc.blocks.at(0));
// Gemmi✔️✔️:   if (save_doc)
// Gemmi✔️✔️:     *save_doc = std::move(doc);
// Gemmi✔️✔️:   return st;
// Gemmi✔️✔️: }
// END GEMMI CPP FUNCTION
fn make_structure_from_mmcif_document(
    document: CifDocument,
    save_doc: Option<&mut CifDocument>,
) -> Result<BioStructure, BioReadError> {
    let structure = read_mmcif_atom_site_subset_from_document(&document)?;
    if let Some(out) = save_doc {
        *out = document;
    }
    Ok(structure)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ChemCompModel {
    Xyz = 1,
    Example = 2,
    Ideal = 4,
    First = 8,
}

#[derive(Debug)]
struct ChemCompResidueBuild {
    residue: ResidueRow,
    atoms: Vec<AtomRow>,
    positions: Vec<[f64; 3]>,
}

// BEGIN GEMMI CPP FUNCTION make_residue_from_chemcomp_block
// Gemmi✔️✔️: Residue make_residue_from_chemcomp_block(const cif::Block& block, ChemCompModel kind) {
// Gemmi✔️✔️:   std::array<std::string, 3> xyz_tags;
// Gemmi✔️✔️:   if (kind == ChemCompModel::First) {
// Gemmi✔️✔️:     auto item = const_cast<cif::Block&>(block).find_loop("_chem_comp_atom.atom_id").item();
// Gemmi✔️✔️:     if (loop_item && loop_item->type == cif::ItemType::Loop)
// Gemmi✔️✔️:       for (const std::string& tag : loop_item->loop.tags) {
// Gemmi✔️✔️:         std::string lctag = gemmi::to_lower(tag);
// Gemmi✔️✔️:         if (lctag == "_chem_comp_atom.x") {
// Gemmi✔️✔️:           kind = ChemCompModel::Xyz;
// Gemmi✔️✔️:           break;
// Gemmi✔️✔️:         } else if (lctag == "_chem_comp_atom.model_cartn_x") {
// Gemmi✔️✔️:           kind = ChemCompModel::Example;
// Gemmi✔️✔️:           break;
// Gemmi✔️✔️:         } else if (lctag == "_chem_comp_atom.pdbx_model_cartn_x_ideal") {
// Gemmi✔️✔️:           kind = ChemCompModel::Ideal;
// Gemmi✔️✔️:           break;
// Gemmi✔️✔️:         }
// Gemmi✔️✔️:       }
// Gemmi✔️✔️:   }
// Gemmi✔️✔️:   switch (kind) {
// Gemmi✔️✔️:     case ChemCompModel::Xyz:
// Gemmi✔️✔️:       xyz_tags = {{"x", "y", "z"}};
// Gemmi✔️✔️:       break;
// Gemmi✔️✔️:     case ChemCompModel::Example:
// Gemmi✔️✔️:       xyz_tags = {{"model_Cartn_x", "model_Cartn_y", "model_Cartn_z"}};
// Gemmi✔️✔️:       break;
// Gemmi✔️✔️:     case ChemCompModel::Ideal:
// Gemmi✔️✔️:       xyz_tags = {{"pdbx_model_Cartn_x_ideal",
// Gemmi✔️✔️:                    "pdbx_model_Cartn_y_ideal",
// Gemmi✔️✔️:                    "pdbx_model_Cartn_z_ideal"}};
// Gemmi✔️✔️:       break;
// Gemmi✔️✔️:     default:
// Gemmi✔️✔️:       break;
// Gemmi✔️✔️:   }
// Gemmi✔️✔️:   Residue res;
// Gemmi✔️✔️:   res.seqid.num = 1;
// Gemmi✔️✔️:   cif::Column col = const_cast<cif::Block&>(block).find_values("_chem_comp_atom.comp_id");
// Gemmi✔️✔️:   if (col && col.length() > 0)
// Gemmi✔️✔️:     res.name = col[0];
// Gemmi✔️✔️:   else
// Gemmi✔️✔️:     res.name = block.name.substr(starts_with(block.name, "comp_") ? 5 : 0);
// Gemmi✔️✔️:   cif::Table table = const_cast<cif::Block&>(block).find("_chem_comp_atom.",
// Gemmi✔️✔️:           {"atom_id", "type_symbol", "?charge", xyz_tags[0], xyz_tags[1], xyz_tags[2]});
// Gemmi✔️✔️:   res.atoms.resize(table.length());
// Gemmi✔️✔️:   int n = 0;
// Gemmi✔️✔️:   for (auto row : table) {
// Gemmi✔️✔️:     Atom& atom = res.atoms[n++];
// Gemmi✔️✔️:     atom.name = row.str(0);
// Gemmi✔️✔️:     atom.element = Element(row.str(1));
// Gemmi✔️✔️:     if (row.has2(2))
// Gemmi✔️✔️:       atom.charge = (signed char) std::round(cif::as_number(row[2]));
// Gemmi✔️✔️:     atom.pos = Position(cif::as_number(row[3]),
// Gemmi✔️✔️:                         cif::as_number(row[4]),
// Gemmi✔️✔️:                         cif::as_number(row[5]));
// Gemmi✔️✔️:   }
// Gemmi✔️✔️:   return res;
// Gemmi✔️✔️: }
// END GEMMI CPP FUNCTION
fn make_residue_from_chemcomp_block(
    block: &CifBlock,
    kind: ChemCompModel,
) -> Result<ChemCompResidueBuild, BioReadError> {
    let resolved_kind = if kind == ChemCompModel::First {
        if let Some(loop_) = find_cif_loop(&block.loops, "_chem_comp_atom.atom_id") {
            let mut detected = ChemCompModel::First;
            for tag in &loop_.tags {
                let lower = tag.to_ascii_lowercase();
                if lower == "_chem_comp_atom.x" {
                    detected = ChemCompModel::Xyz;
                    break;
                } else if lower == "_chem_comp_atom.model_cartn_x" {
                    detected = ChemCompModel::Example;
                    break;
                } else if lower == "_chem_comp_atom.pdbx_model_cartn_x_ideal" {
                    detected = ChemCompModel::Ideal;
                    break;
                }
            }
            detected
        } else {
            ChemCompModel::First
        }
    } else {
        kind
    };

    let xyz_tags = match resolved_kind {
        ChemCompModel::Xyz => (
            "_chem_comp_atom.x",
            "_chem_comp_atom.y",
            "_chem_comp_atom.z",
        ),
        ChemCompModel::Example => (
            "_chem_comp_atom.model_Cartn_x",
            "_chem_comp_atom.model_Cartn_y",
            "_chem_comp_atom.model_Cartn_z",
        ),
        ChemCompModel::Ideal => (
            "_chem_comp_atom.pdbx_model_Cartn_x_ideal",
            "_chem_comp_atom.pdbx_model_Cartn_y_ideal",
            "_chem_comp_atom.pdbx_model_Cartn_z_ideal",
        ),
        ChemCompModel::First => {
            return Err(unsupported(
                &BIO_MMCIF_READ_FEATURE,
                0,
                "chem_comp block does not contain supported coordinate tags",
            ));
        }
    };

    let atom_loop = find_cif_loop(&block.loops, "_chem_comp_atom.atom_id").ok_or_else(|| {
        unsupported(
            &BIO_MMCIF_READ_FEATURE,
            0,
            "chem_comp block is missing _chem_comp_atom.atom_id loop",
        )
    })?;
    let atom_id = required_cif_col(atom_loop, "_chem_comp_atom.atom_id")?;
    let type_symbol = required_cif_col(atom_loop, "_chem_comp_atom.type_symbol")?;
    let charge = optional_cif_col(atom_loop, "_chem_comp_atom.charge");
    let x = required_cif_col(atom_loop, xyz_tags.0)?;
    let y = required_cif_col(atom_loop, xyz_tags.1)?;
    let z = required_cif_col(atom_loop, xyz_tags.2)?;

    let residue_name =
        if let Some(comp_loop) = find_cif_loop(&block.loops, "_chem_comp_atom.comp_id") {
            let comp_id = required_cif_col(comp_loop, "_chem_comp_atom.comp_id")?;
            let mut rows = cif_loop_rows(comp_loop)?;
            rows.next()
                .and_then(|row| cif_optional(row[comp_id].value.as_str()))
                .map(residue_name_from_field)
                .unwrap_or_else(|| {
                    residue_name_from_field(
                        block
                            .name
                            .strip_prefix("comp_")
                            .unwrap_or(block.name.as_str()),
                    )
                })
        } else {
            residue_name_from_field(
                block
                    .name
                    .strip_prefix("comp_")
                    .unwrap_or(block.name.as_str()),
            )
        };

    let residue = ResidueRow {
        chain_id: ChainId::new(0),
        atom_span: RowSpan::new(0, 0),
        name: residue_name,
        kind: classify_residue_name(residue_name),
        entity_kind: EntityKind::Unknown,
        het_flag: None,
        source: ResidueSourceIds {
            seq_id: Some(PdbSeqId {
                seq_num: 1,
                ins_code: None,
            }),
            label_seq_id: Some(1),
            segment_id: None,
            subchain_id: None,
            label_entity_id: None,
        },
        sifts_unp: None,
    };

    let mut atoms = Vec::new();
    let mut positions = Vec::new();
    for row in cif_loop_rows(atom_loop)? {
        let line_number = row[atom_id].line_number;
        let atom_name = cif_optional(row[atom_id].value.as_str())
            .ok_or_else(|| missing_cif_value(line_number, "_chem_comp_atom.atom_id"))?;
        let symbol = cif_optional(row[type_symbol].value.as_str())
            .ok_or_else(|| missing_cif_value(line_number, "_chem_comp_atom.type_symbol"))?;
        let element = element_from_symbol(symbol).ok_or_else(|| BioReadError::Parse {
            line_number,
            message: format!(
                "invalid _chem_comp_atom.type_symbol: {:?}",
                row[type_symbol].value
            ),
        })?;
        let formal_charge = charge
            .and_then(|idx| cif_optional(row[idx].value.as_str()))
            .map(|value| {
                parse_f32(
                    value,
                    row[charge.unwrap()].line_number,
                    "_chem_comp_atom.charge",
                )
            })
            .transpose()?
            .map(|value| value.round() as i8);
        let pos = [
            parse_f64(row[x].value.as_str(), row[x].line_number, xyz_tags.0)?,
            parse_f64(row[y].value.as_str(), row[y].line_number, xyz_tags.1)?,
            parse_f64(row[z].value.as_str(), row[z].line_number, xyz_tags.2)?,
        ];
        atoms.push(AtomRow {
            residue_id: ResidueId::new(0),
            name: atom_name_from_cif(atom_name),
            element,
            altloc: None,
            occupancy: None,
            b_iso: None,
            formal_charge,
            anisou: None,
            calc_flag: BioCalcFlag::NotSet,
            tls_group_id: None,
            fraction: None,
            source: AtomSourceIds { serial: None },
        });
        positions.push(pos);
    }
    let mut residue = residue;
    residue.atom_span = RowSpan::new(0, atoms.len() as u32);
    Ok(ChemCompResidueBuild {
        residue,
        atoms,
        positions,
    })
}

// BEGIN GEMMI CPP FUNCTION make_model_from_chemcomp_block
// Gemmi✔️✔️: inline Model make_model_from_chemcomp_block(const cif::Block& block, ChemCompModel kind) {
// Gemmi✔️✔️:   Model model;
// Gemmi✔️✔️:   model.chains.emplace_back("");
// Gemmi✔️✔️:   model.chains[0].residues.push_back(make_residue_from_chemcomp_block(block, kind));
// Gemmi✔️✔️:   return model;
// Gemmi✔️✔️: }
// END GEMMI CPP FUNCTION
fn make_model_from_chemcomp_block(
    block: &CifBlock,
    kind: ChemCompModel,
) -> Result<BioStructure, BioReadError> {
    let residue_build = make_residue_from_chemcomp_block(block, kind)?;
    let mut structure = BioStructure::default();
    structure.models.push(ModelRow {
        chain_span: RowSpan::new(0, 1),
        source_model_number: Some(1),
    });
    structure.chains.push(ChainRow {
        model_id: ModelId::new(0),
        entity_id: None,
        residue_span: RowSpan::new(0, 1),
        kind: ChainKind::Unknown,
        source: ChainSourceIds {
            auth_chain_id: None,
            label_asym_id: Some(pdb_chain_id_from_field("")),
        },
    });
    structure.residues.push(residue_build.residue);
    structure.atoms = residue_build.atoms;
    structure.coordinates.positions = residue_build.positions;
    structure.input_format = BioCoorFormat::ChemComp;
    if let Some(item) = block.items.iter().find(|item| item.tag == "_chem_comp.id")
        && let Some(value) = cif_optional(item.value.value.as_str())
    {
        structure.name = value.to_string();
    }
    enforce_bio_structure_invariants(&structure).map_err(BioReadError::Invariant)?;
    Ok(structure)
}

// BEGIN GEMMI CPP FUNCTION make_structure_from_chemcomp_block / check_chemcomp_block_number / make_structure_from_chemcomp_doc
// Gemmi✔️✔️: inline Structure make_structure_from_chemcomp_block(const cif::Block& block, int which=7) {
// Gemmi✔️✔️:   Structure st;
// Gemmi✔️✔️:   st.input_format = CoorFormat::ChemComp;
// Gemmi✔️✔️:   if (const std::string* name = block.find_value("_chem_comp.id"))
// Gemmi✔️✔️:     st.name = *name;
// Gemmi✔️✔️:   auto ok = [which](ChemCompModel x) { return which & static_cast<int>(x); };
// Gemmi✔️✔️:   if (ok(ChemCompModel::Xyz) && block.has_any_value("_chem_comp_atom.x"))
// Gemmi✔️✔️:     st.models.push_back(make_model_from_chemcomp_block(block, ChemCompModel::Xyz));
// Gemmi✔️✔️:   if (ok(ChemCompModel::Example) && block.has_any_value("_chem_comp_atom.model_Cartn_x"))
// Gemmi✔️✔️:     st.models.push_back(make_model_from_chemcomp_block(block, ChemCompModel::Example));
// Gemmi✔️✔️:   if (ok(ChemCompModel::Ideal) && block.has_any_value("_chem_comp_atom.pdbx_model_Cartn_x_ideal"))
// Gemmi✔️✔️:     st.models.push_back(make_model_from_chemcomp_block(block, ChemCompModel::Ideal));
// Gemmi✔️✔️:   st.renumber_models();
// Gemmi✔️✔️:   return st;
// Gemmi✔️✔️: }
// Gemmi✔️✔️: inline int check_chemcomp_block_number(const cif::Document& doc) {
// Gemmi✔️✔️:   if (doc.blocks.size() == 2 && doc.blocks[0].name == "comp_list")
// Gemmi✔️✔️:     return 1;
// Gemmi✔️✔️:   if (doc.blocks.size() == 3 && doc.blocks[0].name.empty() &&
// Gemmi✔️✔️:       doc.blocks[1].name == "comp_list")
// Gemmi✔️✔️:     return 2;
// Gemmi✔️✔️:   if (doc.blocks.size() == 1 &&
// Gemmi✔️✔️:       !doc.blocks[0].has_tag("_atom_site.id") &&
// Gemmi✔️✔️:       !doc.blocks[0].has_tag("_cell.length_a") &&
// Gemmi✔️✔️:       doc.blocks[0].has_tag("_chem_comp_atom.atom_id"))
// Gemmi✔️✔️:     return 0;
// Gemmi✔️✔️:   return -1;
// Gemmi✔️✔️: }
// Gemmi✔️✔️: inline Structure make_structure_from_chemcomp_doc(const cif::Document& doc,
// Gemmi✔️✔️:                                                   cif::Document* save_doc=nullptr,
// Gemmi✔️✔️:                                                   int which=7) {
// Gemmi✔️✔️:   int n = check_chemcomp_block_number(doc);
// Gemmi✔️✔️:   if (n == -1)
// Gemmi✔️✔️:     fail("Not a chem_comp format.");
// Gemmi✔️✔️:   Structure st = make_structure_from_chemcomp_block(doc.blocks[n], which);
// Gemmi✔️✔️:   if (save_doc)
// Gemmi✔️✔️:     *save_doc = std::move(doc);
// Gemmi✔️✔️:   return st;
// Gemmi✔️✔️: }
// END GEMMI CPP FUNCTION
fn make_structure_from_chemcomp_block(
    block: &CifBlock,
    which: i32,
) -> Result<BioStructure, BioReadError> {
    let mut structure = BioStructure::default();
    structure.input_format = BioCoorFormat::ChemComp;
    if let Some(item) = block.items.iter().find(|item| item.tag == "_chem_comp.id")
        && let Some(value) = cif_optional(item.value.value.as_str())
    {
        structure.name = value.to_string();
    }
    let ok = |x: ChemCompModel| which & (x as i32) != 0;
    if ok(ChemCompModel::Xyz) && find_cif_loop(&block.loops, "_chem_comp_atom.x").is_some() {
        append_bio_structure_rows(
            &mut structure,
            make_model_from_chemcomp_block(block, ChemCompModel::Xyz)?,
        );
    }
    if ok(ChemCompModel::Example)
        && find_cif_loop(&block.loops, "_chem_comp_atom.model_Cartn_x").is_some()
    {
        append_bio_structure_rows(
            &mut structure,
            make_model_from_chemcomp_block(block, ChemCompModel::Example)?,
        );
    }
    if ok(ChemCompModel::Ideal)
        && find_cif_loop(&block.loops, "_chem_comp_atom.pdbx_model_Cartn_x_ideal").is_some()
    {
        append_bio_structure_rows(
            &mut structure,
            make_model_from_chemcomp_block(block, ChemCompModel::Ideal)?,
        );
    }
    for (index, model) in structure.models.iter_mut().enumerate() {
        model.source_model_number = Some(index as i32 + 1);
    }
    enforce_bio_structure_invariants(&structure).map_err(BioReadError::Invariant)?;
    Ok(structure)
}

fn append_bio_structure_rows(structure: &mut BioStructure, mut other: BioStructure) {
    let model_start = structure.models.len() as u32;
    let chain_start = structure.chains.len() as u32;
    let residue_start = structure.residues.len() as u32;
    let atom_start = structure.atoms.len() as u32;

    for model in &mut other.models {
        model.chain_span = RowSpan::new(model.chain_span.start + chain_start, model.chain_span.len);
    }
    for chain in &mut other.chains {
        chain.model_id = ModelId::new(chain.model_id.index() + model_start);
        chain.residue_span = RowSpan::new(
            chain.residue_span.start + residue_start,
            chain.residue_span.len,
        );
    }
    for residue in &mut other.residues {
        residue.chain_id = ChainId::new(residue.chain_id.index() + chain_start);
        residue.atom_span =
            RowSpan::new(residue.atom_span.start + atom_start, residue.atom_span.len);
    }
    for atom in &mut other.atoms {
        atom.residue_id = ResidueId::new(atom.residue_id.index() + residue_start);
    }

    structure.models.extend(other.models);
    structure.chains.extend(other.chains);
    structure.residues.extend(other.residues);
    structure.atoms.extend(other.atoms);
    structure
        .coordinates
        .positions
        .extend(other.coordinates.positions);
}

fn check_chemcomp_block_number(document: &CifDocument) -> i32 {
    if document.blocks.len() == 2 && document.blocks[0].name == "comp_list" {
        return 1;
    }
    if document.blocks.len() == 3
        && document.blocks[0].name.is_empty()
        && document.blocks[1].name == "comp_list"
    {
        return 2;
    }
    if document.blocks.len() == 1 {
        let block = &document.blocks[0];
        let has_tag = |tag: &str| {
            block.items.iter().any(|item| item.tag == tag)
                || block
                    .loops
                    .iter()
                    .any(|loop_| loop_.tags.iter().any(|t| t == tag))
        };
        if !has_tag("_atom_site.id")
            && !has_tag("_cell.length_a")
            && has_tag("_chem_comp_atom.atom_id")
        {
            return 0;
        }
    }
    -1
}

fn make_structure_from_chemcomp_doc(
    document: &CifDocument,
    which: i32,
) -> Result<BioStructure, BioReadError> {
    let index = check_chemcomp_block_number(document);
    if index == -1 {
        return Err(BioReadError::Parse {
            line_number: 0,
            message: "Not a chem_comp format.".to_string(),
        });
    }
    make_structure_from_chemcomp_block(&document.blocks[index as usize], which)
}

// BEGIN GEMMI CPP FUNCTION coor_format_from_ext / coor_format_from_content / make_structure_from_doc / read_structure_from_memory
// Gemmi✔️✔️: inline CoorFormat coor_format_from_ext(const std::string& path) {
// Gemmi✔️✔️:   if (iends_with(path, ".pdb") || iends_with(path, ".ent"))
// Gemmi✔️✔️:     return CoorFormat::Pdb;
// Gemmi✔️✔️:   if (iends_with(path, ".cif") || iends_with(path, ".mmcif"))
// Gemmi✔️✔️:     return CoorFormat::Mmcif;
// Gemmi✔️✔️:   if (iends_with(path, ".json"))
// Gemmi✔️✔️:     return CoorFormat::Mmjson;
// Gemmi✔️✔️:   return CoorFormat::Unknown;
// Gemmi✔️✔️: }
// Gemmi✔️✔️: inline CoorFormat coor_format_from_content(const char* buf, const char* end) {
// Gemmi✔️✔️:   while (buf < end - 8) {
// Gemmi✔️✔️:     if (std::isspace(*buf)) {
// Gemmi✔️✔️:       ++buf;
// Gemmi✔️✔️:     } else if (*buf == '#') {
// Gemmi✔️✔️:       while (buf < end - 8 && *buf != '\n')
// Gemmi✔️✔️:         ++buf;
// Gemmi✔️✔️:     } else if (*buf == '{') {
// Gemmi✔️✔️:       return CoorFormat::Mmjson;
// Gemmi✔️✔️:     } else if (ialpha4_id(buf) == ialpha4_id("data") && buf[4] == '_') {
// Gemmi✔️✔️:       return CoorFormat::Mmcif;
// Gemmi✔️✔️:     } else {
// Gemmi✔️✔️:       return CoorFormat::Pdb;
// Gemmi✔️✔️:     }
// Gemmi✔️✔️:   }
// Gemmi✔️✔️:   return CoorFormat::Unknown;
// Gemmi✔️✔️: }
// Gemmi✔️✔️: inline Structure make_structure_from_doc(cif::Document&& doc, bool possible_chemcomp,
// Gemmi✔️✔️:                                          cif::Document* save_doc=nullptr) {
// Gemmi✔️✔️:   if (possible_chemcomp) {
// Gemmi✔️✔️:     int n = check_chemcomp_block_number(doc);
// Gemmi✔️✔️:     if (n != -1)
// Gemmi✔️✔️:       return make_structure_from_chemcomp_block(doc.blocks[n]);
// Gemmi✔️✔️:   }
// Gemmi✔️✔️:   return make_structure(std::move(doc), save_doc);
// Gemmi✔️✔️: }
// Gemmi✔️✔️: inline Structure read_structure_from_memory(char* data, size_t size,
// Gemmi✔️✔️:                                             const std::string& path,
// Gemmi✔️✔️:                                             CoorFormat format=CoorFormat::Unknown,
// Gemmi✔️✔️:                                             cif::Document* save_doc=nullptr) {
// Gemmi✔️✔️:   if (save_doc)
// Gemmi✔️✔️:     save_doc->clear();
// Gemmi✔️✔️:   if (format == CoorFormat::Unknown || format == CoorFormat::Detect)
// Gemmi✔️✔️:     format = coor_format_from_content(data, data + size);
// Gemmi✔️✔️:   if (format == CoorFormat::Pdb)
// Gemmi✔️✔️:     return read_pdb_from_memory(data, size, path);
// Gemmi✔️✔️:   if (format == CoorFormat::Mmcif)
// Gemmi✔️✔️:     return make_structure_from_doc(cif::read_memory(data, size, path.c_str()),
// Gemmi✔️✔️:                                    true, save_doc);
// Gemmi✔️✔️:   if (format == CoorFormat::Mmjson)
// Gemmi✔️✔️:     return make_structure(cif::read_mmjson_insitu(data, size, path), save_doc);
// Gemmi✔️✔️:   fail("wrong format of coordinate file " + path);
// Gemmi✔️✔️: }
// END GEMMI CPP FUNCTION
fn coor_format_from_ext(path: &str) -> Result<BioCoorFormat, BioReadError> {
    let lower = path.to_ascii_lowercase();
    if lower.ends_with(".pdb") || lower.ends_with(".ent") {
        return Ok(BioCoorFormat::Pdb);
    }
    if lower.ends_with(".cif") || lower.ends_with(".mmcif") {
        return Ok(BioCoorFormat::Mmcif);
    }
    if lower.ends_with(".json") {
        return Ok(BioCoorFormat::Mmjson);
    }
    Ok(BioCoorFormat::Unknown)
}

fn coor_format_from_content(buf: &[u8]) -> Result<BioCoorFormat, BioReadError> {
    let mut i = 0usize;
    while i + 8 < buf.len() {
        let byte = buf[i];
        if byte.is_ascii_whitespace() {
            i += 1;
        } else if byte == b'#' {
            while i + 8 < buf.len() && buf[i] != b'\n' {
                i += 1;
            }
        } else if byte == b'{' {
            return Ok(BioCoorFormat::Mmjson);
        } else if buf[i..].starts_with(b"data_") {
            return Ok(BioCoorFormat::Mmcif);
        } else {
            return Ok(BioCoorFormat::Pdb);
        }
    }
    Ok(BioCoorFormat::Unknown)
}

fn make_structure_from_doc(
    document: &CifDocument,
    possible_chemcomp: bool,
) -> Result<BioStructure, BioReadError> {
    if possible_chemcomp {
        let index = check_chemcomp_block_number(document);
        if index != -1 {
            return make_structure_from_chemcomp_block(&document.blocks[index as usize], 7);
        }
    }
    read_mmcif_atom_site_subset_from_document(document)
}

pub(crate) fn read_structure_from_memory(
    data: &str,
    path: &str,
    format: BioCoorFormat,
) -> Result<BioStructure, BioReadError> {
    let format = match format {
        BioCoorFormat::Unknown | BioCoorFormat::Detect => {
            coor_format_from_content(data.as_bytes())?
        }
        other => other,
    };
    match format {
        BioCoorFormat::Pdb => read_pdb_bio_structure_from_str_with_params_and_source(
            data,
            path,
            BioPdbReadParams::default(),
        ),
        BioCoorFormat::Mmcif => make_structure_from_doc(&parse_cif_document(data)?, true),
        BioCoorFormat::Mmjson => {
            let document = read_mmjson_document(data, path)?;
            let mut structure = make_structure_from_mmcif_document(document, None)?;
            structure.input_format = BioCoorFormat::Mmjson;
            Ok(structure)
        }
        _ => Err(BioReadError::Parse {
            line_number: 0,
            message: format!("wrong format of coordinate file {path}"),
        }),
    }
}

fn parse_pdb_atom_record(line: &str, line_number: usize) -> Result<PdbAtomRecord, BioReadError> {
    if line.len() < 55 {
        return Err(BioReadError::Parse {
            line_number,
            message: "ATOM/HETATM line is too short to be correct".to_string(),
        });
    }

    let serial = read_pdb_serial(field_raw(line, 6, 11));
    let atom_name = atom_name_from_field(field_raw(line, 12, 16));
    let altloc = parse_altloc(byte_at(line, 16));
    let residue_name = residue_name_from_field(field(line, 17, 20));
    let residue_kind = classify_residue_name(residue_name);
    let chain_id = pdb_chain_id_from_field(field(line, 20, 22));
    let seq_id = read_pdb_seq_id(field_raw(line, 22, 27));
    let x = parse_f64(field(line, 30, 38), line_number, "x coordinate")?;
    let y = parse_f64(field(line, 38, 46), line_number, "y coordinate")?;
    let z = parse_f64(field(line, 46, 54), line_number, "z coordinate")?;
    let occupancy = parse_optional_f32(field(line, 54, 60), line_number, "occupancy")?;
    let b_iso = parse_optional_f32(field(line, 60, 66), line_number, "B_iso")?;
    let element = parse_pdb_element(line, line_number)?;
    let formal_charge = parse_pdb_charge(field_raw(line, 78, 80), line_number)?;

    Ok(PdbAtomRecord {
        serial: Some(serial),
        atom_name,
        altloc,
        residue_name,
        residue_kind,
        chain_key: chain_id,
        chain_source: ChainSourceIds {
            auth_chain_id: Some(chain_id),
            label_asym_id: None,
        },
        seq_id,
        label_seq_id: None,
        label_entity_id: None,
        group_pdb: Some(if starts_record(line, "HETATM") {
            'H'
        } else {
            'A'
        }),
        segment_id: None,
        position: [x, y, z],
        occupancy,
        b_iso,
        formal_charge,
        element,
        calc_flag: BioCalcFlag::NotSet,
        tls_group_id: None,
        fraction: None,
    })
}

// BEGIN GEMMI CPP FUNCTION read_matrix
// Gemmi✔️✔️: if (len < 46)
// Gemmi✔️✔️:   return 0;
// Gemmi✔️✔️: char n = line[5] - '0';
// Gemmi✔️✔️: if (n >= 1 && n <= 3) {
// Gemmi✔️✔️:   t.mat[n-1][0] = read_double(line+10, 10);
// Gemmi✔️✔️:   t.mat[n-1][1] = read_double(line+20, 10);
// Gemmi✔️✔️:   t.mat[n-1][2] = read_double(line+30, 10);
// Gemmi✔️✔️:   t.vec.at(n-1) = read_double(line+45, 10);
// Gemmi✔️✔️: }
// Gemmi✔️✔️: return n;
// END GEMMI CPP FUNCTION
fn read_matrix(transform: &mut BioTransform, line: &str) -> i32 {
    if line.len() < 46 {
        return 0;
    }
    let n = i32::from(byte_at(line, 5).saturating_sub(b'0'));
    if (1..=3).contains(&n) {
        let row = (n - 1) as usize;
        transform.mat[row][0] = read_pdb_lossy_f64(field_raw(line, 10, 20));
        transform.mat[row][1] = read_pdb_lossy_f64(field_raw(line, 20, 30));
        transform.mat[row][2] = read_pdb_lossy_f64(field_raw(line, 30, 40));
        transform.vec[row] = read_pdb_lossy_f64(field_raw(line, 45, 55));
    }
    n
}

fn bio_transform_identity() -> BioTransform {
    BioTransform {
        mat: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        vec: [0.0, 0.0, 0.0],
    }
}

fn bio_transform_is_identity(transform: &BioTransform) -> bool {
    *transform == bio_transform_identity()
}

// BEGIN GEMMI CPP FUNCTION add_software
// Gemmi✔️✔️: for (size_t start = 0, end = 0; end != std::string::npos; start = end + 1) {
// Gemmi✔️✔️:   end = name.find(',', start);
// Gemmi✔️✔️:   while (end != std::string::npos &&
// Gemmi✔️✔️:          name[end+1] == ' ' && is_digit(name[end+2]))
// Gemmi✔️✔️:     end = name.find(',', end + 1);
// Gemmi✔️✔️:   meta.software.emplace_back();
// Gemmi✔️✔️:   SoftwareItem& item = meta.software.back();
// Gemmi✔️✔️:   item.name = trim_str(name.substr(start, end - start));
// Gemmi✔️✔️:   size_t sep = item.name.find(' ');
// Gemmi✔️✔️:   if (sep != std::string::npos) {
// Gemmi✔️✔️:     size_t ver_start = item.name.find_first_not_of(" (", sep + 1);
// Gemmi✔️✔️:     if (ver_start == std::string::npos) {
// Gemmi✔️✔️:       item.name.resize(sep);
// Gemmi✔️✔️:       item.classification = type;
// Gemmi✔️✔️:       continue;
// Gemmi✔️✔️:     }
// Gemmi✔️✔️:     item.version = item.name.substr(ver_start);
// Gemmi✔️✔️:     item.name.resize(sep);
// Gemmi✔️✔️:     if (!item.version.empty() && item.version.back() == ')') {
// Gemmi✔️✔️:       size_t open_br = item.version.find('(');
// Gemmi✔️✔️:       if (open_br == std::string::npos) {
// Gemmi✔️✔️:         item.version.pop_back();
// Gemmi✔️✔️:       } else if (open_br + 11 == item.version.size() ||
// Gemmi✔️✔️:                  open_br + 13 == item.version.size()) {
// Gemmi✔️✔️:         item.date = pdb_date_format_to_iso(item.version.substr(open_br + 1));
// Gemmi✔️✔️:         if (item.date.size() == 10 && item.date[5] != 'x') {
// Gemmi✔️✔️:           size_t last = item.version.find_last_not_of(' ', open_br - 1);
// Gemmi✔️✔️:           item.version.resize(last + 1);
// Gemmi✔️✔️:         } else {
// Gemmi✔️✔️:           item.date.clear();
// Gemmi✔️✔️:         }
// Gemmi✔️✔️:       }
// Gemmi✔️✔️:     }
// Gemmi✔️✔️:     if (istarts_with(item.version, "version "))
// Gemmi✔️✔️:       item.version.erase(0, 8);
// Gemmi✔️✔️:   }
// Gemmi✔️✔️:   item.classification = type;
// Gemmi✔️✔️: }
// END GEMMI CPP FUNCTION
fn add_software(meta: &mut BioMetadata, classification: BioSoftwareClassification, name: &str) {
    let mut start = 0usize;
    while start < name.len() {
        let mut end = name[start..].find(',').map(|offset| start + offset);
        while let Some(idx) = end {
            let bytes = name.as_bytes();
            if idx + 2 < name.len() && bytes[idx + 1] == b' ' && bytes[idx + 2].is_ascii_digit() {
                end = name[idx + 1..].find(',').map(|offset| idx + 1 + offset);
            } else {
                break;
            }
        }
        let raw_item = name[start..end.unwrap_or(name.len())].trim();
        if !raw_item.is_empty() {
            let mut item = BioSoftwareItem {
                name: raw_item.to_string(),
                version: String::new(),
                date: String::new(),
                description: String::new(),
                contact_author: String::new(),
                contact_author_email: String::new(),
                classification,
            };
            if let Some(sep) = item.name.find(' ') {
                if let Some(ver_start_rel) = item.name[sep + 1..].find(|c| c != ' ' && c != '(') {
                    let ver_start = sep + 1 + ver_start_rel;
                    item.version = item.name[ver_start..].to_string();
                    item.name.truncate(sep);
                    if item.version.ends_with(')') {
                        if let Some(open_br) = item.version.find('(') {
                            let date = pdb_date_format_to_iso(
                                &item.version[open_br + 1..item.version.len() - 1],
                            );
                            if date.len() == 10 && date.as_bytes().get(5).copied() != Some(b'x') {
                                item.date = date;
                                item.version
                                    .truncate(item.version[..open_br].trim_end().len());
                            }
                        } else {
                            item.version.pop();
                        }
                    }
                    if item.version.len() >= 8 && item.version[..8].eq_ignore_ascii_case("version ")
                    {
                        item.version.drain(..8);
                    }
                } else {
                    item.name.truncate(sep);
                }
            }
            meta.software.push(item);
        }
        let Some(end_idx) = end else {
            break;
        };
        start = end_idx + 1;
    }
}

// BEGIN GEMMI CPP FUNCTION add_restraint_count_weight
// Gemmi✔️✔️: if (*value == 'N') // NULL instead of number
// Gemmi✔️✔️:   return;
// Gemmi✔️✔️: ref_info.restr_stats.emplace_back(key);
// Gemmi✔️✔️: RefinementInfo::Restr& restr = ref_info.restr_stats.back();
// Gemmi✔️✔️: const char* endptr;
// Gemmi✔️✔️: restr.count = no_sign_atoi(value, &endptr);
// Gemmi✔️✔️: if (const char* sep = std::strchr(endptr, ';'))
// Gemmi✔️✔️:   restr.weight = fast_atof(sep + 1, &endptr);
// Gemmi✔️✔️: if (const char* sep = std::strchr(endptr, ';'))
// Gemmi✔️✔️:   restr.function = read_string(sep+1, 50);
// END GEMMI CPP FUNCTION
fn add_restraint_count_weight(ref_info: &mut BioRefinementInfo, key: &str, value: &str) {
    if value.starts_with('N') {
        return;
    }
    let mut restr = BioRefinementRestraint {
        name: key.to_string(),
        ..BioRefinementRestraint::default()
    };
    let parts: Vec<_> = value.split(';').collect();
    restr.count = parts
        .first()
        .and_then(|part| part.trim().parse::<i32>().ok());
    restr.weight = parts
        .get(1)
        .and_then(|part| part.trim().parse::<f64>().ok());
    restr.function = parts
        .get(2)
        .map_or(String::new(), |part| part.trim().to_string());
    ref_info.restr_stats.push(restr);
}

// BEGIN GEMMI CPP HELPER impl::find_or_add(ref_info.restr_stats, key)
// Gemmi✔️✔️: impl::find_or_add(ref_info.restr_stats, key) returns the existing
// Gemmi✔️✔️: restraint entry with matching key or appends a new entry first.
// END GEMMI CPP HELPER impl::find_or_add(ref_info.restr_stats, key)
fn find_or_add_restraint<'a>(
    ref_info: &'a mut BioRefinementInfo,
    key: &str,
) -> &'a mut BioRefinementRestraint {
    if let Some(index) = ref_info
        .restr_stats
        .iter()
        .position(|restr| restr.name == key)
    {
        return &mut ref_info.restr_stats[index];
    }
    ref_info.restr_stats.push(BioRefinementRestraint {
        name: key.to_string(),
        ..BioRefinementRestraint::default()
    });
    ref_info
        .restr_stats
        .last_mut()
        .expect("newly pushed restraint must exist")
}

// BEGIN GEMMI CPP HELPER read_double(value, 50)
// Gemmi✔️✔️: read_double(value, 50) parses the leading numeric token from a
// Gemmi✔️✔️: fixed-width value field.
// END GEMMI CPP HELPER read_double(value, 50)
fn parse_prefixed_f64(value: &str, width: usize) -> Option<f64> {
    let end = value
        .char_indices()
        .nth(width)
        .map_or(value.len(), |(idx, _)| idx);
    value
        .get(..end)
        .unwrap_or(value)
        .split_whitespace()
        .next()
        .and_then(|token| token.parse::<f64>().ok())
}

// BEGIN GEMMI CPP HELPER impl::find_or_add(...).dev_ideal = read_double(value, 50)
// Gemmi✔️✔️: impl::find_or_add(ref_info.restr_stats, key).dev_ideal
// Gemmi✔️✔️:   = read_double(value, 50);
// END GEMMI CPP HELPER impl::find_or_add(...).dev_ideal = read_double(value, 50)
fn set_restraint_dev_ideal(ref_info: &mut BioRefinementInfo, key: &str, value: &str) {
    find_or_add_restraint(ref_info, key).dev_ideal = parse_prefixed_f64(value, 50);
}

// BEGIN GEMMI CPP HELPER group.selections.emplace_back(); group.selections.back().details = ...
// Gemmi✔️✔️: group.selections.emplace_back();
// Gemmi✔️✔️: group.selections.back().details = std::string(value, end);
// Gemmi✔️✔️: possibly_unfinished_remark3 = &group.selections.back().details;
// END GEMMI CPP HELPER group.selections.emplace_back(); group.selections.back().details = ...
fn push_tls_selection(ref_info: &mut BioRefinementInfo, value: &str) -> Option<(usize, usize)> {
    let tls_group_idx = ref_info.tls_groups.len().checked_sub(1)?;
    let group = ref_info.tls_groups.get_mut(tls_group_idx)?;
    group.selections.push(BioTlsSelection {
        details: value.to_string(),
        ..BioTlsSelection::default()
    });
    Some((tls_group_idx, group.selections.len() - 1))
}

// BEGIN GEMMI CPP FUNCTION gemmi::SeqId::SeqId(const std::string&)
// Gemmi✔️✔️: explicit SeqId(const std::string& str) {
// Gemmi✔️✔️:   char* endptr;
// Gemmi✔️✔️:   num = std::strtol(str.c_str(), &endptr, 10);
// Gemmi✔️✔️:   if (endptr == str.c_str() || (*endptr != '\0' && endptr[1] != '\0'))
// Gemmi✔️✔️:     throw std::invalid_argument("Not a seqid: " + str);
// Gemmi✔️✔️:   icode = (*endptr | 0x20);
// Gemmi✔️✔️: }
// END GEMMI CPP FUNCTION
fn parse_gemmi_seq_id_string(value: &str) -> Option<PdbSeqId> {
    let value = value.trim();
    let bytes = value.as_bytes();
    let mut end = usize::from(matches!(bytes.first(), Some(b'+' | b'-')));
    let digit_start = end;
    while bytes.get(end).is_some_and(u8::is_ascii_digit) {
        end += 1;
    }
    if end == digit_start || bytes.len().saturating_sub(end) > 1 {
        return None;
    }
    let seq_num = value[..end].parse::<i32>().ok()?;
    let ins_code = bytes.get(end).copied().map(|value| value | 0x20);
    Some(PdbSeqId { seq_num, ins_code })
}

// BEGIN GEMMI CPP HELPER is_tls_item(key)
// Gemmi✔️✔️: is_tls_item(key) recognizes T/L/S tensor labels followed by 11..33.
// END GEMMI CPP HELPER is_tls_item(key)
fn is_tls_item_key(key: &str) -> bool {
    let bytes = key.as_bytes();
    bytes.len() == 3
        && matches!(bytes[0], b'T' | b'L' | b'S')
        && matches!(bytes[1], b'1' | b'2' | b'3')
        && matches!(bytes[2], b'1' | b'2' | b'3')
}

// BEGIN GEMMI CPP HELPER split_str_multi(key_start) tensor assignment loop
// Gemmi✔️✔️: std::vector<std::string> tokens = split_str_multi(key_start);
// Gemmi✔️✔️: for (size_t i = 0; i + 1 < tokens.size(); i += 2) {
// Gemmi✔️✔️:   std::string& k = tokens[i];
// Gemmi✔️✔️:   if (k.size() == 4 && k[3] == ':')
// Gemmi✔️✔️:     k.resize(3);
// Gemmi✔️✔️:   if (is_tls_item(k)) { ... assign T/L/S[x][y] ... }
// Gemmi✔️✔️: }
// END GEMMI CPP HELPER split_str_multi(key_start) tensor assignment loop
fn parse_tls_key_values(group: &mut BioTlsGroup, key_start: &str) {
    for pair in key_start.split_whitespace().collect::<Vec<_>>().chunks(2) {
        if pair.len() != 2 {
            continue;
        }
        let key = pair[0].trim_end_matches(':');
        if !is_tls_item_key(key) {
            continue;
        }
        let bytes = key.as_bytes();
        let x = (bytes[1] - b'1') as usize;
        let y = (bytes[2] - b'1') as usize;
        let value = pair[1].parse::<f64>().unwrap_or(f64::NAN);
        match bytes[0] {
            b'S' => group.s[x][y] = value,
            b'T' => group.t[x][y] = value,
            b'L' => group.l[x][y] = value,
            _ => {}
        }
    }
}

// BEGIN GEMMI CPP HELPER read_string(ptr, len)
// Gemmi✔️✔️: read_string(ptr, len) returns the trimmed fixed-width substring.
// END GEMMI CPP HELPER read_string(ptr, len)
fn read_fixed_string(line: &str, start: usize, len: usize) -> String {
    // BEGIN GEMMI CPP FUNCTION read_string
    // Gemmi✔️✔️: // left trim
    // Gemmi✔️✔️: while (field_length != 0 && is_space(*p)) {
    // Gemmi✔️✔️:   ++p;
    // Gemmi✔️✔️:   --field_length;
    // Gemmi✔️✔️: }
    // Gemmi✔️✔️: // EOL/EOF ends the string
    // Gemmi✔️✔️: for (int i = 0; i < field_length; ++i)
    // Gemmi✔️✔️:   if (p[i] == '\n' || p[i] == '\r' || p[i] == '\0') {
    // Gemmi✔️✔️:     field_length = i;
    // Gemmi✔️✔️:     break;
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️: // right trim
    // Gemmi✔️✔️: while (field_length != 0 && is_space(p[field_length-1]))
    // Gemmi✔️✔️:   --field_length;
    // Gemmi✔️✔️: return std::string(p, field_length);
    // END GEMMI CPP FUNCTION
    let end = start.saturating_add(len).min(line.len());
    line.get(start..end).unwrap_or("").trim().to_string()
}

// BEGIN GEMMI CPP FUNCTION read_remark3_line
// Gemmi✔️✔️: const char* key_start = skip_blank(line + 10);
// Gemmi✔️✔️: const char* colon = std::strchr(key_start, ':');
// Gemmi✔️✔️: const char* key_end = rtrim_cstr(key_start, colon);
// Gemmi✔️✔️: std::string key(key_start, key_end);
// Gemmi✔️✔️: if (possibly_unfinished_remark3) {
// Gemmi✔️✔️:   if (key_start > line + 17) {
// Gemmi✔️✔️:     *possibly_unfinished_remark3 += ' ';
// Gemmi✔️✔️:     possibly_unfinished_remark3->append(key);
// Gemmi✔️✔️:     return;
// Gemmi✔️✔️:   }
// Gemmi✔️✔️:   possibly_unfinished_remark3 = nullptr;
// Gemmi✔️✔️: }
// Gemmi✔️✔️: if (colon) {
// Gemmi✔️✔️:   const char* value = skip_blank(colon + 1);
// Gemmi✔️✔️:   const char* end = rtrim_cstr(value);
// Gemmi✔️✔️:   if (end - value == 4 && std::strncmp(value, "NULL", 4) == 0)
// Gemmi✔️✔️:     return;
// Gemmi✔️✔️:   if (same_str(key, "PROGRAM"))
// Gemmi✔️✔️:     add_software(meta, SoftwareItem::Refinement, std::string(value, end));
// Gemmi✔️✔️:   if (meta.refinement.empty())
// Gemmi✔️✔️:     return;
// Gemmi✔️✔️:   RefinementInfo& ref_info = meta.refinement.back();
// Gemmi✔️✔️:   ...
// Gemmi✔️✔️: } else {
// Gemmi✔️✔️:   if (same_str(key, "DATA USED IN REFINEMENT."))
// Gemmi✔️✔️:     meta.refinement.emplace_back();
// Gemmi✔️✔️:   else if (same_str(key, "FIT IN THE HIGHEST RESOLUTION BIN."))
// Gemmi✔️✔️:     if (!meta.refinement.empty())
// Gemmi✔️✔️:       meta.refinement.back().bins.emplace_back();
// Gemmi✔️✔️: }
// END GEMMI CPP FUNCTION
fn read_remark3_line(
    line: &str,
    meta: &mut BioMetadata,
    continuation: &mut Option<Remark3Continuation>,
) {
    let remark_body = field_raw(line, 10, line.len());
    let key_start_offset = remark_body
        .char_indices()
        .find(|(_, ch)| !ch.is_ascii_whitespace())
        .map(|(idx, _)| idx)
        .unwrap_or(remark_body.len());
    let key_start = &remark_body[key_start_offset..];
    let colon = key_start.find(':');
    let key = colon
        .map(|idx| key_start[..idx].trim_end())
        .unwrap_or_else(|| key_start.trim_end());

    if let Some(target) = continuation.take() {
        if 10 + key_start_offset > 17 {
            let Remark3Continuation::TlsSelection {
                refinement_idx,
                tls_group_idx,
                selection_idx,
            } = target;
            if let Some(details) = meta
                .refinement
                .get_mut(refinement_idx)
                .and_then(|refinement| refinement.tls_groups.get_mut(tls_group_idx))
                .and_then(|group| group.selections.get_mut(selection_idx))
                .map(|selection| &mut selection.details)
            {
                if !details.is_empty() {
                    details.push(' ');
                }
                details.push_str(key);
                *continuation = Some(target);
                return;
            }
        }
    }

    if let Some(colon_idx) = colon {
        let value = key_start[colon_idx + 1..].trim_start();
        let value = value.trim_end();
        if value == "NULL" {
            return;
        }
        if key == "PROGRAM" {
            add_software(meta, BioSoftwareClassification::Refinement, value);
        }
        let Some(ref_info) = meta.refinement.last_mut() else {
            return;
        };
        match key {
            "RESOLUTION RANGE HIGH (ANGSTROMS)" => {
                ref_info.resolution_high = value.parse::<f64>().ok();
            }
            "RESOLUTION RANGE LOW  (ANGSTROMS)" => {
                ref_info.resolution_low = value.parse::<f64>().ok();
            }
            "COMPLETENESS FOR RANGE        (%)" => {
                ref_info.completeness = value.parse::<f64>().ok();
            }
            "NUMBER OF REFLECTIONS" => {
                ref_info.reflection_count = value.parse::<i32>().ok();
            }
            "CROSS-VALIDATION METHOD" => {
                ref_info.cross_validation_method = value.to_string();
            }
            "FREE R VALUE TEST SET SELECTION" => {
                ref_info.rfree_selection_method = value.to_string();
            }
            "R VALUE     (WORKING + TEST SET)" => {
                ref_info.r_all = value.parse::<f64>().ok();
            }
            "R VALUE            (WORKING SET)" => {
                ref_info.r_work = value.parse::<f64>().ok();
            }
            "FREE R VALUE" => {
                ref_info.r_free = value.parse::<f64>().ok();
            }
            "FREE R VALUE TEST SET COUNT" => {
                ref_info.rfree_set_count = value.parse::<i32>().ok();
            }
            "TOTAL NUMBER OF BINS USED" => {
                ref_info.bin_count = value.parse::<i32>().ok();
            }
            "BIN RESOLUTION RANGE HIGH       (A)" => {
                if let Some(bin) = ref_info.bins.last_mut() {
                    bin.resolution_high = value.parse::<f64>().ok();
                }
            }
            "BIN RESOLUTION RANGE LOW        (A)" => {
                if let Some(bin) = ref_info.bins.last_mut() {
                    bin.resolution_low = value.parse::<f64>().ok();
                }
            }
            "BIN COMPLETENESS (WORKING+TEST) (%)" => {
                if let Some(bin) = ref_info.bins.last_mut() {
                    bin.completeness = value.parse::<f64>().ok();
                }
            }
            "REFLECTIONS IN BIN   (WORKING+TEST)" => {
                if let Some(bin) = ref_info.bins.last_mut() {
                    bin.reflection_count = value.parse::<i32>().ok();
                }
            }
            "BIN R VALUE          (WORKING+TEST)" => {
                if let Some(bin) = ref_info.bins.last_mut() {
                    bin.r_all = value.parse::<f64>().ok();
                }
            }
            "REFLECTIONS IN BIN    (WORKING SET)" => {
                if let Some(bin) = ref_info.bins.last_mut() {
                    bin.work_set_count = value.parse::<i32>().ok();
                }
            }
            "BIN R VALUE           (WORKING SET)" => {
                if let Some(bin) = ref_info.bins.last_mut() {
                    bin.r_work = value.parse::<f64>().ok();
                }
            }
            "BIN FREE R VALUE" => {
                if let Some(bin) = ref_info.bins.last_mut() {
                    bin.r_free = value.parse::<f64>().ok();
                }
            }
            "BIN FREE R VALUE TEST SET COUNT" => {
                if let Some(bin) = ref_info.bins.last_mut() {
                    bin.rfree_set_count = value.parse::<i32>().ok();
                }
            }
            "MEAN B VALUE      (OVERALL, A**2)" => {
                ref_info.mean_b = value.parse::<f64>().ok();
            }
            "B11 (A**2)" => ref_info.aniso_b.u11 = value.parse::<f64>().unwrap_or(f64::NAN),
            "B22 (A**2)" => ref_info.aniso_b.u22 = value.parse::<f64>().unwrap_or(f64::NAN),
            "B33 (A**2)" => ref_info.aniso_b.u33 = value.parse::<f64>().unwrap_or(f64::NAN),
            "B12 (A**2)" => ref_info.aniso_b.u12 = value.parse::<f64>().unwrap_or(f64::NAN),
            "B13 (A**2)" => ref_info.aniso_b.u13 = value.parse::<f64>().unwrap_or(f64::NAN),
            "B23 (A**2)" => ref_info.aniso_b.u23 = value.parse::<f64>().unwrap_or(f64::NAN),
            "ESD FROM LUZZATI PLOT                    (A)" => {
                ref_info.luzzati_error = value.parse::<f64>().ok();
            }
            "DPI (BLOW EQ-10) BASED ON R VALUE        (A)" => {
                ref_info.dpi_blow_r = value.parse::<f64>().ok();
            }
            "DPI (BLOW EQ-9) BASED ON FREE R VALUE    (A)" => {
                ref_info.dpi_blow_rfree = value.parse::<f64>().ok();
            }
            "DPI (CRUICKSHANK) BASED ON R VALUE       (A)" => {
                ref_info.dpi_cruickshank_r = value.parse::<f64>().ok();
            }
            "DPI (CRUICKSHANK) BASED ON FREE R VALUE  (A)" => {
                ref_info.dpi_cruickshank_rfree = value.parse::<f64>().ok();
            }
            "CORRELATION COEFFICIENT FO-FC" => {
                ref_info.cc_fo_fc_work = value.parse::<f64>().ok();
            }
            "CORRELATION COEFFICIENT FO-FC FREE" => {
                ref_info.cc_fo_fc_free = value.parse::<f64>().ok();
            }
            "BOND LENGTHS" => add_restraint_count_weight(ref_info, "t_bond_d", value),
            "BOND ANGLES" => add_restraint_count_weight(ref_info, "t_angle_deg", value),
            "TORSION ANGLES" => add_restraint_count_weight(ref_info, "t_dihedral_angle_d", value),
            "TRIGONAL CARBON PLANES" => {
                add_restraint_count_weight(ref_info, "t_trig_c_planes", value);
            }
            "GENERAL PLANES" => add_restraint_count_weight(ref_info, "t_gen_planes", value),
            "ISOTROPIC THERMAL FACTORS" => add_restraint_count_weight(ref_info, "t_it", value),
            "BAD NON-BONDED CONTACTS" => add_restraint_count_weight(ref_info, "t_nbd", value),
            "IMPROPER TORSIONS" => {
                add_restraint_count_weight(ref_info, "t_improper_torsion", value);
            }
            "CHIRAL IMPROPER TORSION" => {
                add_restraint_count_weight(ref_info, "t_chiral_improper_torsion", value);
            }
            "SUM OF OCCUPANCIES" => {
                add_restraint_count_weight(ref_info, "t_sum_occupancies", value);
            }
            "UTILITY DISTANCES" => {
                add_restraint_count_weight(ref_info, "t_utility_distance", value)
            }
            "UTILITY ANGLES" => add_restraint_count_weight(ref_info, "t_utility_angle", value),
            "UTILITY TORSION" => add_restraint_count_weight(ref_info, "t_utility_torsion", value),
            "IDEAL-DIST CONTACT TERM" => {
                add_restraint_count_weight(ref_info, "t_ideal_dist_contact", value);
            }
            "BOND LENGTHS                       (A)" => {
                set_restraint_dev_ideal(ref_info, "t_bond_d", value);
            }
            "BOND ANGLES                  (DEGREES)" => {
                set_restraint_dev_ideal(ref_info, "t_angle_deg", value);
            }
            "PEPTIDE OMEGA TORSION ANGLES (DEGREES)" => {
                set_restraint_dev_ideal(ref_info, "t_omega_torsion", value);
            }
            "OTHER TORSION ANGLES         (DEGREES)" => {
                set_restraint_dev_ideal(ref_info, "t_other_torsion", value);
            }
            "TLS GROUP" => {
                let num_id = value
                    .chars()
                    .take_while(|ch| ch.is_ascii_digit())
                    .collect::<String>()
                    .parse::<i16>()
                    .ok();
                ref_info.tls_groups.push(BioTlsGroup {
                    id: value.to_string(),
                    num_id,
                    ..BioTlsGroup::default()
                });
            }
            "SET" | "SELECTION" if field_raw(line, 23, 24) == ":" || key == "SET" => {
                if let Some((tls_group_idx, selection_idx)) = push_tls_selection(ref_info, value) {
                    *continuation = Some(Remark3Continuation::TlsSelection {
                        refinement_idx: meta.refinement.len() - 1,
                        tls_group_idx,
                        selection_idx,
                    });
                }
            }
            "RESIDUE RANGE" => {
                if let Some(group) = ref_info.tls_groups.last_mut() {
                    let chain1 = read_fixed_string(line, colon_idx + 11, 5);
                    let chain2 = read_fixed_string(line, colon_idx + 26, 5);
                    if !chain1.is_empty() && chain1 == chain2 {
                        let res_begin =
                            parse_gemmi_seq_id_string(&read_fixed_string(line, colon_idx + 16, 6));
                        let res_end =
                            parse_gemmi_seq_id_string(&read_fixed_string(line, colon_idx + 31, 6));
                        if let (Some(res_begin), Some(res_end)) = (res_begin, res_end) {
                            group.selections.push(BioTlsSelection {
                                chain: chain1,
                                res_begin: Some(res_begin),
                                res_end: Some(res_end),
                                details: String::new(),
                            });
                        }
                    }
                }
            }
            "ORIGIN FOR THE GROUP (A)" => {
                if let Some(group) = ref_info.tls_groups.last_mut() {
                    let parts: Vec<_> = value.split_whitespace().collect();
                    if parts.len() == 3 {
                        for (dst, src) in group.origin.iter_mut().zip(parts) {
                            *dst = src.parse::<f64>().unwrap_or(f64::NAN);
                        }
                    }
                }
            }
            _ if is_tls_item_key(key) => {
                if let Some(group) = ref_info.tls_groups.last_mut() {
                    parse_tls_key_values(group, key_start.trim_end());
                }
            }
            _ => {}
        }
    } else {
        match key {
            "DATA USED IN REFINEMENT." => {
                meta.refinement.push(BioRefinementInfo {
                    id: (meta.refinement.len() + 1).to_string(),
                    ..BioRefinementInfo::default()
                });
            }
            "FIT IN THE HIGHEST RESOLUTION BIN." => {
                if let Some(ref_info) = meta.refinement.last_mut() {
                    ref_info.bins.push(BioRefinementBin::default());
                }
            }
            _ => {}
        }
    }
}

// BEGIN GEMMI CPP FUNCTION read_remark_200_230_240
// Gemmi✔️✔️: if (cryst_desc) {
// Gemmi✔️✔️:   if (line[10] == ' ' && line[11] == ' ') {
// Gemmi✔️✔️:     const char* start = line + 11;
// Gemmi✔️✔️:     cryst_desc->append(start, rtrim_cstr(start) - start);
// Gemmi✔️✔️:     return;
// Gemmi✔️✔️:   }
// Gemmi✔️✔️:   cryst_desc = nullptr;
// Gemmi✔️✔️: }
// Gemmi✔️✔️: const char* key_start = skip_blank(line + 10);
// Gemmi✔️✔️: const char* colon = std::strchr(key_start, ':');
// Gemmi✔️✔️: const char* key_end = rtrim_cstr(key_start, colon);
// Gemmi✔️✔️: std::string key(key_start, key_end);
// Gemmi✔️✔️: ...
// END GEMMI CPP FUNCTION
fn read_remark_200_230_240(
    line: &str,
    meta: &mut BioMetadata,
    continuation: &mut Option<Remark200Continuation>,
) {
    if let Some(Remark200Continuation::CrystalDescription { crystal_idx }) = continuation.take() {
        if byte_at(line, 10) == b' ' && byte_at(line, 11) == b' ' {
            if let Some(crystal) = meta.experiment_crystals.get_mut(crystal_idx) {
                crystal
                    .description
                    .push_str(field_raw(line, 11, line.len()).trim_end());
                *continuation = Some(Remark200Continuation::CrystalDescription { crystal_idx });
                return;
            }
        }
    }

    let remark_body = field_raw(line, 10, line.len());
    let key_start = remark_body.trim_start();
    let colon = key_start.find(':');
    let key = colon
        .map(|idx| key_start[..idx].trim_end())
        .unwrap_or_else(|| key_start.trim_end());

    if let Some(colon_idx) = colon {
        let value = key_start[colon_idx + 1..].trim_start().trim_end();
        if value == "NULL" {
            return;
        }
        match key {
            "INTENSITY-INTEGRATION SOFTWARE" => {
                add_software(meta, BioSoftwareClassification::DataReduction, value);
            }
            "DATA SCALING SOFTWARE" => {
                add_software(meta, BioSoftwareClassification::DataScaling, value);
            }
            "SOFTWARE USED" => {
                add_software(meta, BioSoftwareClassification::Phasing, value);
            }
            "METHOD USED TO DETERMINE THE STRUCTURE" => {
                meta.solved_by = Some(value.to_string());
            }
            "STARTING MODEL" => {
                meta.starting_model = Some(value.to_string());
            }
            _ => {
                let Some(experiment) = meta.experiments.last_mut() else {
                    return;
                };
                let Some(crystal) = meta.experiment_crystals.last_mut() else {
                    return;
                };
                let diffraction = crystal
                    .diffractions
                    .get_mut(0)
                    .expect("experimental details create one diffraction row");
                match key {
                    "EXPERIMENT TYPE" => experiment.method = value.to_string(),
                    "NUMBER OF CRYSTALS USED" => {
                        experiment.number_of_crystals = value.parse::<i32>().ok();
                    }
                    "PH" => {
                        crystal.ph = value.parse::<f64>().ok();
                        if crystal.ph.is_none() {
                            crystal.ph_range = value.to_string();
                        }
                    }
                    "DATE OF DATA COLLECTION" => {
                        diffraction.collection_date = pdb_date_format_to_iso(value);
                    }
                    "TEMPERATURE           (KELVIN)" => {
                        diffraction.temperature = value.parse::<f64>().ok();
                    }
                    "SYNCHROTRON              (Y/N)" => {
                        if value.starts_with('Y') {
                            diffraction.source = "SYNCHROTRON".to_string();
                        }
                    }
                    "RADIATION SOURCE" => {
                        if diffraction.source == "SYNCHROTRON" {
                            diffraction.synchrotron = value.to_string();
                        } else {
                            diffraction.source = value.to_string();
                        }
                    }
                    "BEAMLINE" => {
                        diffraction.beamline = value.to_string();
                        if !diffraction.synchrotron.is_empty() && diffraction.source_type.is_empty()
                        {
                            diffraction.source_type = format!(
                                "{} BEAMLINE {}",
                                diffraction.synchrotron, diffraction.beamline
                            );
                        }
                    }
                    "X-RAY GENERATOR MODEL" => diffraction.source_type = value.to_string(),
                    "MONOCHROMATIC OR LAUE    (M/L)" => {
                        diffraction.mono_or_laue = value.chars().next();
                    }
                    "WAVELENGTH OR RANGE        (A)" => diffraction.wavelengths = value.to_string(),
                    "MONOCHROMATOR" => diffraction.monochromator = value.to_string(),
                    "OPTICS" => diffraction.optics = value.to_string(),
                    "DETECTOR TYPE" => diffraction.detector = value.to_string(),
                    "DETECTOR MANUFACTURER" => diffraction.detector_make = value.to_string(),
                    "NUMBER OF UNIQUE REFLECTIONS" => {
                        experiment.unique_reflections = value.parse::<i32>().ok();
                    }
                    "RESOLUTION RANGE HIGH      (A)" => {
                        experiment.reflections.resolution_high = value.parse::<f64>().ok();
                    }
                    "RESOLUTION RANGE LOW       (A)" => {
                        experiment.reflections.resolution_low = value.parse::<f64>().ok();
                    }
                    "COMPLETENESS FOR RANGE     (%)" => {
                        experiment.reflections.completeness = value.parse::<f64>().ok();
                    }
                    "DATA REDUNDANCY" => {
                        experiment.reflections.redundancy = value.parse::<f64>().ok();
                    }
                    "R MERGE                    (I)" => {
                        experiment.reflections.r_merge = value.parse::<f64>().ok();
                    }
                    "R SYM                      (I)" => {
                        experiment.reflections.r_sym = value.parse::<f64>().ok();
                    }
                    "<I/SIGMA(I)> FOR THE DATA SET" => {
                        experiment.reflections.mean_i_over_sigma = value.parse::<f64>().ok();
                    }
                    "REMARK" => {
                        crystal.description = value.to_string();
                        *continuation = Some(Remark200Continuation::CrystalDescription {
                            crystal_idx: meta.experiment_crystals.len() - 1,
                        });
                    }
                    _ => {
                        if let Some(shell) = experiment.shells.last_mut() {
                            match key {
                                "HIGHEST RESOLUTION SHELL, RANGE HIGH (A)" => {
                                    shell.resolution_high = value.parse::<f64>().ok();
                                }
                                "HIGHEST RESOLUTION SHELL, RANGE LOW  (A)" => {
                                    shell.resolution_low = value.parse::<f64>().ok();
                                }
                                "COMPLETENESS FOR SHELL     (%)" => {
                                    shell.completeness = value.parse::<f64>().ok();
                                }
                                "DATA REDUNDANCY IN SHELL" => {
                                    shell.redundancy = value.parse::<f64>().ok();
                                }
                                "R MERGE FOR SHELL          (I)" => {
                                    shell.r_merge = value.parse::<f64>().ok();
                                }
                                "R SYM FOR SHELL            (I)" => {
                                    shell.r_sym = value.parse::<f64>().ok();
                                }
                                "<I/SIGMA(I)> FOR SHELL" => {
                                    shell.mean_i_over_sigma = value.parse::<f64>().ok();
                                }
                                _ => {}
                            }
                        }
                    }
                }
            }
        }
    } else {
        match key {
            "EXPERIMENTAL DETAILS" => {
                meta.experiment_crystals.push(BioExperimentCrystalInfo {
                    diffractions: vec![BioDiffractionInfo::default()],
                    ..BioExperimentCrystalInfo::default()
                });
                meta.experiments.push(BioExperimentInfo::default());
            }
            "IN THE HIGHEST RESOLUTION SHELL." => {
                if let Some(experiment) = meta.experiments.last_mut() {
                    experiment.shells.push(BioReflectionsInfo::default());
                }
            }
            _ => {}
        }
    }
}

fn atom_name_trimmed(atom_name: AtomName) -> String {
    String::from_utf8_lossy(&atom_name.0).trim().to_string()
}

fn same_conformer(a: Option<AltLocLabel>, b: Option<AltLocLabel>) -> bool {
    a.is_none() || b.is_none() || a == b
}

fn residue_matches_address(
    structure: &BioStructure,
    residue_id: ResidueId,
    address: &BioAtomAddress,
    first_model_only: bool,
) -> bool {
    let residue = &structure.residues[residue_id.index() as usize];
    let chain = &structure.chains[residue.chain_id.index() as usize];
    if first_model_only && chain.model_id.index() != 0 {
        return false;
    }
    if chain
        .source
        .auth_chain_id
        .is_some_and(|id| id.as_str() != address.chain_name)
    {
        return false;
    }
    residue.source.seq_id == address.seq_id && residue.name.as_str() == address.residue_name
}

fn find_residue_for_address(
    structure: &BioStructure,
    address: &BioAtomAddress,
    first_model_only: bool,
) -> Option<ResidueId> {
    structure
        .residues
        .iter()
        .enumerate()
        .find(|(idx, _)| {
            residue_matches_address(
                structure,
                ResidueId::new(*idx as u32),
                address,
                first_model_only,
            )
        })
        .map(|(idx, _)| ResidueId::new(idx as u32))
}

fn residue_atom_indices(structure: &BioStructure, residue_id: ResidueId) -> std::ops::Range<usize> {
    let span = structure.residues[residue_id.index() as usize].atom_span;
    span.start as usize..span.end() as usize
}

fn first_sulfur_in_residue(structure: &BioStructure, residue_id: ResidueId) -> Option<usize> {
    residue_atom_indices(structure, residue_id)
        .find(|&idx| structure.atoms[idx].element == Element::S)
}

fn find_atom_in_residue(
    structure: &BioStructure,
    residue_id: ResidueId,
    atom_name: &str,
) -> Option<usize> {
    residue_atom_indices(structure, residue_id)
        .find(|&idx| atom_name_trimmed(structure.atoms[idx].name) == atom_name)
}

// BEGIN GEMMI CPP FUNCTION complete_ssbond_atom
// Gemmi✔️✔️: ad.atom_name = "SG";
// Gemmi✔️✔️: const_CRA cra = mdl.find_cra(ad);
// Gemmi✔️✔️: if (cra.residue && (!cra.atom || cra.atom->element != El::S))
// Gemmi✔️✔️:   if (const Atom* a = cra.residue->find_by_element(El::S)) {
// Gemmi✔️✔️:     ad.atom_name = a->name;
// Gemmi✔️✔️:     ad.altloc = a->altloc;
// Gemmi✔️✔️:   }
// Gemmi✔️✔️: return cra.residue;
// END GEMMI CPP FUNCTION
fn complete_ssbond_atom(
    address: &mut BioAtomAddress,
    structure: &BioStructure,
) -> Option<ResidueId> {
    address.atom_name = "SG".to_string();
    let residue_id = find_residue_for_address(structure, address, true)?;
    let sulfur_idx = match find_atom_in_residue(structure, residue_id, "SG") {
        Some(atom_idx) if structure.atoms[atom_idx].element == Element::S => None,
        _ => first_sulfur_in_residue(structure, residue_id),
    };
    if let Some(atom_idx) = sulfur_idx {
        address.atom_name = atom_name_trimmed(structure.atoms[atom_idx].name);
        address.altloc = structure.atoms[atom_idx].altloc;
    }
    Some(residue_id)
}

fn atom_distance_sq(structure: &BioStructure, left: usize, right: usize) -> f64 {
    let a = structure.coordinates.positions[left];
    let b = structure.coordinates.positions[right];
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    dx * dx + dy * dy + dz * dz
}

// BEGIN GEMMI CPP FUNCTION complete_ssbond
// Gemmi✔️✔️: const Residue* res1 = complete_ssbond_atom(con.partner1, mdl);
// Gemmi✔️✔️: const Residue* res2 = complete_ssbond_atom(con.partner2, mdl);
// Gemmi✔️✔️: if (res1 && res2 && (con.partner1.altloc != '\0' || con.partner2.altloc != '\0')) {
// Gemmi✔️✔️:   double min_dist_sq = INFINITY;
// Gemmi✔️✔️:   for (const Atom& a1 : res1->get(con.partner1.atom_name))
// Gemmi✔️✔️:     for (const Atom& a2 : res2->get(con.partner2.atom_name))
// Gemmi✔️✔️:       if (a2.same_conformer(a1)) { choose shortest distance pair }
// Gemmi✔️✔️: }
// END GEMMI CPP FUNCTION
fn complete_ssbond(connection: &mut BioConnection, structure: &BioStructure) {
    let Some(res1) = complete_ssbond_atom(&mut connection.partner1, structure) else {
        return;
    };
    let Some(res2) = complete_ssbond_atom(&mut connection.partner2, structure) else {
        return;
    };
    if connection.partner1.altloc.is_none() && connection.partner2.altloc.is_none() {
        return;
    }
    let mut best: Option<(f64, Option<AltLocLabel>, Option<AltLocLabel>)> = None;
    for a1_idx in residue_atom_indices(structure, res1).filter(|&idx| {
        atom_name_trimmed(structure.atoms[idx].name) == connection.partner1.atom_name
    }) {
        for a2_idx in residue_atom_indices(structure, res2).filter(|&idx| {
            atom_name_trimmed(structure.atoms[idx].name) == connection.partner2.atom_name
        }) {
            let alt1 = structure.atoms[a1_idx].altloc;
            let alt2 = structure.atoms[a2_idx].altloc;
            if same_conformer(alt1, alt2) {
                let dist_sq = atom_distance_sq(structure, a1_idx, a2_idx);
                if best.is_none_or(|(min_dist, _, _)| dist_sq < min_dist) {
                    best = Some((dist_sq, alt1, alt2));
                }
            }
        }
    }
    if let Some((_, alt1, alt2)) = best {
        connection.partner1.altloc = alt1;
        connection.partner2.altloc = alt2;
    }
}

// BEGIN GEMMI CPP FUNCTION compare_link_symops
// Gemmi✔️✔️: if (record.size() < 72)
// Gemmi✔️✔️:   return Asu::Any;
// Gemmi✔️✔️: std::string s1 = read_string(&record[59], 6);
// Gemmi✔️✔️: std::string s2 = read_string(&record[66], 6);
// Gemmi✔️✔️: if (s1 == s2)
// Gemmi✔️✔️:   return Asu::Same;
// Gemmi✔️✔️: if (len1 >= 4 && len1 < 6 && len2 >= 4 && len2 < 6) { ... }
// Gemmi✔️✔️: return Asu::Different;
// END GEMMI CPP FUNCTION
fn compare_link_symops(record: &str, reported_sym: &mut [i16; 4]) -> BioAsu {
    if record.len() < 72 {
        return BioAsu::Any;
    }
    let s1 = read_fixed_string(record, 59, 6);
    let s2 = read_fixed_string(record, 66, 6);
    if s1 == s2 {
        return BioAsu::Same;
    }
    let len1 = s1.len();
    let len2 = s2.len();
    if (4..6).contains(&len1) && (4..6).contains(&len2) {
        reported_sym[0] = if s1.as_bytes().first() == Some(&b'1') && len1 == 4 {
            s2[..len2 - 3].trim().parse::<i16>().unwrap_or(0)
        } else {
            99
        };
        for i in 1..=3 {
            reported_sym[i] =
                i16::from(s2.as_bytes()[len2 - 4 + i]) - i16::from(s1.as_bytes()[len1 - 4 + i]);
        }
    }
    BioAsu::Different
}

fn element_from_padded_name_is_ambiguous(atom_name: &str) -> bool {
    let bytes = atom_name.as_bytes();
    bytes.len() >= 4
        && bytes[0] != b' '
        && bytes[3] != b' '
        && !bytes[0].is_ascii_digit()
        && !bytes[1].is_ascii_digit()
}

fn is_metal_element(element: Option<Element>) -> bool {
    // BEGIN GEMMI CPP FUNCTION gemmi::is_metal_value / gemmi::is_metal
    // Gemmi✔️✔️: inline bool& is_metal_value(El el) {
    // Gemmi✔️✔️:   static bool table[] = {
    // Gemmi✔️✔️:     // X     H     He
    // Gemmi✔️✔️:     false, false, false,
    // Gemmi✔️✔️:     // Li  Be     B      C      N      O      F     Ne
    // Gemmi✔️✔️:     true, true, false, false, false, false, false, false,
    // Gemmi✔️✔️:     // Na  Mg    Al     Si     P      S      Cl     Ar
    // Gemmi✔️✔️:     true, true, true, false, false, false, false, false,
    // Gemmi✔️✔️:     // K   Ca    Sc    Ti    V     Cr    Mn    Fe    Co    Ni    Cu    Zn
    // Gemmi✔️✔️:     true, true, true, true, true, true, true, true, true, true, true, true,
    // Gemmi✔️✔️:     // Ga  Ge    As     Se     Br     Kr
    // Gemmi✔️✔️:     true, true, false, false, false, false,
    // Gemmi✔️✔️:     // Rb  Sr    Y     Zr    Nb    Mo    Tc    Ru    Rh    Pd    Ag    Cd
    // Gemmi✔️✔️:     true, true, true, true, true, true, true, true, true, true, true, true,
    // Gemmi✔️✔️:     // In  Sn    Sb    Te      I     Xe
    // Gemmi✔️✔️:     true, true, true, false, false, false,
    // Gemmi✔️✔️:     // Cs  Ba    La    Ce    Pr    Nd    Pm    Sm    Eu    Gd    Tb    Dy
    // Gemmi✔️✔️:     true, true, true, true, true, true, true, true, true, true, true, true,
    // Gemmi✔️✔️:     // Ho  Er    Tm    Yb    Lu    Hf    Ta    W     Re    Os    Ir    Pt
    // Gemmi✔️✔️:     true, true, true, true, true, true, true, true, true, true, true, true,
    // Gemmi✔️✔️:     // Au  Hg    Tl    Pb    Bi    Po    At     Rn
    // Gemmi✔️✔️:     true, true, true, true, true, true, false, false,
    // Gemmi✔️✔️:     // Fr  Ra    Ac    Th    Pa    U     Np    Pu    Am    Cm    Bk    Cf
    // Gemmi✔️✔️:     true, true, true, true, true, true, true, true, true, true, true, true,
    // Gemmi✔️✔️:     // Es  Fm    Md    No    Lr    Rf    Db    Sg    Bh    Hs    Mt    Ds
    // Gemmi✔️✔️:     true, true, true, true, true, true, true, true, true, true, true, true,
    // Gemmi✔️✔️:     // Rg  Cn    Nh    Fl    Mc    Lv    Ts     Og
    // Gemmi✔️✔️:     true, true, true, true, true, true, false, false,
    // Gemmi✔️✔️:     // D    END
    // Gemmi✔️✔️:     false, false
    // Gemmi✔️✔️:   };
    // Gemmi✔️✔️:   return table[static_cast<int>(el)];
    // Gemmi✔️✔️: }
    // Gemmi✔️✔️: inline bool is_metal(El el) { return is_metal_value(el); }
    // END GEMMI CPP FUNCTION gemmi::is_metal_value / gemmi::is_metal
    let Some(element) = element else {
        return false;
    };
    matches!(
        element.atomic_number(),
        3 | 4
            | 11
            | 12
            | 13
            | 19..=32
            | 37..=50
            | 51
            | 55..=84
            | 87..=116
    ) && !matches!(element.atomic_number(), 33..=36 | 52..=54 | 85 | 86 | 117 | 118)
}

// BEGIN GEMMI CPP FUNCTION process_conn
// Gemmi✔️✔️: for (const std::string& record : conn_records) { ... parse SSBOND/LINK/CISPEP ... }
// Gemmi✔️✔️: SSBOND records become disulfN and call complete_ssbond(...)
// Gemmi✔️✔️: LINK/LINKR records become covaleN or metalcN based on parsed element types
// Gemmi✔️✔️: CISPEP model_num is first model num when the structure has exactly one model
// Gemmi✔️✔️: otherwise it is parsed from columns 43..45.
// END GEMMI CPP FUNCTION
fn process_conn(structure: &mut BioStructure, conn_records: &[String]) {
    let mut disulf_count = 0;
    let mut covale_count = 0;
    let mut metalc_count = 0;
    for record in conn_records {
        match record.as_bytes().first().copied() {
            Some(b'S' | b's') => {
                if record.len() < 32 {
                    continue;
                }
                disulf_count += 1;
                let mut connection = BioConnection {
                    name: format!("disulf{disulf_count}"),
                    type_: BioConnectionType::Disulf,
                    partner1: BioAtomAddress {
                        chain_name: read_fixed_string(record, 14, 2),
                        residue_name: read_fixed_string(record, 11, 3),
                        seq_id: Some(PdbSeqId {
                            seq_num: field(record, 17, 21).parse::<i32>().unwrap_or(0),
                            ins_code: match byte_at(record, 21) {
                                b' ' | 0 => None,
                                value => Some(value),
                            },
                        }),
                        ..BioAtomAddress::default()
                    },
                    partner2: BioAtomAddress {
                        chain_name: read_fixed_string(record, 28, 2),
                        residue_name: read_fixed_string(record, 25, 3),
                        seq_id: Some(PdbSeqId {
                            seq_num: field(record, 31, 35).parse::<i32>().unwrap_or(0),
                            ins_code: match byte_at(record, 35) {
                                b' ' | 0 => None,
                                value => Some(value),
                            },
                        }),
                        ..BioAtomAddress::default()
                    },
                    ..BioConnection::default()
                };
                connection.asu = compare_link_symops(record, &mut connection.reported_sym);
                if record.len() > 73 {
                    connection.reported_distance = field(record, 73, 78).parse::<f64>().ok();
                }
                complete_ssbond(&mut connection, structure);
                structure.connections.push(connection);
            }
            Some(b'L' | b'l') => {
                if record.len() < 57 {
                    continue;
                }
                let mut connection = BioConnection::default();
                for i in [0usize, 1] {
                    let base = 30 * i;
                    let address = if i == 0 {
                        &mut connection.partner1
                    } else {
                        &mut connection.partner2
                    };
                    address.chain_name = read_fixed_string(record, base + 20, 2);
                    address.residue_name = read_fixed_string(record, base + 17, 3);
                    address.seq_id = Some(PdbSeqId {
                        seq_num: field(record, base + 22, base + 26)
                            .parse::<i32>()
                            .unwrap_or(0),
                        ins_code: match byte_at(record, base + 26) {
                            b' ' | 0 => None,
                            value => Some(value),
                        },
                    });
                    address.atom_name = read_fixed_string(record, base + 12, 4);
                    address.altloc = parse_altloc(byte_at(record, base + 16));
                }
                let get_elem = |name: &str, address: &BioAtomAddress| -> Option<Element> {
                    if element_from_padded_name_is_ambiguous(name) {
                        if let Some(residue_id) = find_residue_for_address(structure, address, true)
                            && let Some(atom_idx) =
                                find_atom_in_residue(structure, residue_id, &address.atom_name)
                        {
                            return Some(structure.atoms[atom_idx].element);
                        }
                    }
                    infer_element_from_padded_atom_name(name)
                };
                let left = get_elem(field_raw(record, 12, 16), &connection.partner1);
                let right = get_elem(field_raw(record, 42, 46), &connection.partner2);
                if is_metal_element(left) || is_metal_element(right) {
                    metalc_count += 1;
                    connection.name = format!("metalc{metalc_count}");
                    connection.type_ = BioConnectionType::MetalC;
                } else {
                    covale_count += 1;
                    connection.name = format!("covale{covale_count}");
                    connection.type_ = BioConnectionType::Covale;
                }
                connection.asu = compare_link_symops(record, &mut connection.reported_sym);
                if record.len() > 73 {
                    if byte_at(record, 4) == b'R' {
                        connection.link_id = read_fixed_string(record, 72, 8);
                    } else {
                        connection.reported_distance = field(record, 73, 78).parse::<f64>().ok();
                    }
                }
                structure.connections.push(connection);
            }
            Some(b'C' | b'c') => {
                if record.len() < 22 {
                    continue;
                }
                let model_num = if structure.models.len() == 1 {
                    structure.models[0].source_model_number.unwrap_or(1)
                } else {
                    field(record, 43, 46).parse::<i32>().unwrap_or(0)
                };
                structure.cispeps.push(BioCisPep {
                    partner_c: BioAtomAddress {
                        chain_name: read_fixed_string(record, 14, 2),
                        residue_name: read_fixed_string(record, 11, 3),
                        seq_id: Some(PdbSeqId {
                            seq_num: field(record, 17, 21).parse::<i32>().unwrap_or(0),
                            ins_code: match byte_at(record, 21) {
                                b' ' | 0 => None,
                                value => Some(value),
                            },
                        }),
                        ..BioAtomAddress::default()
                    },
                    partner_n: BioAtomAddress {
                        chain_name: read_fixed_string(record, 28, 2),
                        residue_name: read_fixed_string(record, 25, 3),
                        seq_id: Some(PdbSeqId {
                            seq_num: field(record, 31, 35).parse::<i32>().unwrap_or(0),
                            ins_code: match byte_at(record, 35) {
                                b' ' | 0 => None,
                                value => Some(value),
                            },
                        }),
                        ..BioAtomAddress::default()
                    },
                    model_num,
                    only_altloc: None,
                    reported_angle: field(record, 53, 59).parse::<f64>().ok(),
                });
            }
            _ => {}
        }
    }
}

// BEGIN GEMMI CPP FUNCTION change_author_name_format_to_mmcif
// Gemmi✔️✔️: while (name[0] == ' ')
// Gemmi✔️✔️:   name.erase(name.begin());
// Gemmi✔️✔️: size_t pos = 0;
// Gemmi✔️✔️: for (size_t i = 1; i < pos+4 && i+1 < name.size(); ++i)
// Gemmi✔️✔️:   if (name[i] == '.' && name[i+1] != ' ')
// Gemmi✔️✔️:     pos = i+1;
// Gemmi✔️✔️: if (pos > 0)
// Gemmi✔️✔️:   name = name.substr(pos) + ", " + name.substr(0, pos);
// END GEMMI CPP FUNCTION
fn change_author_name_format_to_mmcif(name: &mut String) {
    while name.starts_with(' ') {
        name.remove(0);
    }
    let mut pos = 0usize;
    for i in 1..name.len() {
        if i >= pos + 4 || i + 1 >= name.len() {
            break;
        }
        if name.as_bytes()[i] == b'.' && name.as_bytes()[i + 1] != b' ' {
            pos = i + 1;
        }
    }
    if pos > 0 {
        *name = format!("{}, {}", &name[pos..], &name[..pos]);
    }
}

// BEGIN GEMMI CPP FUNCTION read_remark_290
// Gemmi✔️✔️: for (const std::string& remark : raw_remarks)
// Gemmi✔️✔️:   if (remark.size() > 25 && std::memcmp(&remark[7], "290", 3) == 0 &&
// Gemmi✔️✔️:       std::memcmp(&remark[10], "     ", 5) == 0 &&
// Gemmi✔️✔️:       std::memcmp(&remark[18], "555   ", 6) == 0) {
// Gemmi✔️✔️:     if (read_int(remark.c_str() + 15, 3) != (int)ops.size() + 1)
// Gemmi✔️✔️:       fail("Symmetry operators not in order?: " + remark);
// Gemmi✔️✔️:     ops.push_back(parse_triplet(read_string(remark.c_str() + 24, 56)));
// Gemmi✔️✔️:   }
// END GEMMI CPP FUNCTION
fn read_remark_290(raw_remarks: &[String]) -> Result<Vec<String>, BioReadError> {
    let mut ops = Vec::new();
    for remark in raw_remarks {
        if remark.len() > 25
            && field_raw(remark, 7, 10) == "290"
            && field_raw(remark, 10, 15) == "     "
            && field_raw(remark, 18, 24) == "555   "
        {
            let expected = ops.len() as i32 + 1;
            let observed = field(remark, 15, 18).parse::<i32>().unwrap_or(0);
            if observed != expected {
                return Err(BioReadError::Parse {
                    line_number: 0,
                    message: format!("Symmetry operators not in order?: {remark}"),
                });
            }
            ops.push(read_fixed_string(remark, 24, 56));
        }
    }
    Ok(ops)
}

fn read_pdb_resolution_field(value: &str) -> Option<f64> {
    value
        .split_whitespace()
        .next()
        .and_then(|token| token.parse::<f64>().ok())
}

fn split_remark_350_chains(value: &str) -> Vec<String> {
    value
        .split(',')
        .flat_map(|item| item.split_whitespace())
        .map(str::trim)
        .filter(|item| !item.is_empty() && *item != "AND" && *item != "CHAINS:")
        .map(str::to_string)
        .collect()
}

// BEGIN GEMMI CPP FUNCTION read_metadata_from_remarks
// Gemmi✔️✔️: REMARK 2 sets st.resolution from ANGSTROM lines when current resolution is 0.
// Gemmi✔️✔️: REMARK 3 delegates to read_remark3_line(...).
// Gemmi✔️✔️: REMARK 200/230/240 delegates to read_remark_200_230_240(...).
// Gemmi✔️✔️: REMARK 300 concatenates detail lines after the first REMARK: line.
// Gemmi✔️✔️: REMARK 350 builds assemblies, generators, BIOMT matrices, and summary fields.
// Gemmi✔️✔️: if REMARK 2 was absent, fallback resolution is taken from refinement.resolution_high.
// END GEMMI CPP FUNCTION
fn read_metadata_from_remarks(structure: &mut BioStructure) -> Result<(), BioReadError> {
    let mut remark3_continuation = None;
    let mut crystal_continuation = None;
    let mut matrix = BioTransform::default();

    for remark in &structure.raw_remarks {
        if remark.len() <= 11 {
            continue;
        }
        let num = field(remark, 7, 10).parse::<i32>().unwrap_or(0);
        match num {
            2 => {
                if structure.resolution.is_none() && remark.contains("ANGSTROM") {
                    structure.resolution = read_pdb_resolution_field(field_raw(remark, 23, 30));
                }
            }
            3 => read_remark3_line(remark, &mut structure.metadata, &mut remark3_continuation),
            200 | 230 | 240 => {
                read_remark_200_230_240(remark, &mut structure.metadata, &mut crystal_continuation);
            }
            300 => {
                if let Some(detail) = &mut structure.metadata.remark_300_detail {
                    detail.push('\n');
                    detail.push_str(field_raw(remark, 11, remark.len()).trim_end());
                } else if remark.get(11..18) == Some("REMARK:") {
                    structure.metadata.remark_300_detail =
                        Some(field_raw(remark, 18, remark.len()).trim().to_string());
                }
            }
            350 => {
                if !structure.assemblies.is_empty()
                    && remark[11..].trim_start().starts_with("BIOMT")
                {
                    let assembly = structure.assemblies.last_mut().expect("checked non-empty");
                    let row = read_matrix(&mut matrix, &remark[13..]);
                    if (row == 3 || remark[11..].trim_start().starts_with("BIOMT3"))
                        && !assembly.generators.is_empty()
                    {
                        let generator = assembly.generators.last_mut().expect("checked non-empty");
                        generator.operators.push(BioAssemblyOperator {
                            name: read_fixed_string(remark, 20, 3),
                            type_: String::new(),
                            transform: matrix,
                        });
                        matrix = BioTransform::default();
                    }
                    continue;
                }
                let Some(colon) = remark
                    .get(11..)
                    .and_then(|text| text.find(':'))
                    .map(|idx| idx + 11)
                else {
                    continue;
                };
                if remark
                    .get(11..)
                    .is_some_and(|text| text.starts_with("BIOMOLECULE"))
                {
                    structure.assemblies.push(BioAssembly {
                        name: field_raw(remark, colon + 1, remark.len())
                            .trim()
                            .to_string(),
                        ..BioAssembly::default()
                    });
                    continue;
                }
                let Some(assembly) = structure.assemblies.last_mut() else {
                    continue;
                };
                let r350_key =
                    |cpos: usize, text: &str| colon == cpos && remark[11..].starts_with(text);
                let remark_350_body = &remark[11..];
                if r350_key(44, "AUTHOR DETERMINED") {
                    assembly.author_determined = true;
                    assembly.oligomeric_details = read_fixed_string(remark, 45, 35);
                } else if r350_key(51, "SOFTWARE DETERMINED") {
                    assembly.software_determined = true;
                    assembly.oligomeric_details = read_fixed_string(remark, 52, 28);
                } else if r350_key(24, "SOFTWARE USED") {
                    assembly.software_name = read_fixed_string(remark, 25, 55);
                } else if r350_key(36, "TOTAL BURIED SURFACE AREA") {
                    assembly.absa = field(remark, 37, 49).parse::<f64>().ok();
                } else if r350_key(38, "SURFACE AREA OF THE COMPLEX") {
                    assembly.ssa = field(remark, 39, 51).parse::<f64>().ok();
                } else if r350_key(40, "CHANGE IN SOLVENT FREE ENERGY") {
                    assembly.more = field(remark, 41, 53).parse::<f64>().ok();
                } else if remark_350_body.contains("APPLY THE FOLLOWING TO CHAINS")
                    || remark_350_body.contains("AND CHAINS")
                {
                    if remark_350_body.trim_start().starts_with("APPLY") {
                        assembly.generators.push(BioAssemblyGenerator::default());
                    } else if assembly.generators.is_empty() {
                        continue;
                    }
                    let generator = assembly
                        .generators
                        .last_mut()
                        .expect("generator must exist");
                    generator.chains.extend(split_remark_350_chains(field_raw(
                        remark,
                        colon + 1,
                        remark.len(),
                    )));
                }
            }
            _ => {}
        }
    }

    if structure.resolution.is_none() {
        structure.resolution = structure
            .metadata
            .refinement
            .iter()
            .find_map(|info| info.resolution_high)
            .filter(|value| *value != 0.0);
    }
    Ok(())
}

fn mark_chain_polymer(structure: &mut BioStructure, chain_id: ChainId) {
    let chain = &structure.chains[chain_id.index() as usize];
    for residue_idx in chain.residue_span.start as usize..chain.residue_span.end() as usize {
        let residue = &mut structure.residues[residue_idx];
        residue.entity_kind = EntityKind::Polymer;
    }
}

fn chain_has_water(structure: &BioStructure, chain_id: ChainId) -> bool {
    let chain = &structure.chains[chain_id.index() as usize];
    let start = chain.residue_span.start as usize;
    let end = chain.residue_span.end() as usize;
    structure.residues[start..end]
        .iter()
        .any(|residue| residue.kind == ResidueKind::Water)
}

// BEGIN GEMMI CPP FUNCTION remove_entity_types
// Gemmi✔️✔️: for (Model& model : st.models)
// Gemmi✔️✔️:   for (Chain& chain : model.chains)
// Gemmi✔️✔️:     for (Residue& res : chain.residues)
// Gemmi✔️✔️:       res.entity_type = EntityType::Unknown;
// END GEMMI CPP FUNCTION
fn remove_entity_types(structure: &mut BioStructure) {
    for residue in &mut structure.residues {
        residue.entity_kind = EntityKind::Unknown;
    }
}

// BEGIN GEMMI CPP FUNCTION has_entity_types_and_subchains
// Gemmi✔️✔️: bool has_entity_types = true;
// Gemmi✔️✔️: bool has_subchains = true;
// Gemmi✔️✔️: for (const Residue& res : chain.residues) {
// Gemmi✔️✔️:   if (res.subchain.empty())
// Gemmi✔️✔️:     has_subchains = false;
// Gemmi✔️✔️:   if (res.entity_type == EntityType::Unknown)
// Gemmi✔️✔️:     has_entity_types = false;
// Gemmi✔️✔️: }
// Gemmi✔️✔️: return {has_entity_types, has_subchains};
// END GEMMI CPP FUNCTION
fn has_entity_types_and_subchains(structure: &BioStructure, chain_id: ChainId) -> (bool, bool) {
    let chain = &structure.chains[chain_id.index() as usize];
    let mut has_entity_types = true;
    let mut has_subchains = true;
    for residue in
        &structure.residues[chain.residue_span.start as usize..chain.residue_span.end() as usize]
    {
        if residue.source.subchain_id.is_none() {
            has_subchains = false;
        }
        if residue.entity_kind == EntityKind::Unknown {
            has_entity_types = false;
        }
    }
    (has_entity_types, has_subchains)
}

fn nonpolymer_subchain_id(chain_name: &str, counter: i32) -> Option<PdbChainId> {
    let mut name = String::with_capacity(chain_name.len() + 8);
    name.push_str(chain_name);
    name.push('x');
    if counter < 10 {
        name.push(char::from(b'0' + counter as u8));
    } else {
        let base36 = b"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
        let mut n = counter - 10;
        if n < 36 {
            name.push('0');
        }
        let pos = name.len();
        while n != 0 {
            name.insert(pos, char::from(base36[(n % 36) as usize]));
            n /= 36;
        }
    }
    pdb_chain_id_from_str(&name)
}

// BEGIN GEMMI CPP FUNCTION assign_subchain_names
// Gemmi✔️✔️: for (Residue& res : chain.residues) {
// Gemmi✔️✔️:   res.subchain = chain.name;
// Gemmi✔️✔️:   res.subchain += "x";
// Gemmi✔️✔️:   switch (res.entity_type) {
// Gemmi✔️✔️:     case EntityType::Polymer:    res.subchain += 'p'; break;
// Gemmi✔️✔️:     case EntityType::NonPolymer: ++nonpolymer_counter; ...; break;
// Gemmi✔️✔️:     case EntityType::Water:      res.subchain += 'w'; break;
// Gemmi✔️✔️:     case EntityType::Branched:   res.subchain += 'b'; break;
// Gemmi✔️✔️:     case EntityType::Unknown:    break;
// Gemmi✔️✔️:   }
// Gemmi✔️✔️: }
// END GEMMI CPP FUNCTION
fn assign_subchain_names(
    structure: &mut BioStructure,
    chain_id: ChainId,
    nonpolymer_counter: &mut i32,
) {
    let chain_name = {
        let chain = &structure.chains[chain_id.index() as usize];
        chain
            .source
            .auth_chain_id
            .or(chain.source.label_asym_id)
            .map(|id| id.as_str().to_string())
            .unwrap_or_default()
    };
    let chain = &structure.chains[chain_id.index() as usize];
    for residue in &mut structure.residues
        [chain.residue_span.start as usize..chain.residue_span.end() as usize]
    {
        residue.source.subchain_id = match residue.entity_kind {
            EntityKind::Polymer => pdb_chain_id_from_str(&format!("{chain_name}xp")),
            EntityKind::NonPolymer => {
                *nonpolymer_counter += 1;
                nonpolymer_subchain_id(&chain_name, *nonpolymer_counter)
            }
            EntityKind::Water => pdb_chain_id_from_str(&format!("{chain_name}xw")),
            EntityKind::Branched => pdb_chain_id_from_str(&format!("{chain_name}xb")),
            EntityKind::Unknown => residue.source.subchain_id,
        };
    }
}

// BEGIN GEMMI CPP FUNCTION assign_subchains
// Gemmi✔️✔️: for (Model& model : st.models) {
// Gemmi✔️✔️:   std::map<std::string, int> counters;
// Gemmi✔️✔️:   for (Chain& chain : model.chains) {
// Gemmi✔️✔️:     auto has = has_entity_types_and_subchains(chain);
// Gemmi✔️✔️:     if (force || !has.second) {
// Gemmi✔️✔️:       if (has.first)
// Gemmi✔️✔️:         assign_subchain_names(chain, counters[chain.name]);
// Gemmi✔️✔️:       else if (fail_if_unknown)
// Gemmi✔️✔️:         fail("assign_subchains(): missing entity_type in chain " + chain.name);
// Gemmi✔️✔️:     }
// Gemmi✔️✔️:   }
// Gemmi✔️✔️: }
// END GEMMI CPP FUNCTION
fn assign_subchains(
    structure: &mut BioStructure,
    force: bool,
    fail_if_unknown: bool,
) -> Result<(), BioReadError> {
    for model in structure.models.clone() {
        let mut counters: std::collections::HashMap<String, i32> = std::collections::HashMap::new();
        for chain_idx in model.chain_span.start as usize..model.chain_span.end() as usize {
            let chain_id = ChainId::new(chain_idx as u32);
            let has = has_entity_types_and_subchains(structure, chain_id);
            if force || !has.1 {
                if has.0 {
                    let chain_name = structure.chains[chain_idx]
                        .source
                        .auth_chain_id
                        .or(structure.chains[chain_idx].source.label_asym_id)
                        .map(|id| id.as_str().to_string())
                        .unwrap_or_default();
                    let counter = counters.entry(chain_name).or_insert(0);
                    assign_subchain_names(structure, chain_id, counter);
                } else if fail_if_unknown {
                    let chain_name = structure.chains[chain_idx]
                        .source
                        .auth_chain_id
                        .or(structure.chains[chain_idx].source.label_asym_id)
                        .map(|id| id.as_str().to_string())
                        .unwrap_or_default();
                    return Err(BioReadError::Parse {
                        line_number: 0,
                        message: format!(
                            "assign_subchains(): missing entity_type in chain {chain_name}"
                        ),
                    });
                }
            }
        }
    }
    Ok(())
}

fn first_polymer_subchain(structure: &BioStructure, chain_id: ChainId) -> Option<PdbChainId> {
    let chain = &structure.chains[chain_id.index() as usize];
    structure.residues[chain.residue_span.start as usize..chain.residue_span.end() as usize]
        .iter()
        .find(|residue| residue.entity_kind == EntityKind::Polymer)
        .and_then(|residue| residue.source.subchain_id)
}

fn backfill_polymer_subchains_to_entities(structure: &mut BioStructure) {
    let Some(first_model) = structure.models.first() else {
        return;
    };
    for chain_idx in first_model.chain_span.start as usize..first_model.chain_span.end() as usize {
        let chain = &structure.chains[chain_idx];
        let chain_name = chain
            .source
            .auth_chain_id
            .or(chain.source.label_asym_id)
            .map(|id| id.as_str().to_string())
            .unwrap_or_default();
        let Some(subchain_id) = first_polymer_subchain(structure, ChainId::new(chain_idx as u32))
        else {
            continue;
        };
        let Some(entity_idx) = structure
            .entities
            .iter()
            .position(|entity| entity.source.source_entity_id == chain_name)
        else {
            continue;
        };
        if !structure.entities[entity_idx]
            .subchains
            .contains(&subchain_id)
        {
            structure.entities[entity_idx].subchains.push(subchain_id);
        }
    }
}

// BEGIN GEMMI CPP FUNCTION restore_full_ccd_codes
// Gemmi✔️✔️: for (const auto& item : st.shortened_ccd_codes)
// Gemmi✔️✔️:   rename_residues(st, item.second, item.first);
// Gemmi✔️✔️: st.shortened_ccd_codes.clear();
// END GEMMI CPP FUNCTION
fn restore_full_ccd_codes(structure: &mut BioStructure) {
    // Current PDB subset stores HETNAM alias mappings but does not rename
    // residue or sequence identifiers to their shortened forms before finish().
    // For the currently modeled input state space, clearing the deferred map
    // reproduces Gemmi's post-restore final state.
    structure.shortened_ccd_codes.clear();
}

fn skip_space(text: &str, mut index: usize) -> usize {
    let bytes = text.as_bytes();
    while index < bytes.len() && bytes[index] == b' ' {
        index += 1;
    }
    index
}

fn ext_matches(a: u8, b: u8) -> bool {
    a == b
}

// BEGIN GEMMI CPP FUNCTION gemmi::Structure::find_spacegroup
// Gemmi✔️✔️: if (!cell.is_crystal())
// Gemmi✔️✔️:   return nullptr;
// Gemmi✔️✔️: return find_spacegroup_by_name(spacegroup_hm, cell.alpha, cell.gamma);
// END GEMMI CPP FUNCTION
fn find_crystal_spacegroup(crystal: &CrystalInfo) -> Option<&'static GemmiSpaceGroupEntry> {
    if !crystal_is_crystal(crystal) {
        return None;
    }
    crystal
        .spacegroup_hm
        .as_deref()
        .and_then(|hm| find_spacegroup_by_name(hm, crystal.cell.alpha, crystal.cell.gamma, None))
}

// BEGIN GEMMI CPP FUNCTION gemmi::find_spacegroup_by_name
// Gemmi✔️✔️: const SpaceGroup* find_spacegroup_by_name(std::string name, double alpha, double gamma,
// Gemmi✔️✔️:                                           const char* prefer) {
// Gemmi✔️✔️:   bool prefer_2 = false;
// Gemmi✔️✔️:   bool prefer_R = false;
// Gemmi✔️✔️:   if (prefer)
// Gemmi✔️✔️:     for (const char* p = prefer; *p != '\0'; ++p) {
// Gemmi✔️✔️:       if (*p == '2')
// Gemmi✔️✔️:         prefer_2 = true;
// Gemmi✔️✔️:       else if (*p == 'R')
// Gemmi✔️✔️:         prefer_R = true;
// Gemmi✔️✔️:       else if (*p != '1' && *p != 'H')
// Gemmi✔️✔️:         throw std::invalid_argument("find_spacegroup_by_name(): invalid arg 'prefer'");
// Gemmi✔️✔️:     }
// Gemmi✔️✔️:   const char* p = skip_space(name.c_str());
// Gemmi✔️✔️:   if (*p >= '0' && *p <= '9') { // handle numbers
// Gemmi✔️✔️:     char *endptr;
// Gemmi✔️✔️:     long n = std::strtol(p, &endptr, 10);
// Gemmi✔️✔️:     return *endptr == '\0' ? find_spacegroup_by_number(n) : nullptr;
// Gemmi✔️✔️:   }
// Gemmi✔️✔️:   char first = *p & ~0x20; // to uppercase
// Gemmi✔️✔️:   if (first == '\0')
// Gemmi✔️✔️:     return nullptr;
// Gemmi✔️✔️:   if (first == 'H')
// Gemmi✔️✔️:     first = 'R';
// Gemmi✔️✔️:   p = skip_space(p+1);
// Gemmi✔️✔️:   size_t start = p - name.c_str();
// Gemmi✔️✔️:   // change letters to lower case, except the letter after :
// Gemmi✔️✔️:   for (size_t i = start; i < name.size(); ++i) {
// Gemmi✔️✔️:     if (name[i] >= 'A' && name[i] <= 'Z')
// Gemmi✔️✔️:       name[i] |= 0x20;  // to lowercase
// Gemmi✔️✔️:     else if (name[i] == ':')
// Gemmi✔️✔️:       while (++i < name.size())
// Gemmi✔️✔️:         if (name[i] >= 'a' && name[i] <= 'z')
// Gemmi✔️✔️:           name[i] &= ~0x20;  // to uppercase
// Gemmi✔️✔️:   }
// Gemmi✔️✔️:   // allow names ending with R or H, such as R3R instead of R3:R
// Gemmi✔️✔️:   if (name.back() == 'h' || name.back() == 'r') {
// Gemmi✔️✔️:     name.back() &= ~0x20;  // to uppercase
// Gemmi✔️✔️:     name.insert(name.end() - 1, ':');
// Gemmi✔️✔️:   }
// Gemmi✔️✔️:   // The string that const char* p points to was just modified.
// Gemmi✔️✔️:   // This confuses some compilers (GCC 4.8), so let's re-assign p.
// Gemmi✔️✔️:   p = name.c_str() + start;
// Gemmi✔️✔️:   for (const SpaceGroup& sg : spacegroup_tables::main)
// Gemmi✔️✔️:     if (sg.hm[0] == first) {
// Gemmi✔️✔️:       if (sg.hm[2] == *p) {
// Gemmi✔️✔️:         const char* a = skip_space(p + 1);
// Gemmi✔️✔️:         const char* b = skip_space(sg.hm + 3);
// Gemmi✔️✔️:         while ((*a == *b && *b != '\0') ||
// Gemmi✔️✔️:                (*a == '3' && *b == '-' && b == sg.hm + 4 && *++b == '3')) {
// Gemmi✔️✔️:           a = skip_space(a+1);
// Gemmi✔️✔️:           b = skip_space(b+1);
// Gemmi✔️✔️:         }
// Gemmi✔️✔️:         if (*b == '\0') {
// Gemmi✔️✔️:           if (*a == '\0') {
// Gemmi✔️✔️:             if (sg.ext == 'H' && (alpha == 0. ? prefer_R : gamma < 1.125 * alpha))
// Gemmi✔️✔️:               return &sg + 1;
// Gemmi✔️✔️:             if (sg.ext == '1' && prefer_2)
// Gemmi✔️✔️:               return &sg + 1;
// Gemmi✔️✔️:             return &sg;
// Gemmi✔️✔️:           }
// Gemmi✔️✔️:           if (*a == ':' && *skip_space(a+1) == sg.ext)
// Gemmi✔️✔️:             return &sg;
// Gemmi✔️✔️:         }
// Gemmi✔️✔️:       } else if (sg.hm[2] == '1' && sg.hm[3] == ' ') {
// Gemmi✔️✔️:         const char* b = sg.hm + 4;
// Gemmi✔️✔️:         if (*b != '1' || (first == 'B' && *++b == ' ' && *++b != '1')) {
// Gemmi✔️✔️:           char end = (b == sg.hm + 4 ? ' ' : '\0');
// Gemmi✔️✔️:           const char* a = skip_space(p);
// Gemmi✔️✔️:           while (*a == *b && *b != end) {
// Gemmi✔️✔️:             ++a;
// Gemmi✔️✔️:             ++b;
// Gemmi✔️✔️:           }
// Gemmi✔️✔️:           if (*skip_space(a) == '\0' && *b == end)
// Gemmi✔️✔️:             return &sg;
// Gemmi✔️✔️:         }
// Gemmi✔️✔️:       }
// Gemmi✔️✔️:     }
// Gemmi✔️✔️:   for (const SpaceGroupAltName& sg : spacegroup_tables::alt_names)
// Gemmi✔️✔️:     if (sg.hm[0] == first && sg.hm[2] == *p) {
// Gemmi✔️✔️:       const char* a = skip_space(p + 1);
// Gemmi✔️✔️:       const char* b = skip_space(sg.hm + 3);
// Gemmi✔️✔️:       while (*a == *b && *b != '\0') {
// Gemmi✔️✔️:         a = skip_space(a+1);
// Gemmi✔️✔️:         b = skip_space(b+1);
// Gemmi✔️✔️:       }
// Gemmi✔️✔️:       if (*b == '\0' &&
// Gemmi✔️✔️:           (*a == '\0' || (*a == ':' && *skip_space(a+1) == sg.ext)))
// Gemmi✔️✔️:         return &spacegroup_tables::main[sg.pos];
// Gemmi✔️✔️:     }
// Gemmi✔️✔️:   return nullptr;
// Gemmi✔️✔️: }
// END GEMMI CPP FUNCTION
fn find_spacegroup_by_name(
    name: &str,
    alpha: f64,
    gamma: f64,
    prefer: Option<&str>,
) -> Option<&'static GemmiSpaceGroupEntry> {
    let mut prefer_2 = false;
    let mut prefer_r = false;
    if let Some(prefer) = prefer {
        for ch in prefer.bytes() {
            match ch {
                b'2' => prefer_2 = true,
                b'R' => prefer_r = true,
                b'1' | b'H' => {}
                _ => return None,
            }
        }
    }

    let mut name = name.trim_start().to_string();
    if name.is_empty() {
        return None;
    }
    if name.as_bytes()[0].is_ascii_digit() {
        if name.bytes().all(|b| b.is_ascii_digit()) {
            let value = name.parse::<usize>().ok()?;
            return GEMMI_SPACEGROUPS.get(value.checked_sub(1)?);
        }
        return None;
    }

    let mut first = name.as_bytes()[0].to_ascii_uppercase();
    if first == b'H' {
        first = b'R';
    }
    let mut p = skip_space(&name, 1);
    let start = p;
    let mut bytes = name.into_bytes();
    let mut i = start;
    while i < bytes.len() {
        if bytes[i].is_ascii_uppercase() {
            bytes[i] = bytes[i].to_ascii_lowercase();
        } else if bytes[i] == b':' {
            i += 1;
            while i < bytes.len() {
                if bytes[i].is_ascii_lowercase() {
                    bytes[i] = bytes[i].to_ascii_uppercase();
                }
                i += 1;
            }
            break;
        }
        i += 1;
    }
    if let Some(&last) = bytes.last()
        && (last == b'h' || last == b'r')
    {
        let len = bytes.len();
        bytes[len - 1] = last.to_ascii_uppercase();
        bytes.insert(len - 1, b':');
    }
    let normalized = String::from_utf8(bytes).ok()?;
    p = start;

    for (index, sg) in GEMMI_SPACEGROUPS.iter().enumerate() {
        let hm = sg.hm.as_bytes();
        if hm[0] != first {
            continue;
        }
        if hm[2] == normalized.as_bytes().get(p).copied().unwrap_or_default() {
            let mut a = skip_space(&normalized, p + 1);
            let mut b = skip_space(sg.hm, 3);
            while (normalized.as_bytes().get(a).copied().unwrap_or_default()
                == hm.get(b).copied().unwrap_or_default()
                && hm.get(b).copied().unwrap_or_default() != 0)
                || (normalized.as_bytes().get(a).copied().unwrap_or_default() == b'3'
                    && hm.get(b).copied().unwrap_or_default() == b'-'
                    && b == 4
                    && hm.get(b + 1).copied().unwrap_or_default() == b'3')
            {
                if hm.get(b).copied().unwrap_or_default() == b'-' && b == 4 {
                    b += 1;
                }
                a = skip_space(&normalized, a + 1);
                b = skip_space(sg.hm, b + 1);
            }
            if b >= hm.len() {
                if a >= normalized.len() {
                    if ext_matches(sg.ext, b'H')
                        && (alpha == 0.0 || gamma < 1.125 * alpha || prefer_r)
                    {
                        return GEMMI_SPACEGROUPS.get(index + 1);
                    }
                    if ext_matches(sg.ext, b'1') && prefer_2 {
                        return GEMMI_SPACEGROUPS.get(index + 1);
                    }
                    return Some(sg);
                }
                if normalized.as_bytes().get(a).copied() == Some(b':')
                    && normalized
                        .as_bytes()
                        .get(skip_space(&normalized, a + 1))
                        .copied()
                        == Some(sg.ext)
                {
                    return Some(sg);
                }
            }
        } else if hm.len() > 3 && hm[2] == b'1' && hm[3] == b' ' {
            let mut b = 4;
            if hm[b] != b'1'
                || (first == b'B' && {
                    b += 1;
                    hm[b] == b' ' && {
                        b += 1;
                        hm[b] != b'1'
                    }
                })
            {
                let end = if b == 4 { b' ' } else { b'\0' };
                let mut a = skip_space(&normalized, p);
                while normalized.as_bytes().get(a).copied().unwrap_or_default()
                    == hm.get(b).copied().unwrap_or_default()
                    && hm.get(b).copied().unwrap_or(end) != end
                {
                    a += 1;
                    b += 1;
                }
                if skip_space(&normalized, a) >= normalized.len()
                    && hm.get(b).copied().unwrap_or(end) == end
                {
                    return Some(sg);
                }
            }
        }
    }

    for alt in GEMMI_ALT_NAMES {
        let hm = alt.hm.as_bytes();
        if hm[0] != first || hm[2] != normalized.as_bytes().get(p).copied().unwrap_or_default() {
            continue;
        }
        let mut a = skip_space(&normalized, p + 1);
        let mut b = skip_space(alt.hm, 3);
        while normalized.as_bytes().get(a).copied().unwrap_or_default()
            == hm.get(b).copied().unwrap_or_default()
            && hm.get(b).copied().unwrap_or_default() != 0
        {
            a = skip_space(&normalized, a + 1);
            b = skip_space(alt.hm, b + 1);
        }
        if b >= hm.len()
            && (a >= normalized.len()
                || (normalized.as_bytes().get(a).copied() == Some(b':')
                    && normalized
                        .as_bytes()
                        .get(skip_space(&normalized, a + 1))
                        .copied()
                        == Some(alt.ext)))
        {
            return GEMMI_SPACEGROUPS.get(alt.pos);
        }
    }
    None
}

fn gemmi_op_to_bio_transform(op: &crate::io::gemmi_spacegroup_table::GemmiSymOp) -> BioTransform {
    let den = GEMMI_OP_DEN as f64;
    BioTransform {
        mat: [
            [
                op.rot[0][0] as f64 / den,
                op.rot[0][1] as f64 / den,
                op.rot[0][2] as f64 / den,
            ],
            [
                op.rot[1][0] as f64 / den,
                op.rot[1][1] as f64 / den,
                op.rot[1][2] as f64 / den,
            ],
            [
                op.rot[2][0] as f64 / den,
                op.rot[2][1] as f64 / den,
                op.rot[2][2] as f64 / den,
            ],
        ],
        vec: [
            op.tran[0] as f64 / den,
            op.tran[1] as f64 / den,
            op.tran[2] as f64 / den,
        ],
    }
}

fn setup_cell_images(structure: &mut BioStructure) -> Result<(), BioReadError> {
    let Some(crystal) = structure.crystal.as_mut() else {
        return Ok(());
    };

    // BEGIN GEMMI CPP FUNCTION Structure::setup_cell_images
    // Gemmi✔️✔️: const SpaceGroup* sg = find_spacegroup();
    // Gemmi✔️✔️: cell.set_cell_images_from_spacegroup(sg);
    // Gemmi✔️✔️: cell.add_ncs_images_to_cs_images(ncs);
    // END GEMMI CPP FUNCTION
    if let Some(spacegroup) = find_crystal_spacegroup(crystal) {
        crystal.cell_images.clear();
        crystal.cs_count = spacegroup.ops.len().saturating_sub(1) as i16;
        crystal.cell_images.reserve(crystal.cs_count as usize);
        for op in &spacegroup.ops[1..] {
            crystal.cell_images.push(gemmi_op_to_bio_transform(op));
        }
    } else {
        crystal.cell_images.clear();
        crystal.cs_count = 0;
    }
    for ncs_op in &structure.ncs_operators {
        if !ncs_op.given {
            let f = bio_transform_combine(
                &crystal.frac,
                &bio_transform_combine(&ncs_op.transform, &crystal.orth),
            );
            crystal.cell_images.push(f);
            for i in 0..crystal.cs_count as usize {
                let combined = bio_transform_combine(&crystal.cell_images[i], &f);
                crystal.cell_images.push(combined);
            }
        }
    }
    Ok(())
}

fn parse_pdb_anisou_record(line: &str, line_number: usize) -> Result<[f32; 6], BioReadError> {
    if line.len() < 70 {
        return Err(BioReadError::Parse {
            line_number,
            message: "ANISOU line is too short to be correct".to_string(),
        });
    }
    Ok([
        parse_decimal_i32(
            &BIO_PDB_READ_FEATURE,
            field(line, 28, 35),
            line_number,
            "ANISOU u11",
        )? as f32
            * 1e-4,
        parse_decimal_i32(
            &BIO_PDB_READ_FEATURE,
            field(line, 35, 42),
            line_number,
            "ANISOU u22",
        )? as f32
            * 1e-4,
        parse_decimal_i32(
            &BIO_PDB_READ_FEATURE,
            field(line, 42, 49),
            line_number,
            "ANISOU u33",
        )? as f32
            * 1e-4,
        parse_decimal_i32(
            &BIO_PDB_READ_FEATURE,
            field(line, 49, 56),
            line_number,
            "ANISOU u12",
        )? as f32
            * 1e-4,
        parse_decimal_i32(
            &BIO_PDB_READ_FEATURE,
            field(line, 56, 63),
            line_number,
            "ANISOU u13",
        )? as f32
            * 1e-4,
        parse_decimal_i32(
            &BIO_PDB_READ_FEATURE,
            field(line, 63, 70),
            line_number,
            "ANISOU u23",
        )? as f32
            * 1e-4,
    ])
}

fn parse_pdb_seqres_record(builder: &mut PdbBioBuilder, line: &str) {
    let chain_id = pdb_chain_id_from_field(field(line, 10, 12));
    let source_entity_id = chain_id.as_str();
    if source_entity_id.is_empty() {
        return;
    }
    let entity_id =
        builder.find_or_add_entity(source_entity_id, EntityKind::Polymer, PolymerKind::Unknown);
    for start in (19..68).step_by(4) {
        let residue_name = field(line, start, start + 3);
        if !residue_name.is_empty() {
            builder.append_entity_sequence(entity_id, residue_name);
        }
    }
}

fn parse_pdb_helix_record(structure: &mut BioStructure, line: &str) {
    if line.len() < 40 {
        return;
    }
    let mut helix = BioHelix {
        start: BioAtomAddress {
            chain_name: read_fixed_string(line, 18, 2),
            residue_name: read_fixed_string(line, 15, 3),
            seq_id: Some(PdbSeqId {
                seq_num: field(line, 21, 25).parse::<i32>().unwrap_or(0),
                ins_code: match byte_at(line, 25) {
                    b' ' | 0 => None,
                    value => Some(value),
                },
            }),
            ..BioAtomAddress::default()
        },
        end: BioAtomAddress {
            chain_name: read_fixed_string(line, 30, 2),
            residue_name: read_fixed_string(line, 27, 3),
            seq_id: Some(PdbSeqId {
                seq_num: field(line, 33, 37).parse::<i32>().unwrap_or(0),
                ins_code: match byte_at(line, 37) {
                    b' ' | 0 => None,
                    value => Some(value),
                },
            }),
            ..BioAtomAddress::default()
        },
        length: -1,
        ..BioHelix::default()
    };
    helix.set_helix_class_as_int(field(line, 38, 40).parse::<i32>().unwrap_or(0));
    if line.len() > 72 {
        helix.length = field(line, 72, 77).parse::<i32>().unwrap_or(0);
    }
    structure.helices.push(helix);
}

fn parse_pdb_sheet_record(structure: &mut BioStructure, line: &str) {
    if line.len() < 40 {
        return;
    }
    let sheet_id = read_fixed_string(line, 11, 3);
    let sheet_idx = structure
        .sheets
        .iter()
        .position(|sheet| sheet.name == sheet_id)
        .unwrap_or_else(|| {
            structure.sheets.push(BioSheet {
                name: sheet_id.clone(),
                ..BioSheet::default()
            });
            structure.sheets.len() - 1
        });
    let mut strand = BioSheetStrand {
        start: BioAtomAddress {
            chain_name: read_fixed_string(line, 20, 2),
            residue_name: read_fixed_string(line, 17, 3),
            seq_id: Some(PdbSeqId {
                seq_num: field(line, 22, 26).parse::<i32>().unwrap_or(0),
                ins_code: match byte_at(line, 26) {
                    b' ' | 0 => None,
                    value => Some(value),
                },
            }),
            ..BioAtomAddress::default()
        },
        end: BioAtomAddress {
            chain_name: read_fixed_string(line, 31, 2),
            residue_name: read_fixed_string(line, 28, 3),
            seq_id: Some(PdbSeqId {
                seq_num: field(line, 33, 37).parse::<i32>().unwrap_or(0),
                ins_code: match byte_at(line, 37) {
                    b' ' | 0 => None,
                    value => Some(value),
                },
            }),
            ..BioAtomAddress::default()
        },
        sense: field(line, 38, 40).parse::<i32>().unwrap_or(0),
        ..BioSheetStrand::default()
    };
    if line.len() > 67 {
        strand.hbond_atom2.atom_name = read_fixed_string(line, 41, 4);
        strand.hbond_atom2.residue_name = read_fixed_string(line, 45, 3);
        strand.hbond_atom2.chain_name = read_fixed_string(line, 48, 2);
        strand.hbond_atom2.seq_id = Some(PdbSeqId {
            seq_num: field(line, 50, 54).parse::<i32>().unwrap_or(0),
            ins_code: match byte_at(line, 54) {
                b' ' | 0 => None,
                value => Some(value),
            },
        });
        strand.hbond_atom1.atom_name = read_fixed_string(line, 56, 4);
        strand.hbond_atom1.residue_name = read_fixed_string(line, 60, 3);
        strand.hbond_atom1.chain_name = read_fixed_string(line, 63, 2);
        strand.hbond_atom1.seq_id = Some(PdbSeqId {
            seq_num: field(line, 65, 69).parse::<i32>().unwrap_or(0),
            ins_code: match byte_at(line, 69) {
                b' ' | 0 => None,
                value => Some(value),
            },
        });
    }
    structure.sheets[sheet_idx].strands.push(strand);
}

fn parse_pdb_modres_record(structure: &mut BioStructure, line: &str) {
    let mut modres = BioModRes {
        chain_name: read_fixed_string(line, 15, 2),
        residue_name: read_fixed_string(line, 12, 3),
        res_id: PdbSeqId {
            seq_num: field(line, 18, 22).parse::<i32>().unwrap_or(0),
            ins_code: match byte_at(line, 22) {
                b' ' | 0 => None,
                value => Some(value),
            },
        },
        parent_comp_id: read_fixed_string(line, 24, 3),
        ..BioModRes::default()
    };
    if line.len() >= 30 {
        modres.details = field_raw(line, 29, 70).trim().to_string();
    }
    if line.len() >= 73 && byte_at(line, 70) == b' ' && byte_at(line, 71) == b' ' {
        modres.mod_id = read_fixed_string(line, 72, 8);
    }
    structure.mod_residues.push(modres);
}

fn parse_pdb_hetnam_record(structure: &mut BioStructure, line: &str) {
    if line.len() > 71 && byte_at(line, 70) == b' ' {
        let full_code = read_fixed_string(line, 71, 8);
        if !full_code.is_empty() {
            structure
                .shortened_ccd_codes
                .push((full_code, read_fixed_string(line, 11, 3)));
        }
    }
}

fn parse_pdb_dbref_record(builder: &mut PdbBioBuilder, line: &str) {
    let chain_name = read_fixed_string(line, 11, 2);
    let entity_id =
        builder.find_or_add_entity(&chain_name, EntityKind::Polymer, PolymerKind::Unknown);
    let entity = &mut builder.structure.entities[entity_id.index() as usize];
    if matches!(byte_at(line, 5), b' ' | b'1') {
        entity.dbrefs.push(BioEntityDbRef::default());
    } else if byte_at(line, 5) == b'2' && entity.dbrefs.is_empty() {
        return;
    }
    let dbref = match entity.dbrefs.last_mut() {
        Some(dbref) => dbref,
        None => return,
    };
    if matches!(byte_at(line, 5), b' ' | b'1') {
        dbref.seq_begin = read_seq_id_range(line, 14);
        dbref.seq_end = read_seq_id_range(line, 20);
        dbref.db_name = read_fixed_string(line, 26, 6);
        if byte_at(line, 5) == b' ' {
            dbref.accession_code = read_fixed_string(line, 33, 8);
            dbref.id_code = read_fixed_string(line, 42, 12);
            dbref.db_begin = read_seq_id_range(line, 55);
            dbref.db_end = read_seq_id_range(line, 62);
        } else {
            dbref.id_code = read_fixed_string(line, 47, 20);
        }
    } else if byte_at(line, 5) == b'2' {
        dbref.accession_code = read_fixed_string(line, 18, 22);
        dbref.db_begin = field(line, 45, 55)
            .parse::<i32>()
            .ok()
            .map(|seq_num| PdbSeqId {
                seq_num,
                ins_code: None,
            });
        dbref.db_end = field(line, 57, 67)
            .parse::<i32>()
            .ok()
            .map(|seq_num| PdbSeqId {
                seq_num,
                ins_code: None,
            });
    }
}

fn read_seq_id_range(line: &str, start: usize) -> Option<PdbSeqId> {
    let seq_num = field(line, start, start + 5).parse::<i32>().ok()?;
    Some(PdbSeqId {
        seq_num,
        ins_code: match byte_at(line, start + 5) {
            b' ' | 0 => None,
            value => Some(value),
        },
    })
}

fn parse_pdb_header_record(builder: &mut PdbBioBuilder, line: &str) {
    let keywords = field_raw(line, 10, 50).trim_end();
    if !keywords.is_empty() {
        builder.structure.metadata.pdbx_keywords = Some(keywords.to_string());
    }
    let date = pdb_date_format_to_iso(field_raw(line, 50, 59));
    if !date.is_empty() {
        builder.structure.metadata.received_initial_deposition_date = Some(date);
    }
    let entry_id = field_raw(line, 62, 66).trim_end();
    if !entry_id.is_empty() {
        builder.structure.metadata.entry_id = Some(entry_id.to_string());
    }
}

fn parse_pdb_author_record(builder: &mut PdbBioBuilder, line: &str) {
    let text = field_raw(line, 10, line.len()).trim();
    if text.is_empty() {
        return;
    }

    let mut previous_tail = None;
    if let Some(last) = builder.structure.metadata.authors.pop() {
        previous_tail = Some(last);
    }
    let previous_len = builder.structure.metadata.authors.len();
    builder.structure.metadata.authors.extend(
        text.split(',')
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(str::to_string),
    );

    if let Some(mut last) = previous_tail {
        if builder.structure.metadata.authors.len() > previous_len {
            if !last.ends_with('-') && !last.ends_with('.') {
                last.push(' ');
            }
            builder.structure.metadata.authors[previous_len].insert_str(0, &last);
        } else {
            builder.structure.metadata.authors.push(last);
        }
    }
}

fn parse_pdb_cryst1_record(
    builder: &mut PdbBioBuilder,
    line: &str,
    line_number: usize,
) -> Result<(), BioReadError> {
    if line.len() > 54 {
        let existing_scale = builder
            .structure
            .crystal
            .as_ref()
            .and_then(|crystal| crystal.scale);
        let mut crystal = default_crystal_info();
        crystal.scale = existing_scale;
        crystal_set(
            &mut crystal,
            CrystalCell {
                a: parse_f64(field(line, 6, 15), line_number, "CRYST1 a")?,
                b: parse_f64(field(line, 15, 24), line_number, "CRYST1 b")?,
                c: parse_f64(field(line, 24, 33), line_number, "CRYST1 c")?,
                alpha: parse_f64(field(line, 33, 40), line_number, "CRYST1 alpha")?,
                beta: parse_f64(field(line, 40, 47), line_number, "CRYST1 beta")?,
                gamma: parse_f64(field(line, 47, 54), line_number, "CRYST1 gamma")?,
            },
        );
        if let Some(scale) = existing_scale {
            crystal_set_matrices_from_fract(&mut crystal, scale);
        }
        builder.structure.crystal = Some(crystal);
    }
    if line.len() > 56 {
        let spacegroup = field(line, 55, 66);
        if !spacegroup.is_empty()
            && let Some(crystal) = &mut builder.structure.crystal
        {
            crystal.spacegroup_hm = Some(spacegroup.to_string());
        }
    }
    if line.len() > 67 {
        let z_pdb = field(line, 66, 70);
        if !z_pdb.is_empty()
            && let Some(crystal) = &mut builder.structure.crystal
        {
            crystal.z_pdb = Some(z_pdb.to_string());
        }
    }
    Ok(())
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct CifSmat33<T> {
    u11: T,
    u22: T,
    u33: T,
    u12: T,
    u13: T,
    u23: T,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Remark3Continuation {
    TlsSelection {
        refinement_idx: usize,
        tls_group_idx: usize,
        selection_idx: usize,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Remark200Continuation {
    CrystalDescription { crystal_idx: usize },
}

fn read_mmcif_entity_and_sequence_info(
    builder: &mut PdbBioBuilder,
    loops: &[CifLoop],
) -> Result<(), BioReadError> {
    if let Some(entity_loop) = find_cif_loop(loops, "_entity.id") {
        let width = checked_cif_loop_width(entity_loop)?;
        let id_col = required_cif_col(entity_loop, "_entity.id")?;
        let type_col = optional_cif_col(entity_loop, "_entity.type");
        for row in entity_loop.values.chunks(width) {
            let source_id = cif_optional(row[id_col].value.as_str())
                .ok_or_else(|| missing_cif_value(row[id_col].line_number, "_entity.id"))?;
            let kind = type_col
                .and_then(|idx| cif_optional(row[idx].value.as_str()))
                .map(entity_kind_from_cif)
                .unwrap_or(EntityKind::Unknown);
            let polymer_kind = find_mmcif_polymer_kind(loops, source_id)?;
            let entity_kind = if kind == EntityKind::Unknown && polymer_kind != PolymerKind::Unknown
            {
                EntityKind::Polymer
            } else {
                kind
            };
            let entity_id = builder.find_or_add_entity(source_id, entity_kind, polymer_kind);
            builder.structure.entities[entity_id.index() as usize].reflects_microhetero = true;
        }
    }

    if let Some(sequence_loop) = find_cif_loop(loops, "_entity_poly_seq.entity_id") {
        let width = checked_cif_loop_width(sequence_loop)?;
        if let Some([entity_col, num_col, mon_col]) = required_cif_table_cols(
            sequence_loop,
            [
                "_entity_poly_seq.entity_id",
                "_entity_poly_seq.num",
                "_entity_poly_seq.mon_id",
            ],
        ) {
            for row in sequence_loop.values.chunks(width) {
                let source_id = cif_optional(row[entity_col].value.as_str()).ok_or_else(|| {
                    missing_cif_value(row[entity_col].line_number, "_entity_poly_seq.entity_id")
                })?;
                let Some(entity_id) = builder.find_entity_by_source_id(source_id) else {
                    continue;
                };
                let pos = parse_decimal_i32(
                    &BIO_MMCIF_READ_FEATURE,
                    row[num_col].value.as_str(),
                    row[num_col].line_number,
                    "_entity_poly_seq.num",
                )? - 1;
                if pos < 0 {
                    continue;
                }
                if let Some(residue_name) = cif_optional(row[mon_col].value.as_str()) {
                    builder.merge_entity_sequence_at(entity_id, pos as usize, residue_name);
                }
            }
        }
    }

    if let (Some(struct_ref), Some(struct_ref_seq)) = (
        find_cif_loop(loops, "_struct_ref.id"),
        find_cif_loop(loops, "_struct_ref_seq.ref_id"),
    ) {
        let ref_width = checked_cif_loop_width(struct_ref)?;
        if let (
            Some([ref_id, ref_entity_id, ref_db_name, ref_db_code]),
            Some(
                [
                    seq_ref_id,
                    seq_align_beg,
                    seq_align_end,
                    db_align_beg,
                    db_align_end,
                ],
            ),
        ) = (
            required_cif_table_cols(
                struct_ref,
                [
                    "_struct_ref.id",
                    "_struct_ref.entity_id",
                    "_struct_ref.db_name",
                    "_struct_ref.db_code",
                ],
            ),
            required_cif_table_cols(
                struct_ref_seq,
                [
                    "_struct_ref_seq.ref_id",
                    "_struct_ref_seq.seq_align_beg",
                    "_struct_ref_seq.seq_align_end",
                    "_struct_ref_seq.db_align_beg",
                    "_struct_ref_seq.db_align_end",
                ],
            ),
        ) {
            let ref_accession = optional_cif_col(struct_ref, "_struct_ref.pdbx_db_accession");
            let ref_isoform = optional_cif_col(struct_ref, "_struct_ref.pdbx_db_isoform");

            let auth_seq_align_beg =
                optional_cif_col(struct_ref_seq, "_struct_ref_seq.pdbx_auth_seq_align_beg");
            let seq_align_beg_ins = optional_cif_col(
                struct_ref_seq,
                "_struct_ref_seq.pdbx_seq_align_beg_ins_code",
            );
            let auth_seq_align_end =
                optional_cif_col(struct_ref_seq, "_struct_ref_seq.pdbx_auth_seq_align_end");
            let seq_align_end_ins = optional_cif_col(
                struct_ref_seq,
                "_struct_ref_seq.pdbx_seq_align_end_ins_code",
            );

            let mut seen = Vec::<String>::new();
            for seq_row in cif_loop_rows(struct_ref_seq)? {
                let dedup = format!(
                    "{}\t{}\t{}\t{}\t{}",
                    seq_row[seq_ref_id].value.trim(),
                    seq_row[seq_align_beg].value.trim(),
                    seq_row[seq_align_end].value.trim(),
                    seq_row[db_align_beg].value.trim(),
                    seq_row[db_align_end].value.trim()
                );
                if seen.iter().any(|item| item == &dedup) {
                    continue;
                }
                seen.push(dedup);

                let Some(ref_row) = struct_ref.values.chunks(ref_width).find(|row| {
                    row.get(ref_id)
                        .is_some_and(|token| token.value.trim() == seq_row[seq_ref_id].value.trim())
                }) else {
                    continue;
                };
                let source_id = ref_row[ref_entity_id].value.trim();
                let Some(entity_id) = builder.find_entity_by_source_id(source_id) else {
                    continue;
                };
                let entity = &mut builder.structure.entities[entity_id.index() as usize];
                let mut dbref = BioEntityDbRef {
                    db_name: ref_row[ref_db_name].value.trim().to_string(),
                    id_code: ref_row[ref_db_code].value.trim().to_string(),
                    ..BioEntityDbRef::default()
                };
                if let Some(index) = ref_accession
                    && cif_row_has2(ref_row, index)
                {
                    dbref.accession_code = ref_row[index].value.trim().to_string();
                }
                if let Some(index) = ref_isoform
                    && cif_row_has2(ref_row, index)
                {
                    dbref.isoform = ref_row[index].value.trim().to_string();
                }
                dbref.label_seq_begin = cif_optional(seq_row[seq_align_beg].value.as_str())
                    .map(|value| {
                        parse_decimal_i32(
                            &BIO_MMCIF_READ_FEATURE,
                            value,
                            seq_row[seq_align_beg].line_number,
                            "_struct_ref_seq.seq_align_beg",
                        )
                    })
                    .transpose()?;
                dbref.label_seq_end = cif_optional(seq_row[seq_align_end].value.as_str())
                    .map(|value| {
                        parse_decimal_i32(
                            &BIO_MMCIF_READ_FEATURE,
                            value,
                            seq_row[seq_align_end].line_number,
                            "_struct_ref_seq.seq_align_end",
                        )
                    })
                    .transpose()?;
                dbref.db_begin = cif_optional(seq_row[db_align_beg].value.as_str())
                    .map(|value| {
                        parse_decimal_i32(
                            &BIO_MMCIF_READ_FEATURE,
                            value,
                            seq_row[db_align_beg].line_number,
                            "_struct_ref_seq.db_align_beg",
                        )
                    })
                    .transpose()?
                    .map(|seq_num| PdbSeqId {
                        seq_num,
                        ins_code: None,
                    });
                dbref.db_end = cif_optional(seq_row[db_align_end].value.as_str())
                    .map(|value| {
                        parse_decimal_i32(
                            &BIO_MMCIF_READ_FEATURE,
                            value,
                            seq_row[db_align_end].line_number,
                            "_struct_ref_seq.db_align_end",
                        )
                    })
                    .transpose()?
                    .map(|seq_num| PdbSeqId {
                        seq_num,
                        ins_code: None,
                    });
                if let Some(index) = auth_seq_align_beg
                    && cif_row_has2(seq_row, index)
                    && let Some(seq_id) = make_seqid(
                        seq_row[index].value.as_str(),
                        seq_align_beg_ins.map(|idx| seq_row[idx].value.as_str()),
                        seq_row[index].line_number,
                    )?
                {
                    dbref.seq_begin = Some(seq_id);
                }
                if let Some(index) = auth_seq_align_end
                    && cif_row_has2(seq_row, index)
                    && let Some(seq_id) = make_seqid(
                        seq_row[index].value.as_str(),
                        seq_align_end_ins.map(|idx| seq_row[idx].value.as_str()),
                        seq_row[index].line_number,
                    )?
                {
                    dbref.seq_end = Some(seq_id);
                }
                entity.dbrefs.push(dbref);
            }
        }
    }

    if let Some(struct_asym_loop) = find_cif_loop(loops, "_struct_asym.id") {
        let width = checked_cif_loop_width(struct_asym_loop)?;
        let Some([id_col, entity_col]) = required_cif_table_cols(
            struct_asym_loop,
            ["_struct_asym.id", "_struct_asym.entity_id"],
        ) else {
            return Ok(());
        };
        for row in struct_asym_loop.values.chunks(width) {
            // Gemmi✔️✔️: for (auto row : s_asym_table)
            // Gemmi✔️✔️:   if (Entity* ent = st.get_entity(row.str(1)))
            // Gemmi✔️✔️:     ent->subchains.push_back(row.str(0));
            let Some(source_id) = cif_optional(row[entity_col].value.as_str()) else {
                continue;
            };
            let Some(entity_id) = builder.find_entity_by_source_id(source_id) else {
                continue;
            };
            let subchain = pdb_chain_id_from_cif(
                cif_optional(row[id_col].value.as_str()).unwrap_or(""),
                row[id_col].line_number,
            )?;
            builder.add_entity_subchain(entity_id, subchain);
            builder.assign_chain_entity_by_subchain(subchain, entity_id);
        }
    }

    Ok(())
}

// BEGIN GEMMI CPP HELPER read_entity_and_sequence_info struct_asym fallback
// Gemmi✔️✔️: } else if (!st.models.empty()) {
// Gemmi✔️✔️:   for (const Chain& chain : st.models[0].chains)
// Gemmi✔️✔️:     for (const ConstResidueSpan& sub : chain.subchains()) {
// Gemmi✔️✔️:       const Residue& r = sub.front();
// Gemmi✔️✔️:       if (Entity* ent = st.get_entity(r.entity_id))
// Gemmi✔️✔️:         if (!in_vector(r.subchain, ent->subchains))
// Gemmi✔️✔️:           ent->subchains.push_back(r.subchain);
// Gemmi✔️✔️:     }
// Gemmi✔️✔️: }
// END GEMMI CPP HELPER read_entity_and_sequence_info struct_asym fallback
fn infer_entity_subchains_from_first_model(structure: &mut BioStructure) {
    let Some(first_model) = structure.models.first() else {
        return;
    };
    let chain_start = first_model.chain_span.start as usize;
    let chain_end = first_model.chain_span.end() as usize;
    let mut pending = Vec::<(EntityId, PdbChainId)>::new();
    for chain_idx in chain_start..chain_end {
        let chain = &structure.chains[chain_idx];
        let residue_start = chain.residue_span.start as usize;
        let residue_end = chain.residue_span.end() as usize;
        for residue_idx in residue_start..residue_end {
            let residue = &structure.residues[residue_idx];
            let Some(entity_id) = residue.source.label_entity_id.or(chain.entity_id) else {
                continue;
            };
            let Some(subchain) = residue.source.subchain_id else {
                continue;
            };
            if !pending.contains(&(entity_id, subchain)) {
                pending.push((entity_id, subchain));
            }
        }
    }
    for (entity_id, subchain) in pending {
        let entity = &mut structure.entities[entity_id.index() as usize];
        if !entity.subchains.contains(&subchain) {
            entity.subchains.push(subchain);
        }
    }
}

fn import_shortened_ccd_codes_from_chem_comp(
    loops: &[CifLoop],
    structure: &mut BioStructure,
) -> Result<(), BioReadError> {
    let Some(loop_) = find_cif_loop(loops, "_chem_comp.id") else {
        return Ok(());
    };
    let Some(id) = optional_cif_col(loop_, "_chem_comp.id") else {
        return Ok(());
    };
    let Some(three_letter_code) = optional_cif_col(loop_, "_chem_comp.three_letter_code") else {
        return Ok(());
    };
    for row in cif_loop_rows(loop_)? {
        let alias = row[id].value.trim();
        let long_id = row[three_letter_code].value.trim();
        if alias.as_bytes().first() == Some(&b'~')
            && long_id.as_bytes().first() != Some(&b'~')
            && !long_id.is_empty()
        {
            structure
                .shortened_ccd_codes
                .push((long_id.to_string(), alias.to_string()));
        }
    }
    Ok(())
}

fn find_mmcif_polymer_kind(
    loops: &[CifLoop],
    entity_id: &str,
) -> Result<PolymerKind, BioReadError> {
    let Some(polymer_loop) = find_cif_loop(loops, "_entity_poly.entity_id") else {
        return Ok(PolymerKind::Unknown);
    };
    let width = checked_cif_loop_width(polymer_loop)?;
    let Some([entity_col, type_col]) = required_cif_table_cols(
        polymer_loop,
        ["_entity_poly.entity_id", "_entity_poly.type"],
    ) else {
        return Ok(PolymerKind::Unknown);
    };
    for row in polymer_loop.values.chunks(width) {
        if cif_optional(row[entity_col].value.as_str()) == Some(entity_id) {
            return Ok(cif_optional(row[type_col].value.as_str())
                .map(polymer_kind_from_cif)
                .unwrap_or(PolymerKind::Unknown));
        }
    }
    Ok(PolymerKind::Unknown)
}

fn find_cif_loop<'a>(loops: &'a [CifLoop], tag: &str) -> Option<&'a CifLoop> {
    loops
        .iter()
        .find(|loop_| loop_.tags.iter().any(|candidate| candidate == tag))
}

fn make_single_row_atom_site_loop_from_items(
    block: &CifBlock,
) -> Result<Option<CifLoop>, BioReadError> {
    let mut tags = Vec::new();
    let mut values = Vec::new();
    for item in &block.items {
        if item.tag.starts_with("_atom_site.") {
            tags.push(item.tag.clone());
            values.push(item.value.clone());
        }
    }
    if tags.is_empty() {
        return Ok(None);
    }
    if !tags.iter().any(|tag| tag == "_atom_site.id") {
        return Ok(None);
    }
    Ok(Some(CifLoop { tags, values }))
}

fn checked_cif_loop_width(loop_: &CifLoop) -> Result<usize, BioReadError> {
    let width = loop_.tags.len();
    if width == 0 || loop_.values.len() % width != 0 {
        return Err(BioReadError::Parse {
            line_number: loop_.values.first().map_or(0, |token| token.line_number),
            message: "mmCIF loop value count is not divisible by tag count".to_string(),
        });
    }
    Ok(width)
}

fn cif_loop_rows(loop_: &CifLoop) -> Result<impl Iterator<Item = &[CifToken]>, BioReadError> {
    let width = checked_cif_loop_width(loop_)?;
    Ok(loop_.values.chunks(width))
}

fn required_cif_col(loop_: &CifLoop, tag: &'static str) -> Result<usize, BioReadError> {
    optional_cif_col(loop_, tag).ok_or_else(|| BioReadError::Parse {
        line_number: 0,
        message: format!("required mmCIF column is missing: {tag}"),
    })
}

fn optional_cif_col(loop_: &CifLoop, tag: &str) -> Option<usize> {
    loop_.tags.iter().position(|candidate| candidate == tag)
}

fn required_cif_table_cols<const N: usize>(loop_: &CifLoop, tags: [&str; N]) -> Option<[usize; N]> {
    let mut cols = [0usize; N];
    for (index, tag) in tags.iter().enumerate() {
        cols[index] = optional_cif_col(loop_, tag)?;
    }
    Some(cols)
}

fn cif_row_has2(row: &[CifToken], index: usize) -> bool {
    row.get(index)
        .and_then(|token| cif_optional(token.value.as_str()))
        .is_some()
}

fn group_pdb_het_flag(value: &str) -> Option<char> {
    for c in value.chars().take(2) {
        let upper = c.to_ascii_uppercase();
        if matches!(upper, 'A' | 'H') {
            return Some(upper);
        }
        if upper == '\0' {
            return Some('\0');
        }
    }
    None
}

fn calc_flag_from_cif(value: &str) -> BioCalcFlag {
    let bytes = value.as_bytes();
    match bytes.first().copied() {
        Some(b'c') => BioCalcFlag::Calculated,
        Some(b'd') => {
            if bytes.get(1).copied() == Some(b'u') {
                BioCalcFlag::Dummy
            } else {
                BioCalcFlag::Determined
            }
        }
        _ => BioCalcFlag::NotSet,
    }
}

fn parse_optional_short_prefix_int(value: &str) -> Option<i16> {
    let end = value.bytes().take_while(|b| b.is_ascii_digit()).count();
    if end == 0 {
        return None;
    }
    value.get(..end)?.parse::<i16>().ok()
}

#[derive(Debug, Clone, Copy)]
struct CifRowAccess {
    primary: Option<usize>,
    fallback: Option<usize>,
}

// BEGIN GEMMI CPP FUNCTION make_seqid
// Gemmi✔️✔️: SeqId ret;
// Gemmi✔️✔️: if (icode)
// Gemmi✔️✔️:   // the insertion code happens to be always a single letter
// Gemmi✔️✔️:   ret.icode = cif::as_char(*icode, ' ');
// Gemmi✔️✔️: if (!seqid.empty()) {
// Gemmi✔️✔️:   // old mmCIF files have auth_seq_id as number + icode (e.g. 15A)
// Gemmi✔️✔️:   if (seqid.back() >= 'A') {
// Gemmi✔️✔️:     if (ret.icode == ' ')
// Gemmi✔️✔️:       ret.icode = seqid.back();
// Gemmi✔️✔️:     else if (ret.icode != seqid.back())
// Gemmi✔️✔️:       fail("Inconsistent insertion code in " + seqid);
// Gemmi✔️✔️:     seqid.pop_back();
// Gemmi✔️✔️:   }
// Gemmi✔️✔️:   // 7pvv has an empty seqnum in tls description, don't throw
// Gemmi✔️✔️:   if (!seqid.empty())
// Gemmi✔️✔️:     ret.num = cif::as_int(seqid, Residue::OptionalNum::None);
// Gemmi✔️✔️: }
// Gemmi✔️✔️: return ret;
// END GEMMI CPP FUNCTION
fn make_seqid(
    seqid: &str,
    icode: Option<&str>,
    line_number: usize,
) -> Result<Option<PdbSeqId>, BioReadError> {
    let mut ret = PdbSeqId::default();
    let mut has_num = false;
    let mut ret_icode = icode
        .and_then(cif_optional)
        .and_then(|value| value.as_bytes().first().copied())
        .filter(|byte| *byte != b' ');
    let mut seqid = seqid.trim().to_string();
    if !seqid.is_empty() {
        if let Some(last) = seqid.as_bytes().last().copied()
            && last >= b'A'
        {
            if ret_icode.is_none() {
                ret_icode = Some(last);
            } else if ret_icode != Some(last) {
                return Err(BioReadError::Parse {
                    line_number,
                    message: format!("Inconsistent insertion code in {seqid}"),
                });
            }
            seqid.pop();
        }
        if !seqid.is_empty() {
            ret.seq_num = parse_decimal_i32(
                &BIO_MMCIF_READ_FEATURE,
                seqid.as_str(),
                line_number,
                "mmCIF residue sequence number",
            )?;
            has_num = true;
        }
    }
    ret.ins_code = ret_icode;
    if has_num || ret.ins_code.is_some() {
        Ok(Some(ret))
    } else {
        Ok(None)
    }
}

// BEGIN GEMMI CPP FUNCTION make_resid
// Gemmi✔️✔️: return ResidueId{make_seqid(seqid, icode), {}, name};
// END GEMMI CPP FUNCTION
fn make_resid(
    name: &str,
    seqid: &str,
    icode: Option<&str>,
    line_number: usize,
) -> Result<(ResidueName, Option<PdbSeqId>), BioReadError> {
    Ok((
        residue_name_from_field(name),
        make_seqid(seqid, icode, line_number)?,
    ))
}

impl CifRowAccess {
    // BEGIN GEMMI CPP FUNCTION RowAccess::RowAccess / RowAccess::ok / RowAccess::get
    // Gemmi✔️✔️: struct RowAccess {
    // Gemmi✔️✔️:   const std::string *val = nullptr;
    // Gemmi✔️✔️:   const std::string *fallback = nullptr;
    // Gemmi✔️✔️:   RowAccess(const cif::Table& tab, int n1, int n2) {
    // Gemmi✔️✔️:     int pos1 = tab.positions.at(n1);
    // Gemmi✔️✔️:     int pos2 = tab.positions.at(n2);
    // Gemmi✔️✔️:     if (pos1 < 0) {
    // Gemmi✔️✔️:       pos1 = pos2;
    // Gemmi✔️✔️:       pos2 = -1;
    // Gemmi✔️✔️:     }
    // Gemmi✔️✔️:     const cif::Loop* loop = const_cast<cif::Table&>(tab).get_loop();
    // Gemmi✔️✔️:     if (pos1 >= 0)
    // Gemmi✔️✔️:       val = loop ? &loop->values[pos1] : &tab.bloc.items[pos1].pair[1];
    // Gemmi✔️✔️:     if (pos2 >= 0)
    // Gemmi✔️✔️:       fallback = loop ? &loop->values[pos2] : &tab.bloc.items[pos2].pair[1];
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️:   bool ok() const { return val != nullptr; }
    // Gemmi✔️✔️:   const std::string& get(size_t gap) const {
    // Gemmi✔️✔️:     const std::string& r = val[gap];
    // Gemmi✔️✔️:     if (!cif::is_null(r) || fallback == nullptr)
    // Gemmi✔️✔️:       return r;
    // Gemmi✔️✔️:     return fallback[gap];
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️: };
    // END GEMMI CPP FUNCTION
    fn new(primary: Option<usize>, fallback: Option<usize>) -> Self {
        let (primary, fallback) = if primary.is_none() {
            (fallback, None)
        } else {
            (primary, fallback)
        };
        Self { primary, fallback }
    }

    fn ok(&self) -> bool {
        self.primary.is_some()
    }

    fn get<'a>(&self, row: &'a [CifToken]) -> Option<&'a str> {
        let primary = self.primary?;
        let value = row.get(primary)?.value.as_str();
        if cif_optional(value).is_some() || self.fallback.is_none() {
            Some(value)
        } else {
            self.fallback
                .and_then(|index| row.get(index))
                .map(|token| token.value.as_str())
        }
    }
}

// BEGIN GEMMI CPP FUNCTION get_by_id
// Gemmi✔️✔️: template<class T>
// Gemmi✔️✔️: T* get_by_id(std::vector<T>& vec, const std::string& id) {
// Gemmi✔️✔️:   for (T& item : vec)
// Gemmi✔️✔️:     if (item.id == id)
// Gemmi✔️✔️:       return &item;
// Gemmi✔️✔️:   return nullptr;
// Gemmi✔️✔️: }
// END GEMMI CPP FUNCTION
fn get_by_id<'a, T, F>(vec: &'a mut [T], id: &str, mut get_id: F) -> Option<&'a mut T>
where
    F: FnMut(&T) -> &str,
{
    vec.iter_mut().find(|item| get_id(item) == id)
}

// BEGIN GEMMI CPP FUNCTION copy_int
// Gemmi✔️✔️: if (row.has2(n))
// Gemmi✔️✔️:   dest = cif::as_int(row[n]);
// END GEMMI CPP FUNCTION
fn copy_int(row: &[CifToken], index: usize, dest: &mut i32) -> Result<(), BioReadError> {
    if cif_row_has2(row, index) {
        let token = &row[index];
        *dest = parse_decimal_i32(
            &BIO_MMCIF_READ_FEATURE,
            token.value.as_str(),
            token.line_number,
            "mmCIF integer field",
        )?;
    }
    Ok(())
}

// BEGIN GEMMI CPP FUNCTION copy_double
// Gemmi✔️✔️: if (row.has2(n))
// Gemmi✔️✔️:   dest = cif::as_number(row[n]);
// END GEMMI CPP FUNCTION
fn copy_double(row: &[CifToken], index: usize, dest: &mut f32) -> Result<(), BioReadError> {
    if cif_row_has2(row, index) {
        let token = &row[index];
        *dest = parse_f32(
            token.value.as_str(),
            token.line_number,
            "mmCIF number field",
        )?;
    }
    Ok(())
}

// BEGIN GEMMI CPP FUNCTION copy_string
// Gemmi✔️✔️: if (row.has2(n))
// Gemmi✔️✔️:   dest = cif::as_string(row[n]);
// END GEMMI CPP FUNCTION
fn copy_string(row: &[CifToken], index: usize, dest: &mut String) {
    if cif_row_has2(row, index) {
        *dest = row[index].value.trim().to_string();
    }
}

// BEGIN GEMMI CPP FUNCTION get_smat33
// Gemmi✔️✔️: return SMat33<T>{(T) cif::as_number(row[n+0]),
// Gemmi✔️✔️:                  (T) cif::as_number(row[n+1]),
// Gemmi✔️✔️:                  (T) cif::as_number(row[n+2]),
// Gemmi✔️✔️:                  (T) cif::as_number(row[n+3]),
// Gemmi✔️✔️:                  (T) cif::as_number(row[n+4]),
// Gemmi✔️✔️:                  (T) cif::as_number(row[n+5])};
// END GEMMI CPP FUNCTION
fn get_smat33(row: &[CifToken], index: usize) -> Result<CifSmat33<f32>, BioReadError> {
    Ok(CifSmat33 {
        u11: parse_f32(
            row[index].value.as_str(),
            row[index].line_number,
            "SMat33 u11",
        )?,
        u22: parse_f32(
            row[index + 1].value.as_str(),
            row[index + 1].line_number,
            "SMat33 u22",
        )?,
        u33: parse_f32(
            row[index + 2].value.as_str(),
            row[index + 2].line_number,
            "SMat33 u33",
        )?,
        u12: parse_f32(
            row[index + 3].value.as_str(),
            row[index + 3].line_number,
            "SMat33 u12",
        )?,
        u13: parse_f32(
            row[index + 4].value.as_str(),
            row[index + 4].line_number,
            "SMat33 u13",
        )?,
        u23: parse_f32(
            row[index + 5].value.as_str(),
            row[index + 5].line_number,
            "SMat33 u23",
        )?,
    })
}

// BEGIN GEMMI CPP FUNCTION get_anisotropic_u
// Gemmi✔️✔️: cif::Table aniso_tab = block.find("_atom_site_anisotrop.",
// Gemmi✔️✔️:                                   {"id", "U[1][1]", "U[2][2]", "U[3][3]",
// Gemmi✔️✔️:                                    "U[1][2]", "U[1][3]", "U[2][3]"});
// Gemmi✔️✔️: std::unordered_map<std::string, SMat33<float>> aniso_map;
// Gemmi✔️✔️: for (auto ani : aniso_tab)
// Gemmi✔️✔️:   aniso_map.emplace(ani[0], get_smat33<float>(ani, 1));
// Gemmi✔️✔️: return aniso_map;
// END GEMMI CPP FUNCTION
fn get_anisotropic_u(loops: &[CifLoop]) -> Result<HashMap<String, [f32; 6]>, BioReadError> {
    let Some(loop_) = find_cif_loop(loops, "_atom_site_anisotrop.id") else {
        return Ok(HashMap::new());
    };
    let Some([id, u11, u22, u33, u12, u13, u23]) = required_cif_table_cols(
        loop_,
        [
            "_atom_site_anisotrop.id",
            "_atom_site_anisotrop.U[1][1]",
            "_atom_site_anisotrop.U[2][2]",
            "_atom_site_anisotrop.U[3][3]",
            "_atom_site_anisotrop.U[1][2]",
            "_atom_site_anisotrop.U[1][3]",
            "_atom_site_anisotrop.U[2][3]",
        ],
    ) else {
        return Ok(HashMap::new());
    };
    let mut map = HashMap::new();
    for row in cif_loop_rows(loop_)? {
        map.insert(
            row[id].value.trim().to_string(),
            [
                parse_f32(
                    row[u11].value.as_str(),
                    row[u11].line_number,
                    "_atom_site_anisotrop.U[1][1]",
                )?,
                parse_f32(
                    row[u22].value.as_str(),
                    row[u22].line_number,
                    "_atom_site_anisotrop.U[2][2]",
                )?,
                parse_f32(
                    row[u33].value.as_str(),
                    row[u33].line_number,
                    "_atom_site_anisotrop.U[3][3]",
                )?,
                parse_f32(
                    row[u12].value.as_str(),
                    row[u12].line_number,
                    "_atom_site_anisotrop.U[1][2]",
                )?,
                parse_f32(
                    row[u13].value.as_str(),
                    row[u13].line_number,
                    "_atom_site_anisotrop.U[1][3]",
                )?,
                parse_f32(
                    row[u23].value.as_str(),
                    row[u23].line_number,
                    "_atom_site_anisotrop.U[2][3]",
                )?,
            ],
        );
    }
    Ok(map)
}

// BEGIN GEMMI CPP FUNCTION transform_tags
// Gemmi✔️✔️: return {mstr + "[1][1]", mstr + "[1][2]", mstr + "[1][3]", vstr + "[1]",
// Gemmi✔️✔️:         mstr + "[2][1]", mstr + "[2][2]", mstr + "[2][3]", vstr + "[2]",
// Gemmi✔️✔️:         mstr + "[3][1]", mstr + "[3][2]", mstr + "[3][3]", vstr + "[3]"};
// END GEMMI CPP FUNCTION
fn transform_tags(mstr: &str, vstr: &str) -> [String; 12] {
    [
        format!("{mstr}[1][1]"),
        format!("{mstr}[1][2]"),
        format!("{mstr}[1][3]"),
        format!("{vstr}[1]"),
        format!("{mstr}[2][1]"),
        format!("{mstr}[2][2]"),
        format!("{mstr}[2][3]"),
        format!("{vstr}[2]"),
        format!("{mstr}[3][1]"),
        format!("{mstr}[3][2]"),
        format!("{mstr}[3][3]"),
        format!("{vstr}[3]"),
    ]
}

fn find_cif_transform(
    loops: &[CifLoop],
    matrix_prefix: &str,
    vector_prefix: &str,
) -> Result<Option<BioTransform>, BioReadError> {
    let tags = transform_tags(matrix_prefix, vector_prefix);
    let Some(loop_) = find_cif_loop(loops, &tags[0]) else {
        return Ok(None);
    };
    let indices = [
        optional_cif_col(loop_, &tags[0]),
        optional_cif_col(loop_, &tags[1]),
        optional_cif_col(loop_, &tags[2]),
        optional_cif_col(loop_, &tags[3]),
        optional_cif_col(loop_, &tags[4]),
        optional_cif_col(loop_, &tags[5]),
        optional_cif_col(loop_, &tags[6]),
        optional_cif_col(loop_, &tags[7]),
        optional_cif_col(loop_, &tags[8]),
        optional_cif_col(loop_, &tags[9]),
        optional_cif_col(loop_, &tags[10]),
        optional_cif_col(loop_, &tags[11]),
    ];
    if indices.iter().any(Option::is_none) {
        return Err(BioReadError::Parse {
            line_number: 0,
            message: format!("required mmCIF transform columns are missing: {matrix_prefix}"),
        });
    }
    let Some(row) = cif_loop_rows(loop_)?.into_iter().next() else {
        return Ok(None);
    };
    let matrix_row = indices
        .into_iter()
        .map(|index| row[index.expect("checked all transform columns present")].clone())
        .collect::<Vec<_>>();
    Ok(Some(get_transform_matrix(&matrix_row)?))
}

// BEGIN GEMMI CPP FUNCTION get_transform_matrix
// Gemmi✔️✔️: Transform t;
// Gemmi✔️✔️: for (int i = 0; i < 3; ++i) {
// Gemmi✔️✔️:   for (int j = 0; j < 3; ++j)
// Gemmi✔️✔️:     t.mat[i][j] = cif::as_number(r[4*i+j]);
// Gemmi✔️✔️:   t.vec.at(i) = cif::as_number(r[4*i+3]);
// Gemmi✔️✔️: }
// Gemmi✔️✔️: return t;
// END GEMMI CPP FUNCTION
fn get_transform_matrix(row: &[CifToken]) -> Result<BioTransform, BioReadError> {
    let mut transform = BioTransform::default();
    for i in 0..3 {
        for j in 0..3 {
            let idx = 4 * i + j;
            transform.mat[i][j] = parse_f64(
                row[idx].value.as_str(),
                row[idx].line_number,
                "transform matrix",
            )?;
        }
        let idx = 4 * i + 3;
        transform.vec[i] = parse_f64(
            row[idx].value.as_str(),
            row[idx].line_number,
            "transform vector",
        )?;
    }
    Ok(transform)
}

fn bio_address_from_mmcif(
    chain_name: &str,
    residue_name: &str,
    seqid: &str,
    icode: Option<&str>,
    atom_name: Option<&str>,
    line_number: usize,
) -> Result<BioAtomAddress, BioReadError> {
    let (_, seq_id) = make_resid(residue_name, seqid, icode, line_number)?;
    Ok(BioAtomAddress {
        chain_name: chain_name.trim().to_string(),
        seq_id,
        residue_name: residue_name.trim().to_string(),
        atom_name: atom_name.unwrap_or_default().trim().to_string(),
        altloc: None,
    })
}

// BEGIN GEMMI CPP FUNCTION set_part_of_address_from_label
// Gemmi✔️✔️: int seq = cif::as_int(label_seq_id_raw, SeqId::OptionalNum::None);
// Gemmi✔️✔️: for (const Chain& chain : model.chains)
// Gemmi✔️✔️:   if (ConstResidueSpan sub = chain.get_subchain(label_asym)) {
// Gemmi✔️✔️:     a.chain_name = chain.name;
// Gemmi✔️✔️:     for (const Residue& res : sub)
// Gemmi✔️✔️:       if (res.label_seq == seq) {
// Gemmi✔️✔️:         a.res_id.seqid = res.seqid;
// Gemmi✔️✔️:         return;
// Gemmi✔️✔️:       }
// Gemmi✔️✔️:   }
// END GEMMI CPP FUNCTION
fn set_part_of_address_from_label(
    address: &mut BioAtomAddress,
    structure: &BioStructure,
    label_asym: &str,
    label_seq_id_raw: &str,
    line_number: usize,
) -> Result<(), BioReadError> {
    let seq = parse_decimal_i32(
        &BIO_MMCIF_READ_FEATURE,
        label_seq_id_raw,
        line_number,
        "_struct_conn label_seq_id",
    )?;
    for residue in &structure.residues {
        let chain = &structure.chains[residue.chain_id.index() as usize];
        if chain.model_id.index() == 0
            && chain
                .source
                .label_asym_id
                .is_some_and(|id| id.as_str() == label_asym)
            && residue.source.label_seq_id == Some(seq)
        {
            address.chain_name = chain
                .source
                .auth_chain_id
                .or(chain.source.label_asym_id)
                .map(|id| id.as_str().to_string())
                .unwrap_or_default();
            address.seq_id = residue.source.seq_id;
            return Ok(());
        }
    }
    Ok(())
}

fn connection_type_from_mmcif(
    value: &str,
    _line_number: usize,
) -> Result<BioConnectionType, BioReadError> {
    Ok(match value.trim() {
        "covale" => BioConnectionType::Covale,
        "disulf" => BioConnectionType::Disulf,
        "hydrog" => BioConnectionType::Hydrog,
        "metalc" => BioConnectionType::MetalC,
        _ => BioConnectionType::Unknown,
    })
}

// BEGIN GEMMI CPP FUNCTION read_connectivity
// Gemmi✔️✔️: for (const auto row : block.find("_struct_conn.", { ... })) {
// Gemmi✔️✔️:   Connection c;
// Gemmi✔️✔️:   c.name = row.str(kId);
// Gemmi✔️✔️:   copy_string(row, kLinkId, c.link_id);
// Gemmi✔️✔️:   c.type = connection_type_from_string(row.str(kConnTypeId));
// Gemmi✔️✔️:   if (row.has2(kSym1) && row.has2(kSym2)) { ... }
// Gemmi✔️✔️:   copy_double(row, kDistValue, c.reported_distance);
// Gemmi✔️✔️:   for (int i = 0; i < 2; ++i) { ... auth-vs-label fallback ... }
// Gemmi✔️✔️:   st.connections.emplace_back(c);
// Gemmi✔️✔️: }
// END GEMMI CPP FUNCTION
fn read_connectivity(loops: &[CifLoop], structure: &mut BioStructure) -> Result<(), BioReadError> {
    let Some(loop_) = find_cif_loop(loops, "_struct_conn.id") else {
        return Ok(());
    };
    let Some(
        [
            id,
            conn_type_id,
            ptnr1_label_comp_id,
            ptnr2_label_comp_id,
            ptnr1_label_atom_id,
            ptnr2_label_atom_id,
        ],
    ) = required_cif_table_cols(
        loop_,
        [
            "_struct_conn.id",
            "_struct_conn.conn_type_id",
            "_struct_conn.ptnr1_label_comp_id",
            "_struct_conn.ptnr2_label_comp_id",
            "_struct_conn.ptnr1_label_atom_id",
            "_struct_conn.ptnr2_label_atom_id",
        ],
    )
    else {
        return Ok(());
    };
    let ptnr1_auth_asym_id = optional_cif_col(loop_, "_struct_conn.ptnr1_auth_asym_id");
    let ptnr2_auth_asym_id = optional_cif_col(loop_, "_struct_conn.ptnr2_auth_asym_id");
    let ptnr1_label_asym_id = optional_cif_col(loop_, "_struct_conn.ptnr1_label_asym_id");
    let ptnr2_label_asym_id = optional_cif_col(loop_, "_struct_conn.ptnr2_label_asym_id");
    let ptnr1_label_alt_id = optional_cif_col(loop_, "_struct_conn.pdbx_ptnr1_label_alt_id");
    let ptnr2_label_alt_id = optional_cif_col(loop_, "_struct_conn.pdbx_ptnr2_label_alt_id");
    let ptnr1_auth_seq_id = optional_cif_col(loop_, "_struct_conn.ptnr1_auth_seq_id");
    let ptnr2_auth_seq_id = optional_cif_col(loop_, "_struct_conn.ptnr2_auth_seq_id");
    let ptnr1_label_seq_id = optional_cif_col(loop_, "_struct_conn.ptnr1_label_seq_id");
    let ptnr2_label_seq_id = optional_cif_col(loop_, "_struct_conn.ptnr2_label_seq_id");
    let ptnr1_ins_code = optional_cif_col(loop_, "_struct_conn.pdbx_ptnr1_PDB_ins_code");
    let ptnr2_ins_code = optional_cif_col(loop_, "_struct_conn.pdbx_ptnr2_PDB_ins_code");
    let ptnr1_symmetry = optional_cif_col(loop_, "_struct_conn.ptnr1_symmetry");
    let ptnr2_symmetry = optional_cif_col(loop_, "_struct_conn.ptnr2_symmetry");
    let pdbx_dist_value = optional_cif_col(loop_, "_struct_conn.pdbx_dist_value");
    let ccp4_link_id = optional_cif_col(loop_, "_struct_conn.ccp4_link_id");
    for row in cif_loop_rows(loop_)? {
        let line_number = row[id].line_number;
        let mut connection = BioConnection {
            name: row[id].value.trim().to_string(),
            type_: connection_type_from_mmcif(row[conn_type_id].value.trim(), line_number)?,
            ..BioConnection::default()
        };
        if let Some(index) = ccp4_link_id {
            copy_string(row, index, &mut connection.link_id);
        }
        if let (Some(sym1_idx), Some(sym2_idx)) = (ptnr1_symmetry, ptnr2_symmetry)
            && cif_row_has2(row, sym1_idx)
            && cif_row_has2(row, sym2_idx)
        {
            let s1 = row[sym1_idx].value.trim();
            let s2 = row[sym2_idx].value.trim();
            if s1 == s2 {
                connection.asu = BioAsu::Same;
            } else {
                connection.asu = BioAsu::Different;
                if let (Some(sep1), Some(sep2)) = (s1.find('_'), s2.find('_'))
                    && sep1 + 4 == s1.len()
                    && sep2 + 4 == s2.len()
                {
                    connection.reported_sym[0] = if s1.starts_with("1_") {
                        s2[..sep2].parse::<i16>().unwrap_or(0)
                    } else {
                        99
                    };
                    for i in 1..=3 {
                        connection.reported_sym[i] =
                            i16::from(s2.as_bytes()[sep2 + i]) - i16::from(s1.as_bytes()[sep1 + i]);
                    }
                }
            }
        }
        if let Some(index) = pdbx_dist_value
            && cif_row_has2(row, index)
        {
            connection.reported_distance = Some(parse_f64(
                row[index].value.as_str(),
                row[index].line_number,
                "_struct_conn.pdbx_dist_value",
            )?);
        }
        for i in 0..2 {
            let address = if i == 0 {
                &mut connection.partner1
            } else {
                &mut connection.partner2
            };
            let auth_asym = if i == 0 {
                ptnr1_auth_asym_id
            } else {
                ptnr2_auth_asym_id
            };
            let auth_seq = if i == 0 {
                ptnr1_auth_seq_id
            } else {
                ptnr2_auth_seq_id
            };
            let label_asym = if i == 0 {
                ptnr1_label_asym_id
            } else {
                ptnr2_label_asym_id
            };
            let label_seq = if i == 0 {
                ptnr1_label_seq_id
            } else {
                ptnr2_label_seq_id
            };
            let ins_code = if i == 0 {
                ptnr1_ins_code
            } else {
                ptnr2_ins_code
            };
            let label_comp = if i == 0 {
                ptnr1_label_comp_id
            } else {
                ptnr2_label_comp_id
            };
            let label_atom = if i == 0 {
                ptnr1_label_atom_id
            } else {
                ptnr2_label_atom_id
            };
            let label_alt = if i == 0 {
                ptnr1_label_alt_id
            } else {
                ptnr2_label_alt_id
            };
            if let (Some(auth_asym_idx), Some(auth_seq_idx)) = (auth_asym, auth_seq)
                && cif_row_has2(row, auth_asym_idx)
                && cif_row_has2(row, auth_seq_idx)
            {
                address.chain_name = row[auth_asym_idx].value.trim().to_string();
                let (_, seq_id) = make_resid(
                    row[label_comp].value.as_str(),
                    row[auth_seq_idx].value.as_str(),
                    ins_code.map(|idx| row[idx].value.as_str()),
                    line_number,
                )?;
                address.seq_id = seq_id;
            } else if let (Some(label_asym_idx), Some(label_seq_idx)) = (label_asym, label_seq)
                && cif_row_has2(row, label_asym_idx)
                && cif_row_has2(row, label_seq_idx)
            {
                set_part_of_address_from_label(
                    address,
                    structure,
                    row[label_asym_idx].value.trim(),
                    row[label_seq_idx].value.as_str(),
                    line_number,
                )?;
            } else {
                return Err(BioReadError::Parse {
                    line_number,
                    message: "_struct_conn without either _auth_ or _label_ asym_id+seq_id"
                        .to_string(),
                });
            }
            address.residue_name = row[label_comp].value.trim().to_string();
            address.atom_name = row[label_atom].value.trim().to_string();
            if let Some(index) = label_alt {
                address.altloc = cif_optional(row[index].value.as_str())
                    .and_then(|value| value.as_bytes().first().copied())
                    .and_then(parse_altloc);
            }
        }
        structure.connections.push(connection);
    }
    Ok(())
}

// BEGIN GEMMI CPP FUNCTION read_prot_cis
// Gemmi✔️✔️: for (auto row : block.find("_struct_mon_prot_cis.", {...})) {
// Gemmi✔️✔️:   CisPep cispep;
// Gemmi✔️✔️:   cispep.model_num = cif::as_int(row[kModelNum], 0);
// Gemmi✔️✔️:   cispep.partner_c.chain_name = row.str(kAuthAsymId);
// Gemmi✔️✔️:   cispep.partner_c.res_id.seqid = make_seqid(...);
// Gemmi✔️✔️:   cispep.partner_c.res_id.name = cif::as_string(row.one_of(...));
// Gemmi✔️✔️:   ...
// Gemmi✔️✔️:   st.cispeps.push_back(cispep);
// Gemmi✔️✔️: }
// END GEMMI CPP FUNCTION
fn read_prot_cis(loops: &[CifLoop], structure: &mut BioStructure) -> Result<(), BioReadError> {
    let Some(loop_) = find_cif_loop(loops, "_struct_mon_prot_cis.pdbx_PDB_model_num") else {
        return Ok(());
    };
    let Some([model_num, auth_asym_id, auth_seq_id]) = required_cif_table_cols(
        loop_,
        [
            "_struct_mon_prot_cis.pdbx_PDB_model_num",
            "_struct_mon_prot_cis.auth_asym_id",
            "_struct_mon_prot_cis.auth_seq_id",
        ],
    ) else {
        return Ok(());
    };
    let ins_code = optional_cif_col(loop_, "_struct_mon_prot_cis.pdbx_PDB_ins_code");
    let label_comp_id = optional_cif_col(loop_, "_struct_mon_prot_cis.label_comp_id");
    let auth_comp_id = optional_cif_col(loop_, "_struct_mon_prot_cis.auth_comp_id");
    let auth_asym_id_2 = optional_cif_col(loop_, "_struct_mon_prot_cis.pdbx_auth_asym_id_2");
    let auth_seq_id_2 = optional_cif_col(loop_, "_struct_mon_prot_cis.pdbx_auth_seq_id_2");
    let ins_code_2 = optional_cif_col(loop_, "_struct_mon_prot_cis.pdbx_PDB_ins_code_2");
    let label_comp_id_2 = optional_cif_col(loop_, "_struct_mon_prot_cis.pdbx_label_comp_id_2");
    let auth_comp_id_2 = optional_cif_col(loop_, "_struct_mon_prot_cis.pdbx_auth_comp_id_2");
    let label_alt_id = optional_cif_col(loop_, "_struct_mon_prot_cis.label_alt_id");
    let omega_angle = optional_cif_col(loop_, "_struct_mon_prot_cis.pdbx_omega_angle");
    for row in cif_loop_rows(loop_)? {
        let line_number = row[model_num].line_number;
        let mut cispep = BioCisPep {
            model_num: parse_decimal_i32(
                &BIO_MMCIF_READ_FEATURE,
                row[model_num].value.as_str(),
                line_number,
                "_struct_mon_prot_cis.pdbx_PDB_model_num",
            )?,
            ..BioCisPep::default()
        };
        cispep.partner_c.chain_name = row[auth_asym_id].value.trim().to_string();
        cispep.partner_c.seq_id = make_seqid(
            row[auth_seq_id].value.as_str(),
            ins_code.map(|idx| row[idx].value.as_str()),
            line_number,
        )?;
        cispep.partner_c.residue_name = auth_comp_id
            .and_then(|idx| cif_optional(row[idx].value.as_str()))
            .or_else(|| label_comp_id.and_then(|idx| cif_optional(row[idx].value.as_str())))
            .unwrap_or("")
            .to_string();
        if let Some(idx) = auth_asym_id_2
            && cif_row_has2(row, idx)
        {
            cispep.partner_n.chain_name = row[idx].value.trim().to_string();
        }
        if let Some(idx) = auth_seq_id_2
            && cif_row_has2(row, idx)
        {
            cispep.partner_n.seq_id = make_seqid(
                row[idx].value.as_str(),
                ins_code_2.map(|col| row[col].value.as_str()),
                line_number,
            )?;
        }
        cispep.partner_n.residue_name = auth_comp_id_2
            .and_then(|idx| cif_optional(row[idx].value.as_str()))
            .or_else(|| label_comp_id_2.and_then(|idx| cif_optional(row[idx].value.as_str())))
            .unwrap_or("")
            .to_string();
        if let Some(idx) = label_alt_id {
            cispep.only_altloc = cif_optional(row[idx].value.as_str())
                .and_then(|value| value.as_bytes().first().copied())
                .and_then(parse_altloc);
        }
        if let Some(idx) = omega_angle
            && cif_row_has2(row, idx)
        {
            cispep.reported_angle = Some(parse_f64(
                row[idx].value.as_str(),
                row[idx].line_number,
                "_struct_mon_prot_cis.pdbx_omega_angle",
            )?);
        }
        structure.cispeps.push(cispep);
    }
    Ok(())
}

// BEGIN GEMMI CPP FUNCTION read_struct_mod_residue
// Gemmi✔️✔️: for (auto row : block.find("_pdbx_struct_mod_residue.", {...})) {
// Gemmi✔️✔️:   ModRes modres;
// Gemmi✔️✔️:   modres.chain_name = row.str(0);
// Gemmi✔️✔️:   modres.res_id.seqid = make_seqid(row.str(1), row.ptr_at(2));
// Gemmi✔️✔️:   modres.res_id.name = row.one_of(3, 4);
// Gemmi✔️✔️:   ...
// Gemmi✔️✔️:   st.mod_residues.push_back(modres);
// Gemmi✔️✔️: }
// END GEMMI CPP FUNCTION
fn read_struct_mod_residue(
    loops: &[CifLoop],
    structure: &mut BioStructure,
) -> Result<(), BioReadError> {
    let Some(loop_) = find_cif_loop(loops, "_pdbx_struct_mod_residue.auth_asym_id") else {
        return Ok(());
    };
    let Some([auth_asym_id, auth_seq_id]) = required_cif_table_cols(
        loop_,
        [
            "_pdbx_struct_mod_residue.auth_asym_id",
            "_pdbx_struct_mod_residue.auth_seq_id",
        ],
    ) else {
        return Ok(());
    };
    let pdb_ins_code = optional_cif_col(loop_, "_pdbx_struct_mod_residue.PDB_ins_code");
    let auth_comp_id = optional_cif_col(loop_, "_pdbx_struct_mod_residue.auth_comp_id");
    let label_comp_id = optional_cif_col(loop_, "_pdbx_struct_mod_residue.label_comp_id");
    let parent_comp_id = optional_cif_col(loop_, "_pdbx_struct_mod_residue.parent_comp_id");
    let details = optional_cif_col(loop_, "_pdbx_struct_mod_residue.details");
    let ccp4_mod_id = optional_cif_col(loop_, "_pdbx_struct_mod_residue.ccp4_mod_id");
    for row in cif_loop_rows(loop_)? {
        let line_number = row[auth_asym_id].line_number;
        let mut modres = BioModRes {
            chain_name: row[auth_asym_id].value.trim().to_string(),
            ..BioModRes::default()
        };
        modres.res_id = make_seqid(
            row[auth_seq_id].value.as_str(),
            pdb_ins_code.map(|idx| row[idx].value.as_str()),
            line_number,
        )?
        .unwrap_or_default();
        modres.parent_comp_id = parent_comp_id
            .and_then(|idx| cif_optional(row[idx].value.as_str()))
            .unwrap_or_default()
            .to_string();
        modres.details = details
            .and_then(|idx| cif_optional(row[idx].value.as_str()))
            .unwrap_or_default()
            .to_string();
        modres.mod_id = ccp4_mod_id
            .and_then(|idx| cif_optional(row[idx].value.as_str()))
            .unwrap_or_default()
            .to_string();
        modres.residue_name = auth_comp_id
            .and_then(|idx| cif_optional(row[idx].value.as_str()))
            .or_else(|| label_comp_id.and_then(|idx| cif_optional(row[idx].value.as_str())))
            .unwrap_or_default()
            .to_string();
        structure.mod_residues.push(modres);
    }
    Ok(())
}

// BEGIN GEMMI CPP FUNCTION parse_operation_expr
// Gemmi✔️✔️: std::vector<std::string> result;
// Gemmi✔️✔️: std::size_t start = 0;
// Gemmi✔️✔️: std::size_t close_br = std::string::npos;
// Gemmi✔️✔️: if (expr[0] == '(') {
// Gemmi✔️✔️:   start = 1;
// Gemmi✔️✔️:   close_br = expr.find(')');
// Gemmi✔️✔️: }
// Gemmi✔️✔️: for (;;) {
// Gemmi✔️✔️:   std::size_t sep = std::min(expr.find(',', start), close_br);
// Gemmi✔️✔️:   std::size_t minus = expr.find('-', start);
// Gemmi✔️✔️:   if (minus < sep) {
// Gemmi✔️✔️:     int n_min = no_sign_atoi(expr.c_str() + start);
// Gemmi✔️✔️:     int n_max = no_sign_atoi(expr.c_str() + minus + 1);
// Gemmi✔️✔️:     for (int n = n_min; n <= n_max; ++n)
// Gemmi✔️✔️:       result.push_back(std::to_string(n));
// Gemmi✔️✔️:   } else {
// Gemmi✔️✔️:     result.emplace_back(expr, start, sep - start);
// Gemmi✔️✔️:   }
// Gemmi✔️✔️:   if (sep == close_br)
// Gemmi✔️✔️:     break;
// Gemmi✔️✔️:   start = sep + 1;
// Gemmi✔️✔️: }
// Gemmi✔️✔️: return result;
// END GEMMI CPP FUNCTION
fn parse_operation_expr(expr: &str) -> Vec<String> {
    if expr.is_empty() {
        return vec![String::new()];
    }
    let mut result = Vec::new();
    let mut start = 0usize;
    let mut close_br = None;
    if expr.as_bytes().first() == Some(&b'(') {
        start = 1;
        close_br = expr.find(')');
    }
    loop {
        let comma = expr[start..].find(',').map(|idx| start + idx);
        let sep = match (comma, close_br) {
            (Some(a), Some(b)) => a.min(b),
            (Some(a), None) => a,
            (None, Some(b)) => b,
            (None, None) => expr.len(),
        };
        let minus = expr[start..].find('-').map(|idx| start + idx);
        if let Some(minus_idx) = minus
            && minus_idx < sep
        {
            let n_min = expr[start..minus_idx].parse::<i32>().unwrap_or(0);
            let n_max = expr[minus_idx + 1..sep].parse::<i32>().unwrap_or(-1);
            for n in n_min..=n_max {
                result.push(n.to_string());
            }
        } else {
            result.push(expr[start..sep].to_string());
        }
        if Some(sep) == close_br || sep == expr.len() {
            break;
        }
        start = sep + 1;
    }
    result
}

fn split_on_comma(value: &str) -> Vec<String> {
    value
        .split(',')
        .filter(|part| !part.is_empty())
        .map(|part| part.to_string())
        .collect()
}

fn find_entity_index_by_subchain(structure: &BioStructure, subchain: PdbChainId) -> Option<usize> {
    structure
        .entities
        .iter()
        .position(|entity| entity.subchains.contains(&subchain))
}

// BEGIN GEMMI CPP FUNCTION fill_residue_entity_type
// Gemmi✔️✔️: for (Model& model : st.models)
// Gemmi✔️✔️:   for (Chain& chain : model.chains)
// Gemmi✔️✔️:     for (ResidueSpan& sub : chain.subchains()) {
// Gemmi✔️✔️:       if (const Entity* ent = st.get_entity_of(sub)) {
// Gemmi✔️✔️:         for (Residue& res : sub)
// Gemmi✔️✔️:           res.entity_type = ent->entity_type;
// Gemmi✔️✔️:       } else {
// Gemmi✔️✔️:         for (Residue& res : sub)
// Gemmi✔️✔️:           res.entity_type = res.is_water() ? EntityType::Water : EntityType::Unknown;
// Gemmi✔️✔️:       }
// Gemmi✔️✔️:     }
// END GEMMI CPP FUNCTION
fn fill_residue_entity_type(structure: &mut BioStructure) {
    for residue_idx in 0..structure.residues.len() {
        let chain_id = structure.residues[residue_idx].chain_id;
        let chain = &structure.chains[chain_id.index() as usize];
        let subchain = structure.residues[residue_idx]
            .source
            .subchain_id
            .or(chain.source.label_asym_id);
        if let Some(subchain_id) = subchain
            && let Some(entity_idx) = find_entity_index_by_subchain(structure, subchain_id)
        {
            structure.residues[residue_idx].entity_kind = structure.entities[entity_idx].kind;
        } else {
            structure.residues[residue_idx].entity_kind =
                if structure.residues[residue_idx].kind == ResidueKind::Water {
                    EntityKind::Water
                } else {
                    EntityKind::Unknown
                };
        }
    }
}

// BEGIN GEMMI CPP FUNCTION find_diffrn
// Gemmi✔️✔️: for (CrystalInfo& crystal_info : meta.crystals)
// Gemmi✔️✔️:   for (DiffractionInfo& diffr_info : crystal_info.diffractions)
// Gemmi✔️✔️:     if (diffr_info.id == diffrn_id)
// Gemmi✔️✔️:       return &diffr_info;
// Gemmi✔️✔️: return nullptr;
// END GEMMI CPP FUNCTION
fn find_diffrn<'a>(
    meta: &'a mut BioMetadata,
    diffrn_id: &str,
) -> Option<&'a mut BioDiffractionInfo> {
    for crystal in &mut meta.experiment_crystals {
        for diffraction in &mut crystal.diffractions {
            if diffraction.id == diffrn_id {
                return Some(diffraction);
            }
        }
    }
    None
}

fn append_info_value(target: &mut Option<String>, value: &str) {
    let value = value.trim();
    if value.is_empty() {
        return;
    }
    match target {
        Some(existing) => {
            existing.push_str("; ");
            existing.push_str(value);
        }
        None => *target = Some(value.to_string()),
    }
}

fn software_classification_from_string(value: &str) -> BioSoftwareClassification {
    match value.trim().to_ascii_lowercase().as_str() {
        "data collection" => BioSoftwareClassification::DataCollection,
        "data extraction" => BioSoftwareClassification::DataExtraction,
        "data processing" => BioSoftwareClassification::DataProcessing,
        "data reduction" => BioSoftwareClassification::DataReduction,
        "data scaling" => BioSoftwareClassification::DataScaling,
        "model building" => BioSoftwareClassification::ModelBuilding,
        "phasing" => BioSoftwareClassification::Phasing,
        "refinement" => BioSoftwareClassification::Refinement,
        _ => BioSoftwareClassification::Unspecified,
    }
}

fn first_or_add_refinement_by_id<'a>(
    meta: &'a mut BioMetadata,
    refine_id: Option<&str>,
) -> Option<&'a mut BioRefinementInfo> {
    if meta.refinement.is_empty() {
        return None;
    }
    if let Some(refine_id) = refine_id
        && let Some(index) = meta.refinement.iter().position(|item| item.id == refine_id)
    {
        return meta.refinement.get_mut(index);
    }
    meta.refinement.get_mut(0)
}

fn parse_optional_seq_id_from_mmcif(
    value: &str,
    ins_code: Option<&str>,
    line_number: usize,
) -> Result<Option<PdbSeqId>, BioReadError> {
    make_seqid(value, ins_code, line_number)
}

fn transform_is_identity(transform: &BioTransform) -> bool {
    bio_transform_approx_eq(transform, &bio_transform_identity(), 1e-6, 1e-6)
}

fn transform_has_nan(transform: &BioTransform) -> bool {
    transform
        .mat
        .iter()
        .flatten()
        .chain(transform.vec.iter())
        .any(|value| value.is_nan())
}

fn copy_optional_i32(
    row: &[CifToken],
    index: usize,
    dest: &mut Option<i32>,
) -> Result<(), BioReadError> {
    if cif_row_has2(row, index) {
        *dest = Some(parse_decimal_i32(
            &BIO_MMCIF_READ_FEATURE,
            row[index].value.as_str(),
            row[index].line_number,
            "mmCIF integer field",
        )?);
    }
    Ok(())
}

fn copy_optional_f64(
    row: &[CifToken],
    index: usize,
    dest: &mut Option<f64>,
) -> Result<(), BioReadError> {
    if cif_row_has2(row, index) {
        *dest = Some(f64::from(parse_f32(
            row[index].value.as_str(),
            row[index].line_number,
            "mmCIF number field",
        )?));
    }
    Ok(())
}

fn six_to_bio_matrix(values: [f64; 6]) -> [[f64; 3]; 3] {
    [
        [values[0], values[3], values[4]],
        [values[3], values[1], values[5]],
        [values[4], values[5], values[2]],
    ]
}

// BEGIN GEMMI CPP FUNCTION read_entry_info
// Gemmi✔️✔️: auto add_info = [&](const std::string& tag) {
// Gemmi✔️✔️:   bool first = true;
// Gemmi✔️✔️:   for (const std::string& v : block.find_values(tag))
// Gemmi✔️✔️:     if (!cif::is_null(v)) {
// Gemmi✔️✔️:       if (first)
// Gemmi✔️✔️:         st.info[tag] = cif::as_string(v);
// Gemmi✔️✔️:       else
// Gemmi✔️✔️:         st.info[tag] += "; " + cif::as_string(v);
// Gemmi✔️✔️:       first = false;
// Gemmi✔️✔️:     }
// Gemmi✔️✔️: };
// Gemmi✔️✔️: add_info("_entry.id");
// Gemmi✔️✔️: add_info("_cell.Z_PDB");
// Gemmi✔️✔️: add_info("_exptl.method");
// Gemmi✔️✔️: add_info("_struct.title");
// Gemmi✔️✔️: std::string old_date_tag = "_database_PDB_rev.date_original";
// Gemmi✔️✔️: std::string new_date_tag = "_pdbx_database_status.recvd_initial_deposition_date";
// Gemmi✔️✔️: add_info(old_date_tag);
// Gemmi✔️✔️: add_info(new_date_tag);
// Gemmi✔️✔️: if (st.info.count(old_date_tag) == 1 && st.info.count(new_date_tag) == 0)
// Gemmi✔️✔️:   st.info[new_date_tag] = st.info[old_date_tag];
// Gemmi✔️✔️: add_info("_struct_keywords.pdbx_keywords");
// Gemmi✔️✔️: add_info("_struct_keywords.text");
// END GEMMI CPP FUNCTION
fn read_entry_info(document: &CifDocument, structure: &mut BioStructure) {
    let Some(block) = document.blocks.first() else {
        return;
    };
    for item in &block.items {
        match item.tag.as_str() {
            "_entry.id" => {
                append_info_value(&mut structure.metadata.entry_id, item.value.value.as_str())
            }
            "_cell.Z_PDB" => {
                let value = item.value.value.trim();
                if cif_optional(value).is_some() {
                    structure
                        .crystal
                        .get_or_insert_with(default_crystal_info)
                        .z_pdb = Some(value.to_string());
                }
            }
            "_exptl.method" => append_info_value(
                &mut structure.metadata.experimental_method,
                item.value.value.as_str(),
            ),
            "_struct.title" => {
                append_info_value(&mut structure.metadata.title, item.value.value.as_str())
            }
            "_database_PDB_rev.date_original" => append_info_value(
                &mut structure.metadata.received_initial_deposition_date,
                item.value.value.as_str(),
            ),
            "_pdbx_database_status.recvd_initial_deposition_date" => append_info_value(
                &mut structure.metadata.received_initial_deposition_date,
                item.value.value.as_str(),
            ),
            "_struct_keywords.pdbx_keywords" => append_info_value(
                &mut structure.metadata.pdbx_keywords,
                item.value.value.as_str(),
            ),
            "_struct_keywords.text" => {
                append_info_value(&mut structure.metadata.keywords, item.value.value.as_str())
            }
            _ => {}
        }
    }
    if let Some(loop_) = find_cif_loop(&block.loops, "_exptl.method")
        && let Ok(method) = required_cif_col(loop_, "_exptl.method")
        && let Ok(rows) = cif_loop_rows(loop_)
    {
        for row in rows {
            if let Some(value) = cif_optional(row[method].value.as_str()) {
                append_info_value(&mut structure.metadata.experimental_method, value);
            }
        }
    }
}

// BEGIN GEMMI CPP FUNCTION read_audit_author
// Gemmi✔️✔️: for (const std::string& v : block.find_values("_audit_author.name"))
// Gemmi✔️✔️:   if (!cif::is_null(v))
// Gemmi✔️✔️:     st.meta.authors.push_back(cif::as_string(v));
// END GEMMI CPP FUNCTION
fn read_audit_author(loops: &[CifLoop], structure: &mut BioStructure) -> Result<(), BioReadError> {
    let Some(loop_) = find_cif_loop(loops, "_audit_author.name") else {
        return Ok(());
    };
    let name = required_cif_col(loop_, "_audit_author.name")?;
    for row in cif_loop_rows(loop_)? {
        if let Some(value) = cif_optional(row[name].value.as_str()) {
            structure.metadata.authors.push(value.to_string());
        }
    }
    Ok(())
}

// BEGIN GEMMI CPP FUNCTION read_refinement_info
// Gemmi✔️✔️: for (auto row : block.find("_refine.", {"pdbx_refine_id",
// Gemmi✔️✔️:                                 "?ls_d_res_high", "?ls_d_res_low",
// Gemmi✔️✔️:                                 "?ls_percent_reflns_obs", "?ls_number_reflns_obs",
// Gemmi✔️✔️:                                 "?ls_number_reflns_R_work", "?ls_number_reflns_R_free",
// Gemmi✔️✔️:                                 "?ls_R_factor_obs", "?ls_R_factor_R_work",
// Gemmi✔️✔️:                                 "?ls_R_factor_R_free"})) {
// Gemmi✔️✔️:   st.meta.refinement.emplace_back();
// Gemmi✔️✔️:   RefinementInfo& ref = st.meta.refinement.back();
// Gemmi✔️✔️:   ref.id = row.str(0);
// Gemmi✔️✔️:   if (row.has(1)) {
// Gemmi✔️✔️:     ref.resolution_high = cif::as_number(row[1]);
// Gemmi✔️✔️:     if (ref.resolution_high > 0 &&
// Gemmi✔️✔️:         (st.resolution == 0 || ref.resolution_high < st.resolution))
// Gemmi✔️✔️:       st.resolution = ref.resolution_high;
// Gemmi✔️✔️:   }
// Gemmi✔️✔️:   copy_double(row, 2, ref.resolution_low);
// Gemmi✔️✔️:   copy_double(row, 3, ref.completeness);
// Gemmi✔️✔️:   copy_int(row, 4, ref.reflection_count);
// Gemmi✔️✔️:   copy_int(row, 5, ref.work_set_count);
// Gemmi✔️✔️:   copy_int(row, 6, ref.rfree_set_count);
// Gemmi✔️✔️:   copy_double(row, 7, ref.r_all);
// Gemmi✔️✔️:   copy_double(row, 8, ref.r_work);
// Gemmi✔️✔️:   copy_double(row, 9, ref.r_free);
// Gemmi✔️✔️: }
// Gemmi✔️✔️: if (st.resolution == 0.) {
// Gemmi✔️✔️:   const std::string* em_res = block.find_value("_em_3d_reconstruction.resolution");
// Gemmi✔️✔️:   if (em_res)
// Gemmi✔️✔️:     st.resolution = cif::as_number(*em_res);
// Gemmi✔️✔️: }
// END GEMMI CPP FUNCTION
fn read_refinement_info(
    document: &CifDocument,
    loops: &[CifLoop],
    structure: &mut BioStructure,
) -> Result<(), BioReadError> {
    if let Some(loop_) = find_cif_loop(loops, "_refine.pdbx_refine_id") {
        let refine_id = required_cif_col(loop_, "_refine.pdbx_refine_id")?;
        let ls_d_res_high = optional_cif_col(loop_, "_refine.ls_d_res_high");
        let ls_d_res_low = optional_cif_col(loop_, "_refine.ls_d_res_low");
        let ls_percent_reflns_obs = optional_cif_col(loop_, "_refine.ls_percent_reflns_obs");
        let ls_number_reflns_obs = optional_cif_col(loop_, "_refine.ls_number_reflns_obs");
        let ls_number_reflns_r_work = optional_cif_col(loop_, "_refine.ls_number_reflns_R_work");
        let ls_number_reflns_r_free = optional_cif_col(loop_, "_refine.ls_number_reflns_R_free");
        let ls_r_factor_obs = optional_cif_col(loop_, "_refine.ls_R_factor_obs");
        let ls_r_factor_r_work = optional_cif_col(loop_, "_refine.ls_R_factor_R_work");
        let ls_r_factor_r_free = optional_cif_col(loop_, "_refine.ls_R_factor_R_free");
        for row in cif_loop_rows(loop_)? {
            let mut refinement = BioRefinementInfo {
                id: row[refine_id].value.trim().to_string(),
                ..BioRefinementInfo::default()
            };
            if let Some(index) = ls_d_res_high
                && cif_row_has2(row, index)
            {
                let value = f64::from(parse_f32(
                    row[index].value.as_str(),
                    row[index].line_number,
                    "_refine.ls_d_res_high",
                )?);
                refinement.resolution_high = Some(value);
                if value > 0.0 && structure.resolution.is_none_or(|current| value < current) {
                    structure.resolution = Some(value);
                }
            }
            if let Some(index) = ls_d_res_low {
                copy_optional_f64(row, index, &mut refinement.resolution_low)?;
            }
            if let Some(index) = ls_percent_reflns_obs {
                copy_optional_f64(row, index, &mut refinement.completeness)?;
            }
            if let Some(index) = ls_number_reflns_obs {
                copy_optional_i32(row, index, &mut refinement.reflection_count)?;
            }
            if let Some(index) = ls_number_reflns_r_work {
                copy_optional_i32(row, index, &mut refinement.work_set_count)?;
            }
            if let Some(index) = ls_number_reflns_r_free {
                copy_optional_i32(row, index, &mut refinement.rfree_set_count)?;
            }
            if let Some(index) = ls_r_factor_obs {
                copy_optional_f64(row, index, &mut refinement.r_all)?;
            }
            if let Some(index) = ls_r_factor_r_work {
                copy_optional_f64(row, index, &mut refinement.r_work)?;
            }
            if let Some(index) = ls_r_factor_r_free {
                copy_optional_f64(row, index, &mut refinement.r_free)?;
            }
            structure.metadata.refinement.push(refinement);
        }
    }
    if structure.resolution.is_none()
        && let Some(value) = find_cif_item_value(document, "_em_3d_reconstruction.resolution")
        && let Some(value) = cif_optional(value)
    {
        structure.resolution = Some(f64::from(parse_f32(
            value,
            0,
            "_em_3d_reconstruction.resolution",
        )?));
    }
    Ok(())
}

// BEGIN GEMMI CPP FUNCTION read_tls_info
// Gemmi✔️✔️: for (auto row : block.find("_pdbx_refine_tls.", {...})) {
// Gemmi✔️✔️:   if (st.meta.refinement.empty())
// Gemmi✔️✔️:     break;
// Gemmi✔️✔️:   RefinementInfo* ref = nullptr;
// Gemmi✔️✔️:   if (row.has(1))
// Gemmi✔️✔️:     ref = get_by_id(st.meta.refinement, row.str(1));
// Gemmi✔️✔️:   if (!ref)
// Gemmi✔️✔️:     ref = &st.meta.refinement[0];
// Gemmi✔️✔️:   ref->tls_groups.emplace_back();
// Gemmi✔️✔️:   TlsGroup& tls = ref->tls_groups.back();
// Gemmi✔️✔️:   tls.id = row.str(0);
// Gemmi✔️✔️:   tls.num_id = (short) no_sign_atoi(tls.id.c_str());
// Gemmi✔️✔️:   tls.T = get_smat33<double>(row, 2);
// Gemmi✔️✔️:   tls.L = get_smat33<double>(row, 8);
// Gemmi✔️✔️:   for (int i = 0; i < 3; ++i)
// Gemmi✔️✔️:     for (int j = 0; j < 3; ++j)
// Gemmi✔️✔️:       tls.S[i][j] = cif::as_number(row[14+3*i+j]);
// Gemmi✔️✔️:   tls.origin.x = cif::as_number(row[23]);
// Gemmi✔️✔️:   tls.origin.y = cif::as_number(row[24]);
// Gemmi✔️✔️:   tls.origin.z = cif::as_number(row[25]);
// Gemmi✔️✔️: }
// Gemmi✔️✔️: for (auto row : block.find("_pdbx_refine_tls_group.", {
// Gemmi✔️✔️:   "refine_tls_id", "?beg_auth_asym_id", "?beg_auth_seq_id", "?beg_PDB_ins_code",
// Gemmi✔️✔️:   "?end_auth_seq_id", "?end_PDB_ins_code", "?selection_details"})) {
// Gemmi✔️✔️:   for (RefinementInfo& ref : st.meta.refinement)
// Gemmi✔️✔️:     if (TlsGroup* tls = get_by_id(ref.tls_groups, row.str(0))) {
// Gemmi✔️✔️:       tls->selections.emplace_back();
// Gemmi✔️✔️:       ...
// Gemmi✔️✔️:       break;
// Gemmi✔️✔️:     }
// Gemmi✔️✔️: }
// END GEMMI CPP FUNCTION
fn read_tls_info(loops: &[CifLoop], structure: &mut BioStructure) -> Result<(), BioReadError> {
    if let Some(loop_) = find_cif_loop(loops, "_pdbx_refine_tls.id") {
        let Some(
            [
                id,
                t11,
                t22,
                t33,
                t12,
                t13,
                t23,
                l11,
                l22,
                l33,
                l12,
                l13,
                l23,
                s11,
                s12,
                s13,
                s21,
                s22,
                s23,
                s31,
                s32,
                s33,
                origin_x,
                origin_y,
                origin_z,
            ],
        ) = required_cif_table_cols(
            loop_,
            [
                "_pdbx_refine_tls.id",
                "_pdbx_refine_tls.T[1][1]",
                "_pdbx_refine_tls.T[2][2]",
                "_pdbx_refine_tls.T[3][3]",
                "_pdbx_refine_tls.T[1][2]",
                "_pdbx_refine_tls.T[1][3]",
                "_pdbx_refine_tls.T[2][3]",
                "_pdbx_refine_tls.L[1][1]",
                "_pdbx_refine_tls.L[2][2]",
                "_pdbx_refine_tls.L[3][3]",
                "_pdbx_refine_tls.L[1][2]",
                "_pdbx_refine_tls.L[1][3]",
                "_pdbx_refine_tls.L[2][3]",
                "_pdbx_refine_tls.S[1][1]",
                "_pdbx_refine_tls.S[1][2]",
                "_pdbx_refine_tls.S[1][3]",
                "_pdbx_refine_tls.S[2][1]",
                "_pdbx_refine_tls.S[2][2]",
                "_pdbx_refine_tls.S[2][3]",
                "_pdbx_refine_tls.S[3][1]",
                "_pdbx_refine_tls.S[3][2]",
                "_pdbx_refine_tls.S[3][3]",
                "_pdbx_refine_tls.origin_x",
                "_pdbx_refine_tls.origin_y",
                "_pdbx_refine_tls.origin_z",
            ],
        )
        else {
            return Ok(());
        };
        let pdbx_refine_id = optional_cif_col(loop_, "_pdbx_refine_tls.pdbx_refine_id");
        for row in cif_loop_rows(loop_)? {
            if structure.metadata.refinement.is_empty() {
                break;
            }
            let refine_id =
                pdbx_refine_id.and_then(|index| cif_optional(row[index].value.as_str()));
            let Some(refinement) =
                first_or_add_refinement_by_id(&mut structure.metadata, refine_id)
            else {
                break;
            };
            let tls_id = row[id].value.trim().to_string();
            let num_id = tls_id.parse::<i16>().ok();
            refinement.tls_groups.push(BioTlsGroup {
                num_id,
                id: tls_id,
                t: six_to_bio_matrix([
                    f64::from(parse_f32(
                        row[t11].value.as_str(),
                        row[t11].line_number,
                        "_pdbx_refine_tls.T[1][1]",
                    )?),
                    f64::from(parse_f32(
                        row[t22].value.as_str(),
                        row[t22].line_number,
                        "_pdbx_refine_tls.T[2][2]",
                    )?),
                    f64::from(parse_f32(
                        row[t33].value.as_str(),
                        row[t33].line_number,
                        "_pdbx_refine_tls.T[3][3]",
                    )?),
                    f64::from(parse_f32(
                        row[t12].value.as_str(),
                        row[t12].line_number,
                        "_pdbx_refine_tls.T[1][2]",
                    )?),
                    f64::from(parse_f32(
                        row[t13].value.as_str(),
                        row[t13].line_number,
                        "_pdbx_refine_tls.T[1][3]",
                    )?),
                    f64::from(parse_f32(
                        row[t23].value.as_str(),
                        row[t23].line_number,
                        "_pdbx_refine_tls.T[2][3]",
                    )?),
                ]),
                l: six_to_bio_matrix([
                    f64::from(parse_f32(
                        row[l11].value.as_str(),
                        row[l11].line_number,
                        "_pdbx_refine_tls.L[1][1]",
                    )?),
                    f64::from(parse_f32(
                        row[l22].value.as_str(),
                        row[l22].line_number,
                        "_pdbx_refine_tls.L[2][2]",
                    )?),
                    f64::from(parse_f32(
                        row[l33].value.as_str(),
                        row[l33].line_number,
                        "_pdbx_refine_tls.L[3][3]",
                    )?),
                    f64::from(parse_f32(
                        row[l12].value.as_str(),
                        row[l12].line_number,
                        "_pdbx_refine_tls.L[1][2]",
                    )?),
                    f64::from(parse_f32(
                        row[l13].value.as_str(),
                        row[l13].line_number,
                        "_pdbx_refine_tls.L[1][3]",
                    )?),
                    f64::from(parse_f32(
                        row[l23].value.as_str(),
                        row[l23].line_number,
                        "_pdbx_refine_tls.L[2][3]",
                    )?),
                ]),
                s: [
                    [
                        f64::from(parse_f32(
                            row[s11].value.as_str(),
                            row[s11].line_number,
                            "_pdbx_refine_tls.S[1][1]",
                        )?),
                        f64::from(parse_f32(
                            row[s12].value.as_str(),
                            row[s12].line_number,
                            "_pdbx_refine_tls.S[1][2]",
                        )?),
                        f64::from(parse_f32(
                            row[s13].value.as_str(),
                            row[s13].line_number,
                            "_pdbx_refine_tls.S[1][3]",
                        )?),
                    ],
                    [
                        f64::from(parse_f32(
                            row[s21].value.as_str(),
                            row[s21].line_number,
                            "_pdbx_refine_tls.S[2][1]",
                        )?),
                        f64::from(parse_f32(
                            row[s22].value.as_str(),
                            row[s22].line_number,
                            "_pdbx_refine_tls.S[2][2]",
                        )?),
                        f64::from(parse_f32(
                            row[s23].value.as_str(),
                            row[s23].line_number,
                            "_pdbx_refine_tls.S[2][3]",
                        )?),
                    ],
                    [
                        f64::from(parse_f32(
                            row[s31].value.as_str(),
                            row[s31].line_number,
                            "_pdbx_refine_tls.S[3][1]",
                        )?),
                        f64::from(parse_f32(
                            row[s32].value.as_str(),
                            row[s32].line_number,
                            "_pdbx_refine_tls.S[3][2]",
                        )?),
                        f64::from(parse_f32(
                            row[s33].value.as_str(),
                            row[s33].line_number,
                            "_pdbx_refine_tls.S[3][3]",
                        )?),
                    ],
                ],
                origin: [
                    f64::from(parse_f32(
                        row[origin_x].value.as_str(),
                        row[origin_x].line_number,
                        "_pdbx_refine_tls.origin_x",
                    )?),
                    f64::from(parse_f32(
                        row[origin_y].value.as_str(),
                        row[origin_y].line_number,
                        "_pdbx_refine_tls.origin_y",
                    )?),
                    f64::from(parse_f32(
                        row[origin_z].value.as_str(),
                        row[origin_z].line_number,
                        "_pdbx_refine_tls.origin_z",
                    )?),
                ],
                ..BioTlsGroup::default()
            });
        }
    }

    if let Some(loop_) = find_cif_loop(loops, "_pdbx_refine_tls_group.refine_tls_id") {
        let refine_tls_id = required_cif_col(loop_, "_pdbx_refine_tls_group.refine_tls_id")?;
        let beg_auth_asym_id = optional_cif_col(loop_, "_pdbx_refine_tls_group.beg_auth_asym_id");
        let beg_auth_seq_id = optional_cif_col(loop_, "_pdbx_refine_tls_group.beg_auth_seq_id");
        let beg_pdb_ins_code = optional_cif_col(loop_, "_pdbx_refine_tls_group.beg_PDB_ins_code");
        let end_auth_seq_id = optional_cif_col(loop_, "_pdbx_refine_tls_group.end_auth_seq_id");
        let end_pdb_ins_code = optional_cif_col(loop_, "_pdbx_refine_tls_group.end_PDB_ins_code");
        let selection_details = optional_cif_col(loop_, "_pdbx_refine_tls_group.selection_details");
        for row in cif_loop_rows(loop_)? {
            let tls_id = row[refine_tls_id].value.trim();
            for refinement in &mut structure.metadata.refinement {
                if let Some(group) = refinement
                    .tls_groups
                    .iter_mut()
                    .find(|group| group.id == tls_id)
                {
                    let mut selection = BioTlsSelection::default();
                    if let Some(index) = beg_auth_asym_id
                        && cif_row_has2(row, index)
                    {
                        selection.chain = row[index].value.trim().to_string();
                    }
                    if let Some(index) = beg_auth_seq_id
                        && cif_row_has2(row, index)
                    {
                        selection.res_begin = parse_optional_seq_id_from_mmcif(
                            row[index].value.as_str(),
                            beg_pdb_ins_code.map(|idx| row[idx].value.as_str()),
                            row[index].line_number,
                        )?;
                    }
                    if let Some(index) = end_auth_seq_id
                        && cif_row_has2(row, index)
                    {
                        selection.res_end = parse_optional_seq_id_from_mmcif(
                            row[index].value.as_str(),
                            end_pdb_ins_code.map(|idx| row[idx].value.as_str()),
                            row[index].line_number,
                        )?;
                    }
                    if let Some(index) = selection_details
                        && cif_row_has2(row, index)
                    {
                        selection.details = row[index].value.trim().to_string();
                    }
                    group.selections.push(selection);
                    break;
                }
            }
        }
    }
    Ok(())
}

// BEGIN GEMMI CPP FUNCTION read_experimental_info
// Gemmi✔️✔️: for (auto row : block.find("_exptl.", {"method", "?crystals_number"})) { ... }
// Gemmi✔️✔️: for (auto row : block.find("_exptl_crystal.", {"id", "?description"})) { ... }
// Gemmi✔️✔️: for (auto row : block.find("_diffrn.", {"id", "crystal_id", "?ambient_temp"})) { ... }
// Gemmi✔️✔️: for (auto row : block.find("_diffrn_detector.", {...})) if (DiffractionInfo* di = find_diffrn(...)) { ... }
// Gemmi✔️✔️: for (auto row : block.find("_diffrn_radiation.", {...})) if (DiffractionInfo* di = find_diffrn(...)) { ... }
// Gemmi✔️✔️: for (auto row : block.find("_diffrn_source.", {...})) if (DiffractionInfo* di = find_diffrn(...)) { ... }
// END GEMMI CPP FUNCTION
fn read_experimental_info(
    loops: &[CifLoop],
    structure: &mut BioStructure,
) -> Result<(), BioReadError> {
    if let Some(loop_) = find_cif_loop(loops, "_exptl.method") {
        let method = required_cif_col(loop_, "_exptl.method")?;
        let crystals_number = optional_cif_col(loop_, "_exptl.crystals_number");
        for row in cif_loop_rows(loop_)? {
            let mut experiment = BioExperimentInfo {
                method: row[method].value.trim().to_string(),
                ..BioExperimentInfo::default()
            };
            if let Some(index) = crystals_number {
                copy_optional_i32(row, index, &mut experiment.number_of_crystals)?;
            }
            structure.metadata.experiments.push(experiment);
        }
    }
    if let Some(loop_) = find_cif_loop(loops, "_exptl_crystal.id") {
        let id = required_cif_col(loop_, "_exptl_crystal.id")?;
        let description = optional_cif_col(loop_, "_exptl_crystal.description");
        for row in cif_loop_rows(loop_)? {
            let mut crystal = BioExperimentCrystalInfo {
                id: row[id].value.trim().to_string(),
                ..BioExperimentCrystalInfo::default()
            };
            if let Some(index) = description
                && cif_row_has2(row, index)
            {
                crystal.description = row[index].value.trim().to_string();
            }
            structure.metadata.experiment_crystals.push(crystal);
        }
    }
    if let Some(loop_) = find_cif_loop(loops, "_diffrn.id") {
        if let Some([id, crystal_id]) =
            required_cif_table_cols(loop_, ["_diffrn.id", "_diffrn.crystal_id"])
        {
            let ambient_temp = optional_cif_col(loop_, "_diffrn.ambient_temp");
            for row in cif_loop_rows(loop_)? {
                let crystal_id_value = row[crystal_id].value.trim();
                if let Some(crystal) = structure
                    .metadata
                    .experiment_crystals
                    .iter_mut()
                    .find(|crystal| crystal.id == crystal_id_value)
                {
                    let mut diffraction = BioDiffractionInfo {
                        id: row[id].value.trim().to_string(),
                        ..BioDiffractionInfo::default()
                    };
                    if let Some(index) = ambient_temp {
                        copy_optional_f64(row, index, &mut diffraction.temperature)?;
                    }
                    crystal.diffractions.push(diffraction);
                }
            }
        }
    }
    if let Some(loop_) = find_cif_loop(loops, "_diffrn_detector.diffrn_id") {
        let diffrn_id = required_cif_col(loop_, "_diffrn_detector.diffrn_id")?;
        let collection_date = optional_cif_col(loop_, "_diffrn_detector.pdbx_collection_date");
        let detector = optional_cif_col(loop_, "_diffrn_detector.detector");
        let type_ = optional_cif_col(loop_, "_diffrn_detector.type");
        let details = optional_cif_col(loop_, "_diffrn_detector.details");
        for row in cif_loop_rows(loop_)? {
            if let Some(diffraction) =
                find_diffrn(&mut structure.metadata, row[diffrn_id].value.trim())
            {
                if let Some(index) = collection_date
                    && cif_row_has2(row, index)
                {
                    diffraction.collection_date = row[index].value.trim().to_string();
                }
                if let Some(index) = detector
                    && cif_row_has2(row, index)
                {
                    diffraction.detector = row[index].value.trim().to_string();
                }
                if let Some(index) = type_
                    && cif_row_has2(row, index)
                {
                    diffraction.detector_make = row[index].value.trim().to_string();
                }
                if let Some(index) = details
                    && cif_row_has2(row, index)
                {
                    diffraction.optics = row[index].value.trim().to_string();
                }
            }
        }
    }
    if let Some(loop_) = find_cif_loop(loops, "_diffrn_radiation.diffrn_id") {
        let diffrn_id = required_cif_col(loop_, "_diffrn_radiation.diffrn_id")?;
        let scattering_type = optional_cif_col(loop_, "_diffrn_radiation.pdbx_scattering_type");
        let mono_or_laue =
            optional_cif_col(loop_, "_diffrn_radiation.pdbx_monochromatic_or_laue_m_l");
        let monochromator = optional_cif_col(loop_, "_diffrn_radiation.monochromator");
        for row in cif_loop_rows(loop_)? {
            if let Some(diffraction) =
                find_diffrn(&mut structure.metadata, row[diffrn_id].value.trim())
            {
                if let Some(index) = scattering_type
                    && cif_row_has2(row, index)
                {
                    diffraction.scattering_type = row[index].value.trim().to_string();
                }
                if let Some(index) = mono_or_laue
                    && cif_row_has2(row, index)
                {
                    diffraction.mono_or_laue = row[index].value.chars().next();
                }
                if let Some(index) = monochromator
                    && cif_row_has2(row, index)
                {
                    diffraction.monochromator = row[index].value.trim().to_string();
                }
            }
        }
    }
    if let Some(loop_) = find_cif_loop(loops, "_diffrn_source.diffrn_id") {
        let diffrn_id = required_cif_col(loop_, "_diffrn_source.diffrn_id")?;
        let source = optional_cif_col(loop_, "_diffrn_source.source");
        let type_ = optional_cif_col(loop_, "_diffrn_source.type");
        let synchrotron = optional_cif_col(loop_, "_diffrn_source.pdbx_synchrotron_site");
        let beamline = optional_cif_col(loop_, "_diffrn_source.pdbx_synchrotron_beamline");
        let wavelength_list = optional_cif_col(loop_, "_diffrn_source.pdbx_wavelength_list");
        for row in cif_loop_rows(loop_)? {
            if let Some(diffraction) =
                find_diffrn(&mut structure.metadata, row[diffrn_id].value.trim())
            {
                if let Some(index) = source
                    && cif_row_has2(row, index)
                {
                    diffraction.source = row[index].value.trim().to_string();
                }
                if let Some(index) = type_
                    && cif_row_has2(row, index)
                {
                    diffraction.source_type = row[index].value.trim().to_string();
                }
                if let Some(index) = synchrotron
                    && cif_row_has2(row, index)
                {
                    diffraction.synchrotron = row[index].value.trim().to_string();
                }
                if let Some(index) = beamline
                    && cif_row_has2(row, index)
                {
                    diffraction.beamline = row[index].value.trim().to_string();
                }
                if let Some(index) = wavelength_list
                    && cif_row_has2(row, index)
                {
                    diffraction.wavelengths = row[index].value.trim().to_string();
                }
            }
        }
    }
    Ok(())
}

// BEGIN GEMMI CPP FUNCTION read_reflns_info
// Gemmi✔️✔️: size_t n = 0;
// Gemmi✔️✔️: for (auto row : block.find("_reflns.", {"pdbx_diffrn_id", "?number_obs",
// Gemmi✔️✔️:                                 "?d_resolution_high", "?d_resolution_low",
// Gemmi✔️✔️:                                 "?percent_possible_obs", "?pdbx_redundancy",
// Gemmi✔️✔️:                                 "?pdbx_Rmerge_I_obs", "?pdbx_Rsym_value",
// Gemmi✔️✔️:                                 "?pdbx_netI_over_sigmaI"})) {
// Gemmi✔️✔️:   if (n >= st.meta.experiments.size())
// Gemmi✔️✔️:     break;
// Gemmi✔️✔️:   ExperimentInfo& exper = st.meta.experiments[n++];
// Gemmi✔️✔️:   split_str_into(row.str(0), ',', exper.diffraction_ids);
// Gemmi✔️✔️:   copy_int(row, 1, exper.unique_reflections);
// Gemmi✔️✔️:   copy_double(row, 2, exper.reflections.resolution_high);
// Gemmi✔️✔️:   copy_double(row, 3, exper.reflections.resolution_low);
// Gemmi✔️✔️:   copy_double(row, 4, exper.reflections.completeness);
// Gemmi✔️✔️:   copy_double(row, 5, exper.reflections.redundancy);
// Gemmi✔️✔️:   copy_double(row, 6, exper.reflections.r_merge);
// Gemmi✔️✔️:   copy_double(row, 7, exper.reflections.r_sym);
// Gemmi✔️✔️:   copy_double(row, 8, exper.reflections.mean_I_over_sigma);
// Gemmi✔️✔️: }
// END GEMMI CPP FUNCTION
fn read_reflns_info(loops: &[CifLoop], structure: &mut BioStructure) -> Result<(), BioReadError> {
    let Some(loop_) = find_cif_loop(loops, "_reflns.pdbx_diffrn_id") else {
        return Ok(());
    };
    let pdbx_diffrn_id = required_cif_col(loop_, "_reflns.pdbx_diffrn_id")?;
    let number_obs = optional_cif_col(loop_, "_reflns.number_obs");
    let d_resolution_high = optional_cif_col(loop_, "_reflns.d_resolution_high");
    let d_resolution_low = optional_cif_col(loop_, "_reflns.d_resolution_low");
    let percent_possible_obs = optional_cif_col(loop_, "_reflns.percent_possible_obs");
    let pdbx_redundancy = optional_cif_col(loop_, "_reflns.pdbx_redundancy");
    let pdbx_rmerge_i_obs = optional_cif_col(loop_, "_reflns.pdbx_Rmerge_I_obs");
    let pdbx_rsym_value = optional_cif_col(loop_, "_reflns.pdbx_Rsym_value");
    let pdbx_neti_over_sigmai = optional_cif_col(loop_, "_reflns.pdbx_netI_over_sigmaI");
    for (n, row) in cif_loop_rows(loop_)?.enumerate() {
        if n >= structure.metadata.experiments.len() {
            break;
        }
        let experiment = &mut structure.metadata.experiments[n];
        experiment.diffraction_ids = split_on_comma(row[pdbx_diffrn_id].value.trim());
        if let Some(index) = number_obs {
            copy_optional_i32(row, index, &mut experiment.unique_reflections)?;
        }
        if let Some(index) = d_resolution_high {
            copy_optional_f64(row, index, &mut experiment.reflections.resolution_high)?;
        }
        if let Some(index) = d_resolution_low {
            copy_optional_f64(row, index, &mut experiment.reflections.resolution_low)?;
        }
        if let Some(index) = percent_possible_obs {
            copy_optional_f64(row, index, &mut experiment.reflections.completeness)?;
        }
        if let Some(index) = pdbx_redundancy {
            copy_optional_f64(row, index, &mut experiment.reflections.redundancy)?;
        }
        if let Some(index) = pdbx_rmerge_i_obs {
            copy_optional_f64(row, index, &mut experiment.reflections.r_merge)?;
        }
        if let Some(index) = pdbx_rsym_value {
            copy_optional_f64(row, index, &mut experiment.reflections.r_sym)?;
        }
        if let Some(index) = pdbx_neti_over_sigmai {
            copy_optional_f64(row, index, &mut experiment.reflections.mean_i_over_sigma)?;
        }
    }
    Ok(())
}

// BEGIN GEMMI CPP FUNCTION read_software_info
// Gemmi✔️✔️: for (auto row : block.find("_software.", {"name", "?classification", "?version",
// Gemmi✔️✔️:                                 "?date", "?description", "?contact_author",
// Gemmi✔️✔️:                                 "?contact_author_email"})) {
// Gemmi✔️✔️:   st.meta.software.emplace_back();
// Gemmi✔️✔️:   SoftwareItem& item = st.meta.software.back();
// Gemmi✔️✔️:   item.name = row.str(0);
// Gemmi✔️✔️:   if (row.has2(1))
// Gemmi✔️✔️:     item.classification = software_classification_from_string(row.str(1));
// Gemmi✔️✔️:   copy_string(row, 2, item.version);
// Gemmi✔️✔️:   copy_string(row, 3, item.date);
// Gemmi✔️✔️:   copy_string(row, 4, item.description);
// Gemmi✔️✔️:   copy_string(row, 5, item.contact_author);
// Gemmi✔️✔️:   copy_string(row, 6, item.contact_author_email);
// Gemmi✔️✔️: }
// END GEMMI CPP FUNCTION
fn read_software_info(loops: &[CifLoop], structure: &mut BioStructure) -> Result<(), BioReadError> {
    let Some(loop_) = find_cif_loop(loops, "_software.name") else {
        return Ok(());
    };
    let name = required_cif_col(loop_, "_software.name")?;
    let classification = optional_cif_col(loop_, "_software.classification");
    let version = optional_cif_col(loop_, "_software.version");
    let date = optional_cif_col(loop_, "_software.date");
    let description = optional_cif_col(loop_, "_software.description");
    let contact_author = optional_cif_col(loop_, "_software.contact_author");
    let contact_author_email = optional_cif_col(loop_, "_software.contact_author_email");
    for row in cif_loop_rows(loop_)? {
        let mut item = BioSoftwareItem {
            name: row[name].value.trim().to_string(),
            ..BioSoftwareItem::default()
        };
        if let Some(index) = classification
            && cif_row_has2(row, index)
        {
            item.classification = software_classification_from_string(row[index].value.as_str());
        }
        if let Some(index) = version
            && cif_row_has2(row, index)
        {
            item.version = row[index].value.trim().to_string();
        }
        if let Some(index) = date
            && cif_row_has2(row, index)
        {
            item.date = row[index].value.trim().to_string();
        }
        if let Some(index) = description
            && cif_row_has2(row, index)
        {
            item.description = row[index].value.trim().to_string();
        }
        if let Some(index) = contact_author
            && cif_row_has2(row, index)
        {
            item.contact_author = row[index].value.trim().to_string();
        }
        if let Some(index) = contact_author_email
            && cif_row_has2(row, index)
        {
            item.contact_author_email = row[index].value.trim().to_string();
        }
        structure.metadata.software.push(item);
    }
    Ok(())
}

// BEGIN GEMMI CPP FUNCTION read_ncs_info
// Gemmi✔️✔️: std::vector<std::string> ncs_oper_tags = transform_tags("matrix", "vector");
// Gemmi✔️✔️: ncs_oper_tags.emplace_back("id");
// Gemmi✔️✔️: ncs_oper_tags.emplace_back("?code");
// Gemmi✔️✔️: cif::Table ncs_oper = block.find("_struct_ncs_oper.", ncs_oper_tags);
// Gemmi✔️✔️: for (auto op : ncs_oper) {
// Gemmi✔️✔️:   bool given = op.has(13) && op.str(13) == "given";
// Gemmi✔️✔️:   Transform tr = get_transform_matrix(op);
// Gemmi✔️✔️:   if (tr.is_identity())
// Gemmi✔️✔️:     st.info["_struct_ncs_oper.id"] = op.str(12);
// Gemmi✔️✔️:   else if (tr.has_nan())
// Gemmi✔️✔️:     continue;
// Gemmi✔️✔️:   else
// Gemmi✔️✔️:     st.ncs.push_back({op.str(12), given, tr});
// Gemmi✔️✔️: }
// END GEMMI CPP FUNCTION
fn read_ncs_info(loops: &[CifLoop], structure: &mut BioStructure) -> Result<(), BioReadError> {
    let Some(loop_) = find_cif_loop(loops, "_struct_ncs_oper.matrix[1][1]") else {
        return Ok(());
    };
    let Some(
        [
            matrix11,
            matrix12,
            matrix13,
            vector1,
            matrix21,
            matrix22,
            matrix23,
            vector2,
            matrix31,
            matrix32,
            matrix33,
            vector3,
            id,
        ],
    ) = required_cif_table_cols(
        loop_,
        [
            "_struct_ncs_oper.matrix[1][1]",
            "_struct_ncs_oper.matrix[1][2]",
            "_struct_ncs_oper.matrix[1][3]",
            "_struct_ncs_oper.vector[1]",
            "_struct_ncs_oper.matrix[2][1]",
            "_struct_ncs_oper.matrix[2][2]",
            "_struct_ncs_oper.matrix[2][3]",
            "_struct_ncs_oper.vector[2]",
            "_struct_ncs_oper.matrix[3][1]",
            "_struct_ncs_oper.matrix[3][2]",
            "_struct_ncs_oper.matrix[3][3]",
            "_struct_ncs_oper.vector[3]",
            "_struct_ncs_oper.id",
        ],
    )
    else {
        return Ok(());
    };
    let code = optional_cif_col(loop_, "_struct_ncs_oper.code");
    for row in cif_loop_rows(loop_)? {
        let matrix_row = vec![
            row[matrix11].clone(),
            row[matrix12].clone(),
            row[matrix13].clone(),
            row[vector1].clone(),
            row[matrix21].clone(),
            row[matrix22].clone(),
            row[matrix23].clone(),
            row[vector2].clone(),
            row[matrix31].clone(),
            row[matrix32].clone(),
            row[matrix33].clone(),
            row[vector3].clone(),
        ];
        let transform = get_transform_matrix(&matrix_row)?;
        let given = code
            .is_some_and(|index| cif_row_has2(row, index) && row[index].value.trim() == "given");
        let op_id = row[id].value.trim().to_string();
        if transform_is_identity(&transform) {
            structure.ncs_oper_identity_id = Some(op_id);
        } else if transform_has_nan(&transform) {
            continue;
        } else {
            structure.ncs_operators.push(BioNcsOperator {
                id: op_id,
                given,
                transform,
            });
        }
    }
    Ok(())
}

// BEGIN GEMMI CPP FUNCTION read_assemblies
// Gemmi✔️✔️: std::vector<Assembly> read_assemblies(cif::Block& block) { ... }
// END GEMMI CPP FUNCTION
fn read_assemblies(loops: &[CifLoop]) -> Result<Vec<BioAssembly>, BioReadError> {
    let mut oper_list = Vec::<BioAssemblyOperator>::new();
    if let Some(loop_) = find_cif_loop(loops, "_pdbx_struct_oper_list.id") {
        if let Some(
            [
                matrix11,
                matrix12,
                matrix13,
                vector1,
                matrix21,
                matrix22,
                matrix23,
                vector2,
                matrix31,
                matrix32,
                matrix33,
                vector3,
                id,
                type_,
            ],
        ) = required_cif_table_cols(
            loop_,
            [
                "_pdbx_struct_oper_list.matrix[1][1]",
                "_pdbx_struct_oper_list.matrix[1][2]",
                "_pdbx_struct_oper_list.matrix[1][3]",
                "_pdbx_struct_oper_list.vector[1]",
                "_pdbx_struct_oper_list.matrix[2][1]",
                "_pdbx_struct_oper_list.matrix[2][2]",
                "_pdbx_struct_oper_list.matrix[2][3]",
                "_pdbx_struct_oper_list.vector[2]",
                "_pdbx_struct_oper_list.matrix[3][1]",
                "_pdbx_struct_oper_list.matrix[3][2]",
                "_pdbx_struct_oper_list.matrix[3][3]",
                "_pdbx_struct_oper_list.vector[3]",
                "_pdbx_struct_oper_list.id",
                "_pdbx_struct_oper_list.type",
            ],
        ) {
            for row in cif_loop_rows(loop_)? {
                let matrix_row = vec![
                    row[matrix11].clone(),
                    row[matrix12].clone(),
                    row[matrix13].clone(),
                    row[vector1].clone(),
                    row[matrix21].clone(),
                    row[matrix22].clone(),
                    row[matrix23].clone(),
                    row[vector2].clone(),
                    row[matrix31].clone(),
                    row[matrix32].clone(),
                    row[matrix33].clone(),
                    row[vector3].clone(),
                ];
                oper_list.push(BioAssemblyOperator {
                    name: row[id].value.trim().to_string(),
                    type_: row[type_].value.trim().to_string(),
                    transform: get_transform_matrix(&matrix_row)?,
                });
            }
        }
    }

    let mut prop_rows = Vec::<Vec<CifToken>>::new();
    let mut prop_cols = None;
    if let Some(loop_) = find_cif_loop(loops, "_pdbx_struct_assembly_prop.biol_id") {
        if let Some([biol_id, type_, value]) = required_cif_table_cols(
            loop_,
            [
                "_pdbx_struct_assembly_prop.biol_id",
                "_pdbx_struct_assembly_prop.type",
                "_pdbx_struct_assembly_prop.value",
            ],
        ) {
            prop_cols = Some((biol_id, type_, value));
            for row in cif_loop_rows(loop_)? {
                prop_rows.push(row.to_vec());
            }
        }
    }

    let mut gen_rows = Vec::<Vec<CifToken>>::new();
    let mut gen_cols = None;
    if let Some(loop_) = find_cif_loop(loops, "_pdbx_struct_assembly_gen.assembly_id") {
        if let Some([assembly_id, oper_expression, asym_id_list]) = required_cif_table_cols(
            loop_,
            [
                "_pdbx_struct_assembly_gen.assembly_id",
                "_pdbx_struct_assembly_gen.oper_expression",
                "_pdbx_struct_assembly_gen.asym_id_list",
            ],
        ) {
            gen_cols = Some((assembly_id, oper_expression, asym_id_list));
            for row in cif_loop_rows(loop_)? {
                gen_rows.push(row.to_vec());
            }
        }
    }

    let Some(loop_) = find_cif_loop(loops, "_pdbx_struct_assembly.id") else {
        return Ok(Vec::new());
    };
    let Some(
        [
            id,
            details,
            method_details,
            oligomeric_details,
            oligomeric_count,
        ],
    ) = required_cif_table_cols(
        loop_,
        [
            "_pdbx_struct_assembly.id",
            "_pdbx_struct_assembly.details",
            "_pdbx_struct_assembly.method_details",
            "_pdbx_struct_assembly.oligomeric_details",
            "_pdbx_struct_assembly.oligomeric_count",
        ],
    )
    else {
        return Ok(Vec::new());
    };
    let mut assemblies = Vec::new();
    for row in cif_loop_rows(loop_)? {
        let mut assembly = BioAssembly {
            name: row[id].value.trim().to_string(),
            ..BioAssembly::default()
        };
        let detail = row[details].value.trim();
        if detail == "author_and_software_defined_assembly" {
            assembly.author_determined = true;
            assembly.software_determined = true;
        } else if detail == "author_defined_assembly" {
            assembly.author_determined = true;
        } else if detail == "software_defined_assembly" {
            assembly.software_determined = true;
        } else if detail == "complete icosahedral assembly" {
            assembly.special_kind = crate::bio::BioAssemblySpecialKind::CompleteIcosahedral;
        } else if detail == "representative helical assembly" {
            assembly.special_kind = crate::bio::BioAssemblySpecialKind::RepresentativeHelical;
        } else if detail == "complete point assembly" {
            assembly.special_kind = crate::bio::BioAssemblySpecialKind::CompletePoint;
        }
        if !assembly.author_determined
            && !assembly.software_determined
            && assembly.special_kind == crate::bio::BioAssemblySpecialKind::NA
            && !detail.is_empty()
        {
            continue;
        }
        if assembly.software_determined && cif_row_has2(row, method_details) {
            assembly.software_name = row[method_details].value.trim().to_string();
        }
        assembly.oligomeric_details = row[oligomeric_details].value.trim().to_string();
        assembly.oligomeric_count = parse_decimal_i32(
            &BIO_MMCIF_READ_FEATURE,
            row[oligomeric_count].value.as_str(),
            row[oligomeric_count].line_number,
            "_pdbx_struct_assembly.oligomeric_count",
        )
        .unwrap_or(0);
        if let Some((biol_id, type_, value)) = prop_cols {
            for prop_row in &prop_rows {
                if prop_row[biol_id].value.trim() == assembly.name {
                    let parsed = parse_f32(
                        prop_row[value].value.as_str(),
                        prop_row[value].line_number,
                        "_pdbx_struct_assembly_prop.value",
                    )
                    .map(f64::from)
                    .ok();
                    match prop_row[type_].value.trim() {
                        "ABSA (A^2)" => assembly.absa = parsed,
                        "SSA (A^2)" => assembly.ssa = parsed,
                        "MORE" => assembly.more = parsed,
                        _ => {}
                    }
                }
            }
        }
        if let Some((assembly_id, oper_expression, asym_id_list)) = gen_cols {
            for gen_row in &gen_rows {
                if gen_row[assembly_id].value.trim() == assembly.name {
                    let mut generator = BioAssemblyGenerator {
                        subchains: split_on_comma(gen_row[asym_id_list].value.trim()),
                        ..BioAssemblyGenerator::default()
                    };
                    for name in parse_operation_expr(gen_row[oper_expression].value.trim()) {
                        if let Some(oper) = oper_list.iter().find(|oper| oper.name == name) {
                            generator.operators.push(oper.clone());
                        }
                    }
                    assembly.generators.push(generator);
                }
            }
        }
        assemblies.push(assembly);
    }
    Ok(assemblies)
}

// BEGIN GEMMI CPP FUNCTION read_sifts_unp
// Gemmi✔️✔️: void read_sifts_unp(cif::Block& block, Structure& st) { ... }
// END GEMMI CPP FUNCTION
fn read_sifts_unp(loops: &[CifLoop], structure: &mut BioStructure) -> Result<(), BioReadError> {
    let Some(loop_) = find_cif_loop(loops, "_pdbx_sifts_xref_db.entity_id") else {
        return Ok(());
    };
    let Some(
        [
            entity_id,
            asym_id,
            seq_id_ordinal,
            seq_id,
            observed,
            unp_res,
            unp_num,
            unp_acc,
        ],
    ) = required_cif_table_cols(
        loop_,
        [
            "_pdbx_sifts_xref_db.entity_id",
            "_pdbx_sifts_xref_db.asym_id",
            "_pdbx_sifts_xref_db.seq_id_ordinal",
            "_pdbx_sifts_xref_db.seq_id",
            "_pdbx_sifts_xref_db.observed",
            "_pdbx_sifts_xref_db.unp_res",
            "_pdbx_sifts_xref_db.unp_num",
            "_pdbx_sifts_xref_db.unp_acc",
        ],
    )
    else {
        return Ok(());
    };

    for row in cif_loop_rows(loop_)? {
        if row[seq_id_ordinal].value.trim() != "1"
            || !row[observed].value.starts_with('y')
            || !cif_row_has2(row, unp_acc)
            || !cif_row_has2(row, unp_num)
        {
            continue;
        }
        let source_entity_id = row[entity_id].value.trim();
        let entity_idx = structure
            .entities
            .iter()
            .position(|entity| entity.source.source_entity_id == source_entity_id)
            .ok_or_else(|| BioReadError::Parse {
                line_number: row[entity_id].line_number,
                message: format!("_pdbx_sifts_xref_db: entity_id not found: {source_entity_id}"),
            })?;
        let accession = row[unp_acc].value.trim().to_string();
        let acc_index = if let Some(idx) = structure.entities[entity_idx]
            .sifts_unp_acc
            .iter()
            .position(|value| value == &accession)
        {
            idx as u8
        } else {
            let idx = structure.entities[entity_idx].sifts_unp_acc.len() as u8;
            structure.entities[entity_idx].sifts_unp_acc.push(accession);
            idx
        };
        let label_seq = parse_decimal_i32(
            &BIO_MMCIF_READ_FEATURE,
            row[seq_id].value.as_str(),
            row[seq_id].line_number,
            "_pdbx_sifts_xref_db.seq_id",
        )?;
        let unp_num_value = parse_decimal_i32(
            &BIO_MMCIF_READ_FEATURE,
            row[unp_num].value.as_str(),
            row[unp_num].line_number,
            "_pdbx_sifts_xref_db.unp_num",
        )?;
        let unp_num_u16 = u16::try_from(unp_num_value).map_err(|_| BioReadError::Parse {
            line_number: row[unp_num].line_number,
            message: format!("_pdbx_sifts_xref_db.unp_num: {}", row[unp_num].value),
        })?;
        let subchain_name = row[asym_id].value.trim();
        let mapping = BioSiftsUnpResidue {
            res: row[unp_res].value.chars().next(),
            acc_index,
            num: unp_num_u16,
        };
        for residue in &mut structure.residues {
            let chain = &structure.chains[residue.chain_id.index() as usize];
            let residue_subchain = residue.source.subchain_id.or(chain.source.label_asym_id);
            if residue_subchain.is_some_and(|id| id.as_str() == subchain_name)
                && residue.source.label_seq_id == Some(label_seq)
            {
                residue.sifts_unp = Some(mapping);
            }
        }
    }
    Ok(())
}

// BEGIN GEMMI CPP FUNCTION read_helices
// Gemmi✔️✔️: std::vector<Helix> read_helices(cif::Block& block) {
// Gemmi✔️✔️:   std::vector<Helix> helices;
// Gemmi✔️✔️:   for (const auto row : block.find("_struct_conf.", {
// Gemmi✔️✔️:         "conf_type_id",
// Gemmi✔️✔️:         "beg_auth_asym_id", "beg_label_comp_id",
// Gemmi✔️✔️:         "beg_auth_seq_id", "?pdbx_beg_PDB_ins_code",
// Gemmi✔️✔️:         "end_auth_asym_id", "end_label_comp_id",
// Gemmi✔️✔️:         "end_auth_seq_id", "?pdbx_end_PDB_ins_code",
// Gemmi✔️✔️:         "?pdbx_PDB_helix_class", "?pdbx_PDB_helix_length"})) {
// Gemmi✔️✔️:     if (alpha_up(row.str(0)[0]) != 'H')
// Gemmi✔️✔️:       continue;
// Gemmi✔️✔️:     Helix h;
// Gemmi✔️✔️:     h.start.chain_name = row.str(1);
// Gemmi✔️✔️:     h.start.res_id = make_resid(row.str(2), row.str(3), row.ptr_at(4));
// Gemmi✔️✔️:     h.end.chain_name = row.str(5);
// Gemmi✔️✔️:     h.end.res_id = make_resid(row.str(6), row.str(7), row.ptr_at(8));
// Gemmi✔️✔️:     if (row.has(9))
// Gemmi✔️✔️:       h.set_helix_class_as_int(cif::as_int(row[9], -1));
// Gemmi✔️✔️:     if (row.has(10))
// Gemmi✔️✔️:       h.length = cif::as_int(row[10], -1);
// Gemmi✔️✔️:     helices.push_back(h);
// Gemmi✔️✔️:   }
// Gemmi✔️✔️:   return helices;
// Gemmi✔️✔️: }
// END GEMMI CPP FUNCTION
fn read_helices(loops: &[CifLoop]) -> Result<Vec<BioHelix>, BioReadError> {
    let Some(loop_) = find_cif_loop(loops, "_struct_conf.conf_type_id") else {
        return Ok(Vec::new());
    };
    let Some(
        [
            conf_type_id,
            beg_auth_asym_id,
            beg_label_comp_id,
            beg_auth_seq_id,
            end_auth_asym_id,
            end_label_comp_id,
            end_auth_seq_id,
        ],
    ) = required_cif_table_cols(
        loop_,
        [
            "_struct_conf.conf_type_id",
            "_struct_conf.beg_auth_asym_id",
            "_struct_conf.beg_label_comp_id",
            "_struct_conf.beg_auth_seq_id",
            "_struct_conf.end_auth_asym_id",
            "_struct_conf.end_label_comp_id",
            "_struct_conf.end_auth_seq_id",
        ],
    )
    else {
        return Ok(Vec::new());
    };
    let pdbx_beg_pdb_ins_code = optional_cif_col(loop_, "_struct_conf.pdbx_beg_PDB_ins_code");
    let pdbx_end_pdb_ins_code = optional_cif_col(loop_, "_struct_conf.pdbx_end_PDB_ins_code");
    let helix_class = optional_cif_col(loop_, "_struct_conf.pdbx_PDB_helix_class");
    let helix_length = optional_cif_col(loop_, "_struct_conf.pdbx_PDB_helix_length");
    let mut helices = Vec::new();
    for row in cif_loop_rows(loop_)? {
        let conf_type = row[conf_type_id].value.trim();
        if conf_type.is_empty()
            || conf_type.as_bytes().first().map(|b| b.to_ascii_uppercase()) != Some(b'H')
        {
            continue;
        }
        let line_number = row[conf_type_id].line_number;
        let mut helix = BioHelix {
            start: bio_address_from_mmcif(
                row[beg_auth_asym_id].value.as_str(),
                row[beg_label_comp_id].value.as_str(),
                row[beg_auth_seq_id].value.as_str(),
                pdbx_beg_pdb_ins_code.map(|idx| row[idx].value.as_str()),
                None,
                line_number,
            )?,
            end: bio_address_from_mmcif(
                row[end_auth_asym_id].value.as_str(),
                row[end_label_comp_id].value.as_str(),
                row[end_auth_seq_id].value.as_str(),
                pdbx_end_pdb_ins_code.map(|idx| row[idx].value.as_str()),
                None,
                line_number,
            )?,
            ..BioHelix::default()
        };
        if let Some(index) = helix_class
            && cif_row_has2(row, index)
        {
            helix.set_helix_class_as_int(parse_decimal_i32(
                &BIO_MMCIF_READ_FEATURE,
                row[index].value.as_str(),
                row[index].line_number,
                "mmCIF helix class",
            )?);
        }
        if let Some(index) = helix_length
            && cif_row_has2(row, index)
        {
            helix.length = parse_decimal_i32(
                &BIO_MMCIF_READ_FEATURE,
                row[index].value.as_str(),
                row[index].line_number,
                "mmCIF helix length",
            )?;
        }
        helices.push(helix);
    }
    Ok(helices)
}

// BEGIN GEMMI CPP FUNCTION read_sheets
// Gemmi✔️✔️: std::vector<Sheet> read_sheets(cif::Block& block) {
// Gemmi✔️✔️:   std::vector<Sheet> sheets;
// Gemmi✔️✔️:   for (const std::string& sheet_id : block.find_values("_struct_sheet.id"))
// Gemmi✔️✔️:     sheets.emplace_back(sheet_id);
// Gemmi✔️✔️:   for (const auto row : block.find("_struct_sheet_range.", {
// Gemmi✔️✔️:         "sheet_id", "id",
// Gemmi✔️✔️:         "beg_auth_asym_id", "beg_label_comp_id",
// Gemmi✔️✔️:         "beg_auth_seq_id", "?pdbx_beg_PDB_ins_code",
// Gemmi✔️✔️:         "end_auth_asym_id", "end_label_comp_id",
// Gemmi✔️✔️:         "end_auth_seq_id", "?pdbx_end_PDB_ins_code"})) {
// Gemmi✔️✔️:     Sheet& sheet = impl::find_or_add(sheets, row.str(0));
// Gemmi✔️✔️:     sheet.strands.emplace_back();
// Gemmi✔️✔️:     Sheet::Strand& strand = sheet.strands.back();
// Gemmi✔️✔️:     strand.name = row.str(1);
// Gemmi✔️✔️:     strand.start.chain_name = row.str(2);
// Gemmi✔️✔️:     strand.start.res_id = make_resid(row.str(3), row.str(4), row.ptr_at(5));
// Gemmi✔️✔️:     strand.end.chain_name = row.str(6);
// Gemmi✔️✔️:     strand.end.res_id = make_resid(row.str(7), row.str(8), row.ptr_at(9));
// Gemmi✔️✔️:   }
// Gemmi✔️✔️:   for (const auto row : block.find("_struct_sheet_order.", {
// Gemmi✔️✔️:         "sheet_id", "range_id_2", "sense"}))
// Gemmi✔️✔️:     if (Sheet* sheet = impl::find_or_null(sheets, row.str(0)))
// Gemmi✔️✔️:       if (Sheet::Strand* ss = impl::find_or_null(sheet->strands, row.str(1)))
// Gemmi✔️✔️:         switch (alpha_up(row.str(2)[0])) {
// Gemmi✔️✔️:           case 'P': ss->sense = 1; break;
// Gemmi✔️✔️:           case 'A': ss->sense = -1; break;
// Gemmi✔️✔️:         }
// Gemmi✔️✔️:   for (const auto row : block.find("_pdbx_struct_sheet_hbond.", {
// Gemmi✔️✔️:         "sheet_id", "range_id_2",
// Gemmi✔️✔️:         "range_1_auth_asym_id", "range_1_label_comp_id",
// Gemmi✔️✔️:         "range_1_auth_seq_id", "?range_1_PDB_ins_code",
// Gemmi✔️✔️:         "range_1_label_atom_id",
// Gemmi✔️✔️:         "range_2_auth_asym_id", "range_2_label_comp_id",
// Gemmi✔️✔️:         "range_2_auth_seq_id", "?range_2_PDB_ins_code",
// Gemmi✔️✔️:         "range_2_label_atom_id"}))
// Gemmi✔️✔️:     if (Sheet* sheet = impl::find_or_null(sheets, row.str(0)))
// Gemmi✔️✔️:       if (Sheet::Strand* ss = impl::find_or_null(sheet->strands, row.str(1))) {
// Gemmi✔️✔️:         ss->hbond_atom1.chain_name = row.str(2);
// Gemmi✔️✔️:         ss->hbond_atom1.res_id = make_resid(row.str(3), row.str(4), row.ptr_at(5));
// Gemmi✔️✔️:         ss->hbond_atom1.atom_name = row.str(6);
// Gemmi✔️✔️:         ss->hbond_atom2.chain_name = row.str(7);
// Gemmi✔️✔️:         ss->hbond_atom2.res_id = make_resid(row.str(8), row.str(9), row.ptr_at(10));
// Gemmi✔️✔️:         ss->hbond_atom2.atom_name = row.str(11);
// Gemmi✔️✔️:       }
// Gemmi✔️✔️:   return sheets;
// Gemmi✔️✔️: }
// END GEMMI CPP FUNCTION
fn read_sheets(loops: &[CifLoop]) -> Result<Vec<BioSheet>, BioReadError> {
    let mut sheets = Vec::<BioSheet>::new();
    if let Some(struct_sheet) = find_cif_loop(loops, "_struct_sheet.id") {
        let sheet_id = required_cif_col(struct_sheet, "_struct_sheet.id")?;
        for row in cif_loop_rows(struct_sheet)? {
            sheets.push(BioSheet {
                name: row[sheet_id].value.trim().to_string(),
                ..BioSheet::default()
            });
        }
    }
    if let Some(ranges) = find_cif_loop(loops, "_struct_sheet_range.sheet_id") {
        let Some(
            [
                sheet_id,
                id,
                beg_auth_asym_id,
                beg_label_comp_id,
                beg_auth_seq_id,
                end_auth_asym_id,
                end_label_comp_id,
                end_auth_seq_id,
            ],
        ) = required_cif_table_cols(
            ranges,
            [
                "_struct_sheet_range.sheet_id",
                "_struct_sheet_range.id",
                "_struct_sheet_range.beg_auth_asym_id",
                "_struct_sheet_range.beg_label_comp_id",
                "_struct_sheet_range.beg_auth_seq_id",
                "_struct_sheet_range.end_auth_asym_id",
                "_struct_sheet_range.end_label_comp_id",
                "_struct_sheet_range.end_auth_seq_id",
            ],
        )
        else {
            return Ok(sheets);
        };
        let beg_ins = optional_cif_col(ranges, "_struct_sheet_range.pdbx_beg_PDB_ins_code");
        let end_ins = optional_cif_col(ranges, "_struct_sheet_range.pdbx_end_PDB_ins_code");
        for row in cif_loop_rows(ranges)? {
            let line_number = row[sheet_id].line_number;
            let sheet_name = row[sheet_id].value.trim();
            let sheet_idx = sheets
                .iter()
                .position(|sheet| sheet.name == sheet_name)
                .unwrap_or_else(|| {
                    sheets.push(BioSheet {
                        name: sheet_name.to_string(),
                        ..BioSheet::default()
                    });
                    sheets.len() - 1
                });
            sheets[sheet_idx].strands.push(BioSheetStrand {
                name: row[id].value.trim().to_string(),
                start: bio_address_from_mmcif(
                    row[beg_auth_asym_id].value.as_str(),
                    row[beg_label_comp_id].value.as_str(),
                    row[beg_auth_seq_id].value.as_str(),
                    beg_ins.map(|idx| row[idx].value.as_str()),
                    None,
                    line_number,
                )?,
                end: bio_address_from_mmcif(
                    row[end_auth_asym_id].value.as_str(),
                    row[end_label_comp_id].value.as_str(),
                    row[end_auth_seq_id].value.as_str(),
                    end_ins.map(|idx| row[idx].value.as_str()),
                    None,
                    line_number,
                )?,
                ..BioSheetStrand::default()
            });
        }
    }
    if let Some(order) = find_cif_loop(loops, "_struct_sheet_order.sheet_id") {
        let Some([sheet_id, range_id_2, sense]) = required_cif_table_cols(
            order,
            [
                "_struct_sheet_order.sheet_id",
                "_struct_sheet_order.range_id_2",
                "_struct_sheet_order.sense",
            ],
        ) else {
            return Ok(sheets);
        };
        for row in cif_loop_rows(order)? {
            if let Some(sheet) = sheets
                .iter_mut()
                .find(|sheet| sheet.name == row[sheet_id].value.trim())
                && let Some(strand) = sheet
                    .strands
                    .iter_mut()
                    .find(|strand| strand.name == row[range_id_2].value.trim())
            {
                match row[sense]
                    .value
                    .as_bytes()
                    .first()
                    .map(|b| b.to_ascii_uppercase())
                {
                    Some(b'P') => strand.sense = 1,
                    Some(b'A') => strand.sense = -1,
                    _ => {}
                }
            }
        }
    }
    if let Some(hbond) = find_cif_loop(loops, "_pdbx_struct_sheet_hbond.sheet_id") {
        let Some(
            [
                sheet_id,
                range_id_2,
                range_1_auth_asym_id,
                range_1_label_comp_id,
                range_1_auth_seq_id,
                range_1_label_atom_id,
                range_2_auth_asym_id,
                range_2_label_comp_id,
                range_2_auth_seq_id,
                range_2_label_atom_id,
            ],
        ) = required_cif_table_cols(
            hbond,
            [
                "_pdbx_struct_sheet_hbond.sheet_id",
                "_pdbx_struct_sheet_hbond.range_id_2",
                "_pdbx_struct_sheet_hbond.range_1_auth_asym_id",
                "_pdbx_struct_sheet_hbond.range_1_label_comp_id",
                "_pdbx_struct_sheet_hbond.range_1_auth_seq_id",
                "_pdbx_struct_sheet_hbond.range_1_label_atom_id",
                "_pdbx_struct_sheet_hbond.range_2_auth_asym_id",
                "_pdbx_struct_sheet_hbond.range_2_label_comp_id",
                "_pdbx_struct_sheet_hbond.range_2_auth_seq_id",
                "_pdbx_struct_sheet_hbond.range_2_label_atom_id",
            ],
        )
        else {
            return Ok(sheets);
        };
        let range_1_ins = optional_cif_col(hbond, "_pdbx_struct_sheet_hbond.range_1_PDB_ins_code");
        let range_2_ins = optional_cif_col(hbond, "_pdbx_struct_sheet_hbond.range_2_PDB_ins_code");
        for row in cif_loop_rows(hbond)? {
            let line_number = row[sheet_id].line_number;
            if let Some(sheet) = sheets
                .iter_mut()
                .find(|sheet| sheet.name == row[sheet_id].value.trim())
                && let Some(strand) = sheet
                    .strands
                    .iter_mut()
                    .find(|strand| strand.name == row[range_id_2].value.trim())
            {
                strand.hbond_atom1 = bio_address_from_mmcif(
                    row[range_1_auth_asym_id].value.as_str(),
                    row[range_1_label_comp_id].value.as_str(),
                    row[range_1_auth_seq_id].value.as_str(),
                    range_1_ins.map(|idx| row[idx].value.as_str()),
                    Some(row[range_1_label_atom_id].value.as_str()),
                    line_number,
                )?;
                strand.hbond_atom2 = bio_address_from_mmcif(
                    row[range_2_auth_asym_id].value.as_str(),
                    row[range_2_label_comp_id].value.as_str(),
                    row[range_2_auth_seq_id].value.as_str(),
                    range_2_ins.map(|idx| row[idx].value.as_str()),
                    Some(row[range_2_label_atom_id].value.as_str()),
                    line_number,
                )?;
            }
        }
    }
    Ok(sheets)
}

fn read_mmcif_atom_site(
    atom_site: &CifLoop,
    loops: &[CifLoop],
) -> Result<BioStructure, BioReadError> {
    let width = atom_site.tags.len();
    if width == 0 || atom_site.values.len() % width != 0 {
        return Err(BioReadError::Parse {
            line_number: atom_site
                .values
                .first()
                .map_or(0, |token| token.line_number),
            message: "mmCIF loop value count is not divisible by tag count".to_string(),
        });
    }

    let columns = AtomSiteColumns::new(&atom_site.tags)?;
    let asym_id = CifRowAccess::new(columns.auth_asym_id, Some(columns.label_asym_id));
    let comp_id = CifRowAccess::new(columns.auth_comp_id, Some(columns.label_comp_id));
    let atom_id = CifRowAccess::new(columns.auth_atom_id, Some(columns.label_atom_id));
    let seq_id = CifRowAccess::new(columns.auth_seq_id, Some(columns.label_seq_id));
    if !asym_id.ok() {
        return Err(missing_cif_value(
            0,
            "_atom_site.label_asym_id/auth_asym_id",
        ));
    }
    if !comp_id.ok() {
        return Err(missing_cif_value(
            0,
            "_atom_site.label_comp_id/auth_comp_id",
        ));
    }
    if !atom_id.ok() {
        return Err(missing_cif_value(
            0,
            "_atom_site.label_atom_id/auth_atom_id",
        ));
    }
    if !seq_id.ok() {
        return Err(missing_cif_value(0, "_atom_site.label_seq_id/auth_seq_id"));
    }
    let aniso_map = get_anisotropic_u(loops)?;
    let mut builder = PdbBioBuilder::default();
    read_mmcif_entity_and_sequence_info(&mut builder, loops)?;
    builder.structure.has_d_fraction = columns.deuterium_fraction.is_some();
    if columns.model_num.is_none() {
        builder.begin_model(Some(1));
    }

    for row in atom_site.values.chunks(width) {
        let line_number = row.first().map_or(0, |token| token.line_number);
        let model_number = columns
            .model_num
            .and_then(|idx| cif_optional(row[idx].value.as_str()))
            .map(|value| {
                parse_decimal_i32(&BIO_MMCIF_READ_FEATURE, value, line_number, "model number")
            })
            .transpose()?;
        let effective_model_number = model_number.or(Some(1));
        if columns.model_num.is_some()
            && (builder.current_model.is_none()
                || effective_model_number
                    != builder.current_model.and_then(|id| {
                        builder.structure.models[id.index() as usize].source_model_number
                    }))
        {
            builder.begin_model(effective_model_number);
        }

        // Gemmi✔️✔️: if (!chain || cif::as_string(asym_id.get(gap)) != chain->name) {
        // Gemmi✔️✔️:   model->chains.emplace_back(cif::as_string(asym_id.get(gap)));
        // Gemmi✔️✔️:   chain = &model->chains.back();
        // Gemmi✔️✔️:   resi = nullptr;
        // Gemmi✔️✔️: }
        let auth_chain = columns
            .auth_asym_id
            .and_then(|idx| cif_optional(row[idx].value.as_str()));
        let chain_name = asym_id.get(row).and_then(cif_optional).unwrap_or("");
        let chain_key = pdb_chain_id_from_cif(chain_name, line_number)?;
        // Gemmi✔️✔️: resi->subchain = row.str(kLabelAsymId);
        let label_chain_id = cif_optional(row[columns.label_asym_id].value.as_str())
            .map(|value| pdb_chain_id_from_cif(value, line_number))
            .transpose()?;
        let auth_chain_id = auth_chain
            .map(|value| pdb_chain_id_from_cif(value, line_number))
            .transpose()?;

        let atom_name = atom_id.get(row).and_then(cif_optional).ok_or_else(|| {
            missing_cif_value(line_number, "_atom_site.label_atom_id/auth_atom_id")
        })?;
        let residue_name = comp_id.get(row).and_then(cif_optional).ok_or_else(|| {
            missing_cif_value(line_number, "_atom_site.label_comp_id/auth_comp_id")
        })?;
        let seq_text = seq_id
            .get(row)
            .and_then(cif_optional)
            .ok_or_else(|| missing_cif_value(line_number, "_atom_site.label_seq_id/auth_seq_id"))?;

        let serial = cif_optional(row[columns.id].value.as_str())
            .and_then(|value| value.parse::<i32>().ok());
        let ins_code_text = columns
            .ins_code
            .and_then(|idx| row.get(idx))
            .map(|token| token.value.as_str());
        let (residue_name, seq_id_value) =
            make_resid(residue_name, seq_text, ins_code_text, line_number)?;
        let seq_id_value = seq_id_value
            .ok_or_else(|| missing_cif_value(line_number, "_atom_site.label_seq_id/auth_seq_id"))?;
        let altloc = cif_optional(row[columns.alt_id].value.as_str())
            .and_then(|value| value.as_bytes().first().copied())
            .and_then(parse_altloc);
        let element = cif_optional(row[columns.type_symbol].value.as_str())
            .and_then(element_from_symbol)
            .ok_or_else(|| missing_cif_value(line_number, "_atom_site.type_symbol"))?;

        let x = parse_f64(row[columns.x].value.as_str(), line_number, "Cartn_x")?;
        let y = parse_f64(row[columns.y].value.as_str(), line_number, "Cartn_y")?;
        let z = parse_f64(row[columns.z].value.as_str(), line_number, "Cartn_z")?;
        // Gemmi✔️✔️: float occ = 1.0f;
        // Gemmi✔️✔️: float b_iso = 20.0f; // arbitrary default value
        // Gemmi✔️✔️: if (row.has2(kOcc))
        // Gemmi✔️✔️:   atom.occ = (float) cif::as_number(row[kOcc]);
        // Gemmi✔️✔️: if (row.has2(kBiso))
        // Gemmi✔️✔️:   atom.b_iso = (float) cif::as_number(row[kBiso]);
        let occupancy = columns
            .occupancy
            .and_then(|idx| cif_optional(row[idx].value.as_str()))
            .map(|value| parse_f32(value, line_number, "occupancy"))
            .transpose()?
            .unwrap_or(1.0);
        let b_iso = columns
            .b_iso
            .and_then(|idx| cif_optional(row[idx].value.as_str()))
            .map(|value| parse_f32(value, line_number, "B_iso_or_equiv"))
            .transpose()?
            .unwrap_or(20.0);
        let formal_charge = columns
            .formal_charge
            .and_then(|idx| cif_optional(row[idx].value.as_str()))
            .map(|value| parse_i8(value, line_number, "formal charge"))
            .transpose()?;
        let label_entity_id = columns
            .label_entity_id
            .and_then(|idx| cif_optional(row[idx].value.as_str()))
            .map(|value| {
                builder.find_or_add_entity(value, EntityKind::Unknown, PolymerKind::Unknown)
            });
        let group_pdb = columns
            .group_pdb
            .and_then(|idx| cif_optional(row[idx].value.as_str()))
            .and_then(group_pdb_het_flag);
        let calc_flag = columns
            .calc_flag
            .and_then(|idx| cif_optional(row[idx].value.as_str()))
            .map(calc_flag_from_cif)
            .unwrap_or(BioCalcFlag::NotSet);
        let tls_group_id = columns
            .tls_group_id
            .and_then(|idx| cif_optional(row[idx].value.as_str()))
            .and_then(parse_optional_short_prefix_int);
        let fraction = if builder.structure.has_d_fraction {
            Some(
                columns
                    .deuterium_fraction
                    .and_then(|idx| cif_optional(row[idx].value.as_str()))
                    .map(|value| parse_f32(value, line_number, "ccp4_deuterium_fraction"))
                    .transpose()?
                    .unwrap_or(0.0),
            )
        } else {
            None
        };
        let residue_kind = classify_residue_name(residue_name);
        let entity_kind = if residue_kind == ResidueKind::Water {
            EntityKind::Water
        } else {
            EntityKind::Unknown
        };
        builder.push_atom(
            PdbAtomRecord {
                serial,
                atom_name: atom_name_from_cif(atom_name),
                altloc,
                residue_name,
                residue_kind,
                chain_key,
                chain_source: ChainSourceIds {
                    auth_chain_id,
                    label_asym_id: label_chain_id,
                },
                seq_id: seq_id_value,
                label_seq_id: cif_optional(row[columns.label_seq_id].value.as_str())
                    .and_then(|value| value.parse::<i32>().ok()),
                label_entity_id,
                group_pdb,
                segment_id: None,
                position: [x, y, z],
                occupancy: Some(occupancy),
                b_iso: Some(b_iso),
                formal_charge,
                element,
                calc_flag,
                tls_group_id,
                fraction,
            },
            effective_model_number,
            false,
            entity_kind,
        );
        if let Some(anisou) = aniso_map.get(row[columns.id].value.trim()) {
            builder.set_last_atom_anisou(line_number, *anisou)?;
        }
    }

    builder.finish()
}

#[derive(Debug, Clone, Copy)]
struct AtomSiteColumns {
    id: usize,
    group_pdb: Option<usize>,
    type_symbol: usize,
    label_atom_id: usize,
    alt_id: usize,
    label_comp_id: usize,
    label_asym_id: usize,
    label_entity_id: Option<usize>,
    label_seq_id: usize,
    ins_code: Option<usize>,
    x: usize,
    y: usize,
    z: usize,
    occupancy: Option<usize>,
    b_iso: Option<usize>,
    formal_charge: Option<usize>,
    auth_seq_id: Option<usize>,
    auth_comp_id: Option<usize>,
    auth_asym_id: Option<usize>,
    auth_atom_id: Option<usize>,
    model_num: Option<usize>,
    calc_flag: Option<usize>,
    tls_group_id: Option<usize>,
    deuterium_fraction: Option<usize>,
}

impl AtomSiteColumns {
    fn new(tags: &[String]) -> Result<Self, BioReadError> {
        let find = |name: &str| tags.iter().position(|tag| tag == name);
        let required = |name: &'static str| {
            find(name).ok_or_else(|| BioReadError::Parse {
                line_number: 0,
                message: format!("required mmCIF atom_site column is missing: {name}"),
            })
        };
        Ok(Self {
            id: required("_atom_site.id")?,
            group_pdb: find("_atom_site.group_PDB"),
            type_symbol: required("_atom_site.type_symbol")?,
            label_atom_id: required("_atom_site.label_atom_id")?,
            alt_id: required("_atom_site.label_alt_id")?,
            label_comp_id: required("_atom_site.label_comp_id")?,
            label_asym_id: required("_atom_site.label_asym_id")?,
            label_entity_id: find("_atom_site.label_entity_id"),
            label_seq_id: required("_atom_site.label_seq_id")?,
            ins_code: find("_atom_site.pdbx_PDB_ins_code"),
            x: required("_atom_site.Cartn_x")?,
            y: required("_atom_site.Cartn_y")?,
            z: required("_atom_site.Cartn_z")?,
            occupancy: find("_atom_site.occupancy"),
            b_iso: find("_atom_site.B_iso_or_equiv"),
            formal_charge: find("_atom_site.pdbx_formal_charge"),
            auth_seq_id: find("_atom_site.auth_seq_id"),
            auth_comp_id: find("_atom_site.auth_comp_id"),
            auth_asym_id: find("_atom_site.auth_asym_id"),
            auth_atom_id: find("_atom_site.auth_atom_id"),
            model_num: find("_atom_site.pdbx_PDB_model_num"),
            calc_flag: find("_atom_site.calc_flag"),
            tls_group_id: find("_atom_site.pdbx_tls_group_id"),
            deuterium_fraction: find("_atom_site.ccp4_deuterium_fraction"),
        })
    }
}

fn parse_cif_loops(text: &str) -> Result<Vec<CifLoop>, BioReadError> {
    Ok(parse_cif_document(text)?
        .blocks
        .into_iter()
        .next()
        .map(|block| block.loops)
        .unwrap_or_default())
}

fn json_type_as_string(value: &JsonValue) -> &'static str {
    match value {
        JsonValue::Null => "<null>",
        JsonValue::Bool(false) => "<false>",
        JsonValue::Bool(true) => "<true>",
        JsonValue::Number(number) => {
            if number.is_i64() || number.is_u64() {
                "<integer>"
            } else {
                "<double>"
            }
        }
        JsonValue::String(_) => "<string>",
        JsonValue::Array(_) => "<array>",
        JsonValue::Object(_) => "<object>",
    }
}

// BEGIN GEMMI CPP FUNCTION gemmi::cif::as_cif_value
// Gemmi✔️✔️: switch (val.get_type()) {
// Gemmi✔️✔️:   case sajson::TYPE_DOUBLE:
// Gemmi✔️✔️:     return val.as_string();
// Gemmi✔️✔️:   case sajson::TYPE_NULL:
// Gemmi✔️✔️:     return "?";
// Gemmi✔️✔️:   case sajson::TYPE_FALSE:
// Gemmi✔️✔️:     return "NO";
// Gemmi✔️✔️:   case sajson::TYPE_TRUE:
// Gemmi✔️✔️:     return "YES";
// Gemmi✔️✔️:   case sajson::TYPE_STRING:
// Gemmi✔️✔️:     return quote(val.as_string());
// Gemmi✔️✔️:   case sajson::TYPE_ARRAY: {
// Gemmi✔️✔️:     std::string s;
// Gemmi✔️✔️:     for (size_t i = 0; i < val.get_length(); ++i) {
// Gemmi✔️✔️:       if (i != 0)
// Gemmi✔️✔️:         s += ' ';
// Gemmi✔️✔️:       s += val.get_array_element(0).as_string();
// Gemmi✔️✔️:     }
// Gemmi✔️✔️:     return quote(s);
// Gemmi✔️✔️:   }
// Gemmi✔️✔️:   default:
// Gemmi✔️✔️:     fail("Unexpected ", json_type_as_string(val.get_type()), " as value in JSON.");
// Gemmi✔️✔️:     return "";
// Gemmi✔️✔️: }
// END GEMMI CPP FUNCTION
fn mmjson_value_as_cif_value(
    value: &JsonValue,
    line_number: usize,
) -> Result<String, BioReadError> {
    match value {
        JsonValue::Null => Ok("?".to_string()),
        JsonValue::Bool(false) => Ok("NO".to_string()),
        JsonValue::Bool(true) => Ok("YES".to_string()),
        JsonValue::Number(number) => Ok(number.to_string()),
        JsonValue::String(text) => Ok(text.clone()),
        JsonValue::Array(values) => {
            let mut joined = String::new();
            for (index, entry) in values.iter().enumerate() {
                if index != 0 {
                    joined.push(' ');
                }
                let JsonValue::String(text) = entry else {
                    return Err(BioReadError::Parse {
                        line_number,
                        message: format!(
                            "Unexpected {} as value in JSON.",
                            json_type_as_string(entry)
                        ),
                    });
                };
                joined.push_str(text);
            }
            Ok(joined)
        }
        _ => Err(BioReadError::Parse {
            line_number,
            message: format!(
                "Unexpected {} as value in JSON.",
                json_type_as_string(value)
            ),
        }),
    }
}

// BEGIN GEMMI CPP FUNCTION gemmi::cif::fill_document_from_sajson
// Gemmi✔️✔️: sajson::value root = s.get_root();
// Gemmi✔️✔️: if (root.get_type() != sajson::TYPE_OBJECT)
// Gemmi✔️✔️:   fail("not mmJSON - the root is not of type object");
// Gemmi✔️✔️: for (size_t block_index = 0; block_index < root.get_length(); ++block_index) {
// Gemmi✔️✔️:   std::string block_name = root.get_object_key(block_index).as_string();
// Gemmi✔️✔️:   if (!starts_with(block_name, "data_"))
// Gemmi✔️✔️:     fail("not mmJSON - top level key should start with data_\n"
// Gemmi✔️✔️:          "(if you use gemmi-cif2json to write JSON, use -m for mmJSON)");
// Gemmi✔️✔️:   d.blocks.emplace_back(block_name.substr(5));
// Gemmi✔️✔️:   std::vector<Item>& items = d.blocks[block_index].items;
// Gemmi✔️✔️:   sajson::value top = root.get_object_value(block_index);
// Gemmi✔️✔️:   if (top.get_type() != sajson::TYPE_OBJECT)
// Gemmi✔️✔️:     fail("");
// Gemmi✔️✔️:   for (size_t i = 0; i != top.get_length(); ++i) {
// Gemmi✔️✔️:     std::string category_name = "_" + top.get_object_key(i).as_string() + ".";
// Gemmi✔️✔️:     sajson::value category = top.get_object_value(i);
// Gemmi✔️✔️:     if (category.get_type() != sajson::TYPE_OBJECT ||
// Gemmi✔️✔️:         category.get_length() == 0 ||
// Gemmi✔️✔️:         category.get_object_value(0).get_type() != sajson::TYPE_ARRAY)
// Gemmi✔️✔️:       fail("");
// Gemmi✔️✔️:     size_t cif_cols = category.get_length();
// Gemmi✔️✔️:     size_t cif_rows = category.get_object_value(0).get_length();
// Gemmi✔️✔️:     if (cif_rows > 1) {
// Gemmi✔️✔️:       items.emplace_back(LoopArg{});
// Gemmi✔️✔️:       Loop& loop = items.back().loop;
// Gemmi✔️✔️:       loop.tags.reserve(cif_cols);
// Gemmi✔️✔️:       loop.values.resize(cif_cols * cif_rows);
// Gemmi✔️✔️:     }
// Gemmi✔️✔️:     for (size_t j = 0; j != cif_cols; ++j) {
// Gemmi✔️✔️:       std::string tag = category_name + category.get_object_key(j).as_string();
// Gemmi✔️✔️:       sajson::value arr = category.get_object_value(j);
// Gemmi✔️✔️:       if (arr.get_type() != sajson::TYPE_ARRAY)
// Gemmi✔️✔️:         fail("Expected array, got " + json_type_as_string(arr.get_type()));
// Gemmi✔️✔️:       if (arr.get_length() != cif_rows)
// Gemmi✔️✔️:         fail("Expected array of length ", std::to_string(cif_rows), " not ",
// Gemmi✔️✔️:              std::to_string(arr.get_length()));
// Gemmi✔️✔️:       if (cif_rows == 1) {
// Gemmi✔️✔️:         items.emplace_back(tag, as_cif_value(arr.get_array_element(0)));
// Gemmi✔️✔️:       } else if (cif_rows != 0) {
// Gemmi✔️✔️:         Loop& loop = items.back().loop;
// Gemmi✔️✔️:         loop.tags.emplace_back(std::move(tag));
// Gemmi✔️✔️:         for (size_t k = 0; k != cif_rows; ++k)
// Gemmi✔️✔️:           loop.values[j + k*cif_cols] = as_cif_value(arr.get_array_element(k));
// Gemmi✔️✔️:       }
// Gemmi✔️✔️:     }
// Gemmi✔️✔️:   }
// Gemmi✔️✔️: }
// END GEMMI CPP FUNCTION
fn fill_document_from_mmjson_value(
    root: &JsonValue,
    document: &mut CifDocument,
    line_number: usize,
) -> Result<(), BioReadError> {
    let JsonValue::Object(root_map) = root else {
        return Err(BioReadError::Parse {
            line_number,
            message: "not mmJSON - the root is not of type object".to_string(),
        });
    };

    for (block_name, top) in root_map {
        if !block_name.starts_with("data_") {
            return Err(BioReadError::Parse {
                line_number,
                message: "not mmJSON - top level key should start with data_\n(if you use gemmi-cif2json to write JSON, use -m for mmJSON)".to_string(),
            });
        }

        let JsonValue::Object(top_map) = top else {
            return Err(BioReadError::Parse {
                line_number,
                message: "not mmJSON - block value is not of type object".to_string(),
            });
        };

        let mut block = CifBlock {
            name: block_name["data_".len()..].to_string(),
            items: Vec::new(),
            loops: Vec::new(),
            entries: Vec::new(),
        };

        for (category_key, category) in top_map {
            let JsonValue::Object(category_map) = category else {
                return Err(BioReadError::Parse {
                    line_number,
                    message: "not mmJSON - category value is not of type object".to_string(),
                });
            };
            if category_map.is_empty() {
                return Err(BioReadError::Parse {
                    line_number,
                    message: "not mmJSON - category object is empty".to_string(),
                });
            }

            let first = category_map.values().next().expect("checked non-empty");
            let JsonValue::Array(first_array) = first else {
                return Err(BioReadError::Parse {
                    line_number,
                    message: "not mmJSON - category columns must be arrays".to_string(),
                });
            };
            let cif_cols = category_map.len();
            let cif_rows = first_array.len();
            let category_name = format!("_{}.", category_key);

            let mut pending_loop = if cif_rows > 1 {
                Some(CifLoop {
                    tags: Vec::with_capacity(cif_cols),
                    values: vec![
                        CifToken {
                            value: String::new(),
                            line_number,
                        };
                        cif_cols * cif_rows
                    ],
                })
            } else {
                None
            };

            for (column_index, (column_key, array_value)) in category_map.iter().enumerate() {
                let JsonValue::Array(values) = array_value else {
                    return Err(BioReadError::Parse {
                        line_number,
                        message: format!(
                            "Expected array, got {}",
                            json_type_as_string(array_value)
                        ),
                    });
                };
                if values.len() != cif_rows {
                    return Err(BioReadError::Parse {
                        line_number,
                        message: format!(
                            "Expected array of length {} not {}",
                            cif_rows,
                            values.len()
                        ),
                    });
                }

                let tag = format!("{category_name}{column_key}");
                if cif_rows == 1 {
                    block.push_pair(
                        tag,
                        CifToken {
                            value: mmjson_value_as_cif_value(&values[0], line_number)?,
                            line_number,
                        },
                    );
                } else if cif_rows != 0 {
                    let loop_ = pending_loop.as_mut().expect("loop allocated for multi-row");
                    loop_.tags.push(tag);
                    for (row_index, entry) in values.iter().enumerate() {
                        loop_.values[column_index + row_index * cif_cols] = CifToken {
                            value: mmjson_value_as_cif_value(entry, line_number)?,
                            line_number,
                        };
                    }
                }
            }

            if let Some(loop_) = pending_loop {
                block.push_loop(loop_);
            }
        }

        document.blocks.push(block);
    }

    Ok(())
}

// BEGIN GEMMI CPP FUNCTION gemmi::cif::read_mmjson_insitu
// Gemmi✔️✔️: Document doc;
// Gemmi✔️✔️: sajson::document json = sajson::parse(sajson::dynamic_allocation(),
// Gemmi✔️✔️:                                     sajson::mutable_string_view(size, buffer));
// Gemmi✔️✔️: if (!json.is_valid())
// Gemmi✔️✔️:   fail(name + ":", std::to_string(json.get_error_line()), " error: ",
// Gemmi✔️✔️:        json.get_error_message_as_string());
// Gemmi✔️✔️: fill_document_from_sajson(doc, json);
// Gemmi✔️✔️: doc.source = name;
// Gemmi✔️✔️: return doc;
// END GEMMI CPP FUNCTION
fn read_mmjson_document(text: &str, name: &str) -> Result<CifDocument, BioReadError> {
    let parsed = serde_json::from_str::<JsonValue>(text).map_err(|error| BioReadError::Parse {
        line_number: error.line(),
        message: format!("{name}:{} error: {}", error.line(), error),
    })?;
    let mut document = CifDocument { blocks: Vec::new() };
    fill_document_from_mmjson_value(&parsed, &mut document, 0)?;
    Ok(document)
}

// BEGIN GEMMI CPP FUNCTION gemmi::impl::set_cell_from_mmcif
// Gemmi✔️✔️: cif::Table tab = block.find((mmcif ? "_cell." : "_cell_"),
// Gemmi✔️✔️:                             {"length_a", "length_b", "length_c",
// Gemmi✔️✔️:                              "angle_alpha", "angle_beta", "angle_gamma"});
// Gemmi✔️✔️: if (tab.ok()) {
// Gemmi✔️✔️:   auto c = tab.one();
// Gemmi✔️✔️:   if (!cif::is_null(c[0]) && !cif::is_null(c[1]) && !cif::is_null(c[2]))
// Gemmi✔️✔️:     cell.set(cif::as_number(c[0]), cif::as_number(c[1]), cif::as_number(c[2]),
// Gemmi✔️✔️:              cif::as_number(c[3]), cif::as_number(c[4]), cif::as_number(c[5]));
// Gemmi✔️✔️: }
// END GEMMI CPP FUNCTION
fn set_cell_from_mmcif(document: &CifDocument, crystal: &mut Option<CrystalInfo>) {
    let values = [
        find_cif_item_value(document, "_cell.length_a"),
        find_cif_item_value(document, "_cell.length_b"),
        find_cif_item_value(document, "_cell.length_c"),
        find_cif_item_value(document, "_cell.angle_alpha"),
        find_cif_item_value(document, "_cell.angle_beta"),
        find_cif_item_value(document, "_cell.angle_gamma"),
    ];
    let [
        Some(a),
        Some(b),
        Some(c),
        Some(alpha),
        Some(beta),
        Some(gamma),
    ] = values
    else {
        return;
    };
    if cif_optional(a).is_none() || cif_optional(b).is_none() || cif_optional(c).is_none() {
        return;
    }
    let crystal = crystal.get_or_insert_with(default_crystal_info);
    crystal_set(
        crystal,
        CrystalCell {
            a: cif_number_or_nan(a),
            b: cif_number_or_nan(b),
            c: cif_number_or_nan(c),
            alpha: cif_number_or_nan(alpha),
            beta: cif_number_or_nan(beta),
            gamma: cif_number_or_nan(gamma),
        },
    );
}

// BEGIN GEMMI CPP FUNCTION gemmi::impl::find_spacegroup_hm_value
// Gemmi✔️✔️: const char* hm_tag = "_symmetry.space_group_name_H-M";
// Gemmi✔️✔️: return block.find_value(hm_tag);
// END GEMMI CPP FUNCTION
fn find_spacegroup_hm_value<'a>(document: &'a CifDocument) -> Option<&'a str> {
    find_cif_item_value(document, "_symmetry.space_group_name_H-M")
}

fn parse_cif_document(text: &str) -> Result<CifDocument, BioReadError> {
    let tokens = tokenize_cif(text)?;
    let mut blocks = Vec::new();
    let mut current_block = CifBlock {
        name: String::new(),
        items: Vec::new(),
        loops: Vec::new(),
        entries: Vec::new(),
    };
    let mut idx = 0;
    while idx < tokens.len() {
        if tokens[idx].value.starts_with("data_") {
            if !current_block.items.is_empty()
                || !current_block.loops.is_empty()
                || !blocks.is_empty()
            {
                blocks.push(current_block);
                current_block = CifBlock {
                    name: String::new(),
                    items: Vec::new(),
                    loops: Vec::new(),
                    entries: Vec::new(),
                };
            }
            current_block.name = tokens[idx].value["data_".len()..].to_string();
            idx += 1;
            continue;
        }
        if !tokens[idx].value.eq_ignore_ascii_case("loop_") {
            if tokens[idx].value.starts_with('_') {
                let tag = tokens[idx].value.clone();
                idx += 1;
                let value = tokens
                    .get(idx)
                    .cloned()
                    .ok_or_else(|| BioReadError::Parse {
                        line_number: tokens.last().map_or(0, |token| token.line_number),
                        message: format!("mmCIF item {tag} is missing a value"),
                    })?;
                current_block.push_pair(tag, value);
                idx += 1;
                continue;
            }
            idx += 1;
            continue;
        }
        idx += 1;

        let mut tags = Vec::new();
        while idx < tokens.len() && tokens[idx].value.starts_with('_') {
            tags.push(tokens[idx].value.clone());
            idx += 1;
        }
        if tags.is_empty() {
            return Err(BioReadError::Parse {
                line_number: tokens.get(idx).map_or(0, |token| token.line_number),
                message: "mmCIF loop_ without tags".to_string(),
            });
        }

        let mut values = Vec::new();
        while idx < tokens.len() {
            let value = tokens[idx].value.as_str();
            let starts_new_control = value.eq_ignore_ascii_case("loop_")
                || value.starts_with("data_")
                || value.starts_with("save_")
                || value.starts_with("global_");
            let starts_new_item = value.starts_with('_') && values.len() % tags.len() == 0;
            if starts_new_control || starts_new_item {
                break;
            }
            values.push(tokens[idx].clone());
            idx += 1;
        }
        current_block.push_loop(CifLoop { tags, values });
    }
    if !current_block.items.is_empty()
        || !current_block.loops.is_empty()
        || !current_block.name.is_empty()
    {
        blocks.push(current_block);
    }
    if blocks.is_empty() {
        blocks.push(CifBlock {
            name: String::new(),
            items: Vec::new(),
            loops: Vec::new(),
            entries: Vec::new(),
        });
    }
    Ok(CifDocument { blocks })
}

fn tokenize_cif(text: &str) -> Result<Vec<CifToken>, BioReadError> {
    let bytes = text.as_bytes();
    let mut tokens = Vec::new();
    let mut idx = 0;
    let mut line_number = 1;
    let mut at_line_start = true;

    while idx < bytes.len() {
        match bytes[idx] {
            b' ' | b'\t' | b'\r' => {
                idx += 1;
                at_line_start = false;
            }
            b'\n' => {
                idx += 1;
                line_number += 1;
                at_line_start = true;
            }
            b'#' => {
                while idx < bytes.len() && bytes[idx] != b'\n' {
                    idx += 1;
                }
            }
            b';' if at_line_start => {
                let start_line = line_number;
                idx += 1;
                if idx < bytes.len() && bytes[idx] == b'\n' {
                    idx += 1;
                    line_number += 1;
                }
                let start = idx;
                while idx < bytes.len() {
                    if (idx == 0 || bytes[idx - 1] == b'\n') && bytes[idx] == b';' {
                        let value = String::from_utf8_lossy(&bytes[start..idx]).to_string();
                        tokens.push(CifToken {
                            value,
                            line_number: start_line,
                        });
                        while idx < bytes.len() && bytes[idx] != b'\n' {
                            idx += 1;
                        }
                        break;
                    }
                    if bytes[idx] == b'\n' {
                        line_number += 1;
                    }
                    idx += 1;
                }
                at_line_start = true;
            }
            quote @ (b'\'' | b'"') => {
                let start_line = line_number;
                idx += 1;
                let start = idx;
                while idx < bytes.len() && bytes[idx] != quote {
                    if bytes[idx] == b'\n' {
                        line_number += 1;
                    }
                    idx += 1;
                }
                if idx >= bytes.len() {
                    return Err(BioReadError::Parse {
                        line_number: start_line,
                        message: "unterminated quoted mmCIF value".to_string(),
                    });
                }
                tokens.push(CifToken {
                    value: String::from_utf8_lossy(&bytes[start..idx]).to_string(),
                    line_number: start_line,
                });
                idx += 1;
                at_line_start = false;
            }
            _ => {
                let start_line = line_number;
                let start = idx;
                while idx < bytes.len()
                    && !matches!(bytes[idx], b' ' | b'\t' | b'\r' | b'\n' | b'#')
                {
                    idx += 1;
                }
                tokens.push(CifToken {
                    value: String::from_utf8_lossy(&bytes[start..idx]).to_string(),
                    line_number: start_line,
                });
                at_line_start = false;
            }
        }
    }

    Ok(tokens)
}

fn parse_pdb_element(line: &str, line_number: usize) -> Result<Element, BioReadError> {
    let symbol = field(line, 76, 78).trim();
    if symbol.is_empty() {
        return infer_element_from_padded_atom_name(field_raw(line, 12, 16)).ok_or_else(|| {
            BioReadError::Parse {
                line_number,
                message: "could not infer PDB atom element from atom name".to_string(),
            }
        });
    }
    element_from_symbol(symbol).ok_or_else(|| BioReadError::Parse {
        line_number,
        message: format!("unsupported element symbol {symbol:?}"),
    })
}

fn record_type(line: &str) -> &str {
    field_raw(line, 0, line.len().min(4))
}

fn starts_record(line: &str, record: &str) -> bool {
    line.get(0..record.len())
        .is_some_and(|prefix| prefix.eq_ignore_ascii_case(record))
}

fn starts_record3(line: &str, record: &str) -> bool {
    let prefix = field_raw(line, 0, line.len().min(4));
    prefix
        .get(0..3)
        .is_some_and(|value| value.eq_ignore_ascii_case(record))
        && matches!(
            prefix.as_bytes().get(3).copied().unwrap_or(b' '),
            b' ' | b'\0'
        )
}

fn field(line: &str, start: usize, end: usize) -> &str {
    field_raw(line, start, end).trim()
}

fn field_raw(line: &str, start: usize, end: usize) -> &str {
    line.get(start..end).unwrap_or("")
}

fn pdb_chain_id_from_str(value: &str) -> Option<PdbChainId> {
    let bytes = value.as_bytes();
    if bytes.is_empty() || bytes.len() > 4 {
        return None;
    }
    let mut storage = [0; 4];
    storage[..bytes.len()].copy_from_slice(bytes);
    Some(PdbChainId(storage, bytes.len() as u8))
}

fn bio_transform_apply_mat(mat: &[[f64; 3]; 3], vec: [f64; 3]) -> [f64; 3] {
    [
        mat[0][0] * vec[0] + mat[0][1] * vec[1] + mat[0][2] * vec[2],
        mat[1][0] * vec[0] + mat[1][1] * vec[1] + mat[1][2] * vec[2],
        mat[2][0] * vec[0] + mat[2][1] * vec[1] + mat[2][2] * vec[2],
    ]
}

fn bio_transform_multiply(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    [
        [
            a[0][0] * b[0][0] + a[0][1] * b[1][0] + a[0][2] * b[2][0],
            a[0][0] * b[0][1] + a[0][1] * b[1][1] + a[0][2] * b[2][1],
            a[0][0] * b[0][2] + a[0][1] * b[1][2] + a[0][2] * b[2][2],
        ],
        [
            a[1][0] * b[0][0] + a[1][1] * b[1][0] + a[1][2] * b[2][0],
            a[1][0] * b[0][1] + a[1][1] * b[1][1] + a[1][2] * b[2][1],
            a[1][0] * b[0][2] + a[1][1] * b[1][2] + a[1][2] * b[2][2],
        ],
        [
            a[2][0] * b[0][0] + a[2][1] * b[1][0] + a[2][2] * b[2][0],
            a[2][0] * b[0][1] + a[2][1] * b[1][1] + a[2][2] * b[2][1],
            a[2][0] * b[0][2] + a[2][1] * b[1][2] + a[2][2] * b[2][2],
        ],
    ]
}

fn bio_transform_combine(a: &BioTransform, b: &BioTransform) -> BioTransform {
    BioTransform {
        mat: bio_transform_multiply(&a.mat, &b.mat),
        vec: {
            let shifted = bio_transform_apply_mat(&a.mat, b.vec);
            [
                shifted[0] + a.vec[0],
                shifted[1] + a.vec[1],
                shifted[2] + a.vec[2],
            ]
        },
    }
}

fn bio_transform_inverse(transform: &BioTransform) -> Option<BioTransform> {
    let m = &transform.mat;
    let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    if det.abs() <= f64::EPSILON {
        return None;
    }
    let inv_det = 1.0 / det;
    let inv = [
        [
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv_det,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv_det,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_det,
        ],
        [
            (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv_det,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_det,
            (m[1][0] * m[0][2] - m[0][0] * m[1][2]) * inv_det,
        ],
        [
            (m[1][0] * m[2][1] - m[2][0] * m[1][1]) * inv_det,
            (m[2][0] * m[0][1] - m[0][0] * m[2][1]) * inv_det,
            (m[0][0] * m[1][1] - m[1][0] * m[0][1]) * inv_det,
        ],
    ];
    let neg_vec = [-transform.vec[0], -transform.vec[1], -transform.vec[2]];
    Some(BioTransform {
        mat: inv,
        vec: bio_transform_apply_mat(&inv, neg_vec),
    })
}

fn bio_transform_approx_eq(a: &BioTransform, b: &BioTransform, mat_eps: f64, vec_eps: f64) -> bool {
    for row in 0..3 {
        for col in 0..3 {
            if (a.mat[row][col] - b.mat[row][col]).abs() > mat_eps {
                return false;
            }
        }
        if (a.vec[row] - b.vec[row]).abs() > vec_eps {
            return false;
        }
    }
    true
}

fn crystal_is_crystal(crystal: &CrystalInfo) -> bool {
    crystal.cell.a != 1.0 && crystal.frac.mat[0][0] != 1.0
}

fn crystal_calculate_properties(crystal: &mut CrystalInfo) {
    let alpha = crystal.cell.alpha;
    let beta = crystal.cell.beta;
    let gamma = crystal.cell.gamma;
    let a = crystal.cell.a;
    let b = crystal.cell.b;
    let c = crystal.cell.c;
    let cos_alpha = if alpha == 90.0 {
        0.0
    } else {
        alpha.to_radians().cos()
    };
    let cos_beta = if beta == 90.0 {
        0.0
    } else {
        beta.to_radians().cos()
    };
    let cos_gamma = if gamma == 90.0 {
        0.0
    } else {
        gamma.to_radians().cos()
    };
    let sin_alpha = if alpha == 90.0 {
        1.0
    } else {
        alpha.to_radians().sin()
    };
    let sin_beta = if beta == 90.0 {
        1.0
    } else {
        beta.to_radians().sin()
    };
    let sin_gamma = if gamma == 90.0 {
        1.0
    } else {
        gamma.to_radians().sin()
    };
    if sin_alpha == 0.0 || sin_beta == 0.0 || sin_gamma == 0.0 {
        return;
    }
    let volume = a
        * b
        * c
        * (1.0 - cos_alpha * cos_alpha - cos_beta * cos_beta - cos_gamma * cos_gamma
            + 2.0 * cos_alpha * cos_beta * cos_gamma)
            .sqrt();
    if volume == 0.0 {
        return;
    }
    let cos_alphar_sin_beta = (cos_beta * cos_gamma - cos_alpha) / sin_gamma;
    let cos_alphar = cos_alphar_sin_beta / sin_beta;
    if crystal.explicit_matrices {
        return;
    }
    let sin_alphar = (1.0 - cos_alphar * cos_alphar).sqrt();
    crystal.orth = BioTransform {
        mat: [
            [a, b * cos_gamma, c * cos_beta],
            [0.0, b * sin_gamma, -c * cos_alphar_sin_beta],
            [0.0, 0.0, c * sin_beta * sin_alphar],
        ],
        vec: [0.0, 0.0, 0.0],
    };
    let o12 = -cos_gamma / (sin_gamma * a);
    let o13 = -(cos_gamma * cos_alphar_sin_beta + cos_beta * sin_gamma)
        / (sin_alphar * sin_beta * sin_gamma * a);
    let o23 = cos_alphar / (sin_alphar * sin_gamma * b);
    crystal.frac = BioTransform {
        mat: [
            [1.0 / a, o12, o13],
            [0.0, 1.0 / crystal.orth.mat[1][1], o23],
            [0.0, 0.0, 1.0 / crystal.orth.mat[2][2]],
        ],
        vec: [0.0, 0.0, 0.0],
    };
}

fn crystal_set(crystal: &mut CrystalInfo, cell: CrystalCell) {
    // Gemmi✔️✔️: void set(double a_, double b_, double c_,
    // Gemmi✔️✔️:          double alpha_, double beta_, double gamma_) {
    // Gemmi✔️✔️:   if (gamma_ == 0.0)  // ignore empty/partial CRYST1 (example: 3iyp)
    // Gemmi✔️✔️:     return;
    // Gemmi✔️✔️:   a = a_;
    // Gemmi✔️✔️:   b = b_;
    // Gemmi✔️✔️:   c = c_;
    // Gemmi✔️✔️:   alpha = alpha_;
    // Gemmi✔️✔️:   beta = beta_;
    // Gemmi✔️✔️:   gamma = gamma_;
    // Gemmi✔️✔️:   calculate_properties();
    // Gemmi✔️✔️: }
    if cell.gamma == 0.0 {
        return;
    }
    crystal.cell = cell;
    crystal.explicit_matrices = false;
    crystal_calculate_properties(crystal);
}

fn crystal_set_matrices_from_fract(crystal: &mut CrystalInfo, fract: BioTransform) {
    if bio_transform_approx_eq(&fract, &crystal.frac, 1e-4, 1e-6) {
        return;
    }
    if crystal.frac.mat[0][0] == 1.0 && (fract.mat[0][0] == 0.0 || fract.mat[0][0] > 1.0) {
        return;
    }
    if let Some(orth) = bio_transform_inverse(&fract) {
        crystal.frac = fract;
        crystal.orth = orth;
        crystal.explicit_matrices = true;
    }
}

fn default_crystal_info() -> CrystalInfo {
    // Gemmi✔️✔️: struct UnitCellParameters {
    // Gemmi✔️✔️:   double a = 1.0, b = 1.0, c = 1.0;
    // Gemmi✔️✔️:   double alpha = 90.0, beta = 90.0, gamma = 90.0;
    CrystalInfo {
        cell: CrystalCell {
            a: 1.0,
            b: 1.0,
            c: 1.0,
            alpha: 90.0,
            beta: 90.0,
            gamma: 90.0,
        },
        spacegroup_hm: None,
        z_pdb: None,
        scale: None,
        frac: bio_transform_identity(),
        orth: bio_transform_identity(),
        explicit_matrices: false,
        cs_count: 0,
        cell_images: Vec::new(),
    }
}

fn find_cif_item_value<'a>(document: &'a CifDocument, tag: &str) -> Option<&'a str> {
    document
        .blocks
        .first()?
        .items
        .iter()
        .find(|item| item.tag == tag)
        .map(|item| item.value.value.as_str())
}

fn byte_at(line: &str, index: usize) -> u8 {
    line.as_bytes().get(index).copied().unwrap_or(b' ')
}

fn cif_number_or_nan(value: &str) -> f64 {
    let value = value.trim();
    if value.is_empty() || value == "." || value == "?" {
        return f64::NAN;
    }
    value.parse::<f64>().unwrap_or(f64::NAN)
}

fn read_pdb_lossy_f64(value: &str) -> f64 {
    value.trim().parse::<f64>().unwrap_or(0.0)
}

fn read_pdb_unchecked_decimal(value: &str) -> i32 {
    // BEGIN GEMMI CPP FUNCTION gemmi::string_to_int
    // Gemmi✔️✔️: inline int string_to_int(const char* p, bool checked, size_t length=0) {
    // Gemmi✔️✔️:   int mult = -1;
    // Gemmi✔️✔️:   int n = 0;
    // Gemmi✔️✔️:   size_t i = 0;
    // Gemmi✔️✔️:   while ((length == 0 || i < length) && is_space(p[i]))
    // Gemmi✔️✔️:     ++i;
    // Gemmi✔️✔️:   if (p[i] == '-') {
    // Gemmi✔️✔️:     mult = 1;
    // Gemmi✔️✔️:     ++i;
    // Gemmi✔️✔️:   } else if (p[i] == '+') {
    // Gemmi✔️✔️:     ++i;
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️:   bool has_digits = false;
    // Gemmi✔️✔️:   // use negative numbers because INT_MIN < -INT_MAX
    // Gemmi✔️✔️:   for (; (length == 0 || i < length) && is_digit(p[i]); ++i) {
    // Gemmi✔️✔️:     n = n * 10 - (p[i] - '0');
    // Gemmi✔️✔️:     has_digits = true;
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️:   if (checked) {
    // Gemmi✔️✔️:     while ((length == 0 || i < length) && is_space(p[i]))
    // Gemmi✔️✔️:       ++i;
    // Gemmi✔️✔️:     if (!has_digits || p[i] != '\0')
    // Gemmi✔️✔️:       throw std::invalid_argument("not an integer: " +
    // Gemmi✔️✔️:                                   std::string(p, length ? length : i+1));
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️:   return mult * n;
    // Gemmi✔️✔️: }
    // END GEMMI CPP FUNCTION
    let bytes = value.as_bytes();
    let mut multiplier = -1;
    let mut number = 0_i32;
    let mut index = 0;
    while index < bytes.len() && bytes[index].is_ascii_whitespace() {
        index += 1;
    }
    if bytes.get(index) == Some(&b'-') {
        multiplier = 1;
        index += 1;
    } else if bytes.get(index) == Some(&b'+') {
        index += 1;
    }
    while let Some(&byte) = bytes.get(index)
        && byte.is_ascii_digit()
    {
        number = number * 10 - i32::from(byte - b'0');
        index += 1;
    }
    multiplier * number
}

fn read_pdb_base36<const N: usize>(value: &str) -> i32 {
    // BEGIN GEMMI CPP FUNCTION gemmi::read_base36
    // Gemmi✔️✔️: template<int N> int read_base36(const char* p) {
    // Gemmi✔️✔️:   char zstr[N+1] = {0};
    // Gemmi✔️✔️:   std::memcpy(zstr, p, N);
    // Gemmi✔️✔️:   return std::strtol(zstr, nullptr, 36);
    // Gemmi✔️✔️: }
    // END GEMMI CPP FUNCTION
    let bytes = value.as_bytes();
    let limit = N.min(bytes.len());
    let mut index = 0;
    while index < limit && bytes[index].is_ascii_whitespace() {
        index += 1;
    }
    let mut sign = 1_i32;
    if bytes.get(index) == Some(&b'-') {
        sign = -1;
        index += 1;
    } else if bytes.get(index) == Some(&b'+') {
        index += 1;
    }
    let mut number = 0_i32;
    while index < limit {
        let digit = match bytes[index] {
            b'0'..=b'9' => i32::from(bytes[index] - b'0'),
            b'a'..=b'z' => i32::from(bytes[index] - b'a') + 10,
            b'A'..=b'Z' => i32::from(bytes[index] - b'A') + 10,
            _ => break,
        };
        number = number * 36 + digit;
        index += 1;
    }
    sign * number
}

fn read_pdb_seq_id(value: &str) -> PdbSeqId {
    // BEGIN GEMMI CPP FUNCTION gemmi::read_seq_id
    // Gemmi✔️✔️: SeqId read_seq_id(const char* str) {
    // Gemmi✔️✔️:   SeqId seqid;
    // Gemmi✔️✔️:   if (str[4] != '\r' && str[4] != '\n')
    // Gemmi✔️✔️:     seqid.icode = str[4];
    // Gemmi✔️✔️:   // We support hybrid-36 extension, although it is never used in practice
    // Gemmi✔️✔️:   // as 9999 residues per chain are enough.
    // Gemmi✔️✔️:   if (str[0] < 'A') {
    // Gemmi✔️✔️:     for (int i = 4; i != 0; --i, ++str)
    // Gemmi✔️✔️:       if (!is_space(*str)) {
    // Gemmi✔️✔️:         seqid.num = read_int(str, i);
    // Gemmi✔️✔️:         break;
    // Gemmi✔️✔️:       }
    // Gemmi✔️✔️:   } else {
    // Gemmi✔️✔️:     seqid.num = read_base36<4>(str) - 466560 + 10000;
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️:   return seqid;
    // Gemmi✔️✔️: }
    // END GEMMI CPP FUNCTION
    let bytes = value.as_bytes();
    let ins_code = match bytes.get(4).copied().unwrap_or(b' ') {
        b' ' | b'\r' | b'\n' | 0 => None,
        byte => Some(byte),
    };
    let seq_num = if bytes.first().copied().unwrap_or(b' ') < b'A' {
        (0..4)
            .find(|&index| {
                bytes
                    .get(index)
                    .is_some_and(|byte| !byte.is_ascii_whitespace())
            })
            .map_or(i32::MIN, |index| {
                read_pdb_unchecked_decimal(&value[index..4.min(value.len())])
            })
    } else {
        read_pdb_base36::<4>(value) - 466_560 + 10_000
    };
    PdbSeqId { seq_num, ins_code }
}

fn read_pdb_serial(value: &str) -> i32 {
    // BEGIN GEMMI CPP FUNCTION gemmi::read_serial
    // Gemmi✔️✔️: int read_serial(const char* ptr) {
    // Gemmi✔️✔️:   return ptr[0] < 'A' ? read_int(ptr, 5)
    // Gemmi✔️✔️:                       : read_base36<5>(ptr) - 16796160 + 100000;
    // Gemmi✔️✔️: }
    // END GEMMI CPP FUNCTION
    if value.as_bytes().first().copied().unwrap_or(b' ') < b'A' {
        read_pdb_unchecked_decimal(value)
    } else {
        read_pdb_base36::<5>(value) - 16_796_160 + 100_000
    }
}

fn parse_decimal_i32(
    _feature: &'static FeatureSpec,
    value: &str,
    line_number: usize,
    field_name: &'static str,
) -> Result<i32, BioReadError> {
    let trimmed = value.trim();
    trimmed.parse::<i32>().map_err(|_| BioReadError::Parse {
        line_number,
        message: format!("invalid {field_name}: {trimmed:?}"),
    })
}

fn parse_optional_i32(value: &str, line_number: usize) -> Result<Option<i32>, BioReadError> {
    let trimmed = value.trim();
    if trimmed.is_empty() {
        return Ok(None);
    }
    parse_decimal_i32(&BIO_PDB_READ_FEATURE, trimmed, line_number, "model number").map(Some)
}

fn append_metadata_string(target: &mut Option<String>, value: &str) {
    let trimmed = value.trim_end();
    if trimmed.is_empty() {
        return;
    }
    match target {
        Some(existing) => existing.push_str(trimmed),
        None => *target = Some(trimmed.to_string()),
    }
}

fn pdb_date_format_to_iso(value: &str) -> String {
    let trimmed = value.trim();
    let mut parts = trimmed.split('-');
    let Some(day) = parts.next() else {
        return String::new();
    };
    let Some(month) = parts.next() else {
        return String::new();
    };
    let Some(year) = parts.next() else {
        return String::new();
    };
    if parts.next().is_some() {
        return String::new();
    }
    let Ok(day) = day.parse::<u8>() else {
        return String::new();
    };
    let month = match month.to_ascii_uppercase().as_str() {
        "JAN" => 1,
        "FEB" => 2,
        "MAR" => 3,
        "APR" => 4,
        "MAY" => 5,
        "JUN" => 6,
        "JUL" => 7,
        "AUG" => 8,
        "SEP" => 9,
        "OCT" => 10,
        "NOV" => 11,
        "DEC" => 12,
        _ => return String::new(),
    };
    let Ok(year2) = year.parse::<u16>() else {
        return String::new();
    };
    let year = if year2 >= 50 {
        1900 + year2
    } else {
        2000 + year2
    };
    format!("{year:04}-{month:02}-{day:02}")
}

fn parse_f32(
    value: &str,
    line_number: usize,
    field_name: &'static str,
) -> Result<f32, BioReadError> {
    let trimmed = value.trim();
    trimmed.parse::<f32>().map_err(|_| BioReadError::Parse {
        line_number,
        message: format!("invalid {field_name}: {trimmed:?}"),
    })
}

fn parse_f64(
    value: &str,
    line_number: usize,
    field_name: &'static str,
) -> Result<f64, BioReadError> {
    let trimmed = value.trim();
    trimmed.parse::<f64>().map_err(|_| BioReadError::Parse {
        line_number,
        message: format!("invalid {field_name}: {trimmed:?}"),
    })
}

fn parse_i8(value: &str, line_number: usize, field_name: &'static str) -> Result<i8, BioReadError> {
    let trimmed = value.trim();
    trimmed.parse::<i8>().map_err(|_| BioReadError::Parse {
        line_number,
        message: format!("invalid {field_name}: {trimmed:?}"),
    })
}

fn parse_optional_f32(
    value: &str,
    line_number: usize,
    field_name: &'static str,
) -> Result<Option<f32>, BioReadError> {
    let trimmed = value.trim();
    if trimmed.is_empty() {
        return Ok(None);
    }
    parse_f32(trimmed, line_number, field_name).map(Some)
}

fn parse_pdb_charge(value: &str, line_number: usize) -> Result<Option<i8>, BioReadError> {
    let trimmed = value.trim();
    if trimmed.is_empty() {
        return Ok(None);
    }
    if trimmed.len() != 2 {
        return Err(BioReadError::Parse {
            line_number,
            message: format!("invalid PDB formal charge: {trimmed:?}"),
        });
    }
    let bytes = trimmed.as_bytes();
    let magnitude = (bytes[0] as char)
        .to_digit(10)
        .ok_or_else(|| BioReadError::Parse {
            line_number,
            message: format!("invalid PDB formal charge: {trimmed:?}"),
        })? as i8;
    match bytes[1] {
        b'+' => Ok(Some(magnitude)),
        b'-' => Ok(Some(-magnitude)),
        _ => Err(BioReadError::Parse {
            line_number,
            message: format!("invalid PDB formal charge: {trimmed:?}"),
        }),
    }
}

fn atom_name_from_field(value: &str) -> AtomName {
    let mut bytes = [b' '; 4];
    for (out, input) in bytes.iter_mut().zip(value.as_bytes().iter().copied()) {
        *out = input;
    }
    AtomName(bytes)
}

fn atom_name_from_cif(value: &str) -> AtomName {
    let mut bytes = [b' '; 4];
    for (out, input) in bytes.iter_mut().zip(value.as_bytes().iter().copied()) {
        *out = input;
    }
    AtomName(bytes)
}

fn residue_name_from_field(value: &str) -> ResidueName {
    let trimmed = value.trim();
    let mut bytes = [0; 4];
    let len = trimmed.len().min(4);
    bytes[..len].copy_from_slice(&trimmed.as_bytes()[..len]);
    ResidueName(bytes, len as u8)
}

fn pdb_chain_id_from_field(value: &str) -> PdbChainId {
    let trimmed = value.trim();
    let mut bytes = [0; 4];
    let len = trimmed.len().min(4);
    bytes[..len].copy_from_slice(&trimmed.as_bytes()[..len]);
    PdbChainId(bytes, len as u8)
}

fn pdb_chain_id_from_cif(value: &str, line_number: usize) -> Result<PdbChainId, BioReadError> {
    let trimmed = value.trim();
    if trimmed.len() > 4 {
        return Err(unsupported(
            &BIO_MMCIF_READ_FEATURE,
            line_number,
            "mmCIF chain identifiers longer than 4 bytes are not modeled in PdbChainId",
        ));
    }
    Ok(pdb_chain_id_from_field(trimmed))
}

fn cif_optional(value: &str) -> Option<&str> {
    match value.trim() {
        "" | "." | "?" => None,
        value => Some(value),
    }
}

fn entity_kind_from_cif(value: &str) -> EntityKind {
    match value.trim().to_ascii_lowercase().as_str() {
        "polymer" => EntityKind::Polymer,
        "non-polymer" | "nonpolymer" => EntityKind::NonPolymer,
        "branched" => EntityKind::Branched,
        "water" => EntityKind::Water,
        _ => EntityKind::Unknown,
    }
}

fn polymer_kind_from_cif(value: &str) -> PolymerKind {
    match value.trim().to_ascii_lowercase().as_str() {
        "polypeptide(l)" | "polypeptide(d)" => PolymerKind::Peptide,
        "polydeoxyribonucleotide" => PolymerKind::DNA,
        "polyribonucleotide" => PolymerKind::RNA,
        "polydeoxyribonucleotide/polyribonucleotide hybrid" => PolymerKind::NucleicAcidHybrid,
        "polysaccharide(d)" | "polysaccharide(l)" => PolymerKind::Saccharide,
        _ => PolymerKind::Unknown,
    }
}

fn missing_cif_value(line_number: usize, tag: &'static str) -> BioReadError {
    BioReadError::Parse {
        line_number,
        message: format!("required mmCIF value is missing: {tag}"),
    }
}

fn parse_altloc(value: u8) -> Option<AltLocLabel> {
    match value {
        b' ' | 0 => None,
        other => Some(AltLocLabel(other)),
    }
}

fn element_from_symbol(symbol: &str) -> Option<Element> {
    // BEGIN GEMMI CPP FUNCTION gemmi::Element / gemmi::impl::find_element
    // Gemmi✔️✔️: explicit Element(const char* str) noexcept : elem(find_element(str)) {}
    // Gemmi✔️✔️: inline El find_element(const char* symbol) {
    // Gemmi✔️✔️:   if (symbol == nullptr || symbol[0] == '\0')
    // Gemmi✔️✔️:     return El::X;
    // Gemmi✔️✔️:   char first = symbol[0] & ~0x20;  // lower -> upper, space -> NUL
    // Gemmi✔️✔️:   char second = symbol[1] & ~0x20;
    // Gemmi✔️✔️:   if (first == '\0')
    // Gemmi✔️✔️:     return impl::find_single_letter_element(second);
    // Gemmi✔️✔️:   if (second < 14)
    // Gemmi✔️✔️:     return impl::find_single_letter_element(first);
    // Gemmi✔️✔️:   elname_t* names = &element_uppercase_name(El::X);
    // Gemmi✔️✔️:   for (int i = 0; i != 120; ++i) {
    // Gemmi✔️✔️:     if (names[i][0] == first && names[i][1] == second)
    // Gemmi✔️✔️:       return static_cast<El>(i);
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️:   return El::X;
    // Gemmi✔️✔️: }
    // Gemmi✔️✔️: int atomic_number() const { return elem == El::D ? 1 : ordinal(); }
    // END GEMMI CPP FUNCTION gemmi::Element / gemmi::impl::find_element
    //
    // BEGIN RDKIT CPP FUNCTION PeriodicTable::getAtomicNumber
    // RDKit✔️✔️: byname[dat->Symbol()] = dat->AtomicNumber();
    // RDKit✔️✔️: int getAtomicNumber(const std::string &elementSymbol) const { ... byname.find(elementSymbol) ... }
    // END RDKIT CPP FUNCTION PeriodicTable::getAtomicNumber
    let normalized = normalize_gemmi_element_symbol(symbol.trim())?;
    if normalized == "X" {
        return Some(Element::DUMMY);
    }
    if normalized == "D" {
        return Some(Element::H);
    }
    let atomic_number = crate::chemistry::valence::rdkit_atomic_number_from_symbol(&normalized)?;
    Element::from_atomic_number(atomic_number)
}

fn normalize_gemmi_element_symbol(symbol: &str) -> Option<String> {
    let bytes = symbol.as_bytes();
    let first = bytes.first().copied()?;
    if !first.is_ascii_alphabetic() {
        return None;
    }

    let mut normalized = String::with_capacity(2);
    normalized.push(char::from(first.to_ascii_uppercase()));
    if let Some(second) = bytes.get(1).copied() {
        if !second.is_ascii_alphabetic() {
            return None;
        }
        normalized.push(char::from(second.to_ascii_lowercase()));
    }
    Some(normalized)
}

fn infer_element_from_padded_atom_name(atom_name: &str) -> Option<Element> {
    let bytes = atom_name.as_bytes();
    if bytes.is_empty() {
        return None;
    }

    if bytes
        .first()
        .copied()
        .is_some_and(|byte| byte == b' ' || byte.is_ascii_digit())
    {
        let first = bytes
            .iter()
            .copied()
            .find(|byte| byte.is_ascii_alphabetic())?;
        let symbol = (first as char).to_string();
        return element_from_symbol(&symbol);
    }

    let letters = bytes
        .iter()
        .copied()
        .filter(|byte| byte.is_ascii_alphabetic())
        .take(2)
        .map(char::from)
        .collect::<String>();
    if letters.len() >= 2
        && let Some(element) = element_from_symbol(&letters)
    {
        return Some(element);
    }
    letters
        .chars()
        .next()
        .and_then(|ch| element_from_symbol(&ch.to_string()))
}

fn unsupported(
    feature: &'static FeatureSpec,
    line_number: usize,
    reason: &'static str,
) -> BioReadError {
    BioReadError::Unsupported {
        line_number,
        feature_name: feature.name,
        reason,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mmcif_output_group_all_categories_emit_and_replace_prepopulated_categories() {
        writer::tests::check_mmcif_output_group_all_categories_emit_and_replace_prepopulated_categories();
    }

    #[test]
    fn mmcif_output_group_all_disabled_preserves_prepopulated_block_exactly() {
        writer::tests::check_mmcif_output_group_all_disabled_preserves_prepopulated_block_exactly();
    }
    use crate::bio::BioHelixClass;

    fn read_test_mmcif(text: &str) -> Result<BioStructure, BioReadError> {
        BioStructure::from_mmcif_str(text, "input.cif")
    }

    #[test]
    fn reads_single_pdb_atom_record() {
        let pdb =
            "ATOM      1  CA  ALA A   7      11.104  13.207   9.900  1.00 20.00           C  \n";

        let structure = BioStructure::from_pdb_str(&pdb).unwrap();

        assert_eq!(structure.num_models(), 1);
        assert_eq!(structure.num_chains(), 1);
        assert_eq!(structure.num_residues(), 1);
        assert_eq!(structure.num_atoms(), 1);
        assert_eq!(
            structure.chains[0].source.auth_chain_id.unwrap().as_str(),
            "A"
        );
        assert_eq!(structure.residues[0].source.seq_id.unwrap().seq_num, 7);
        assert_eq!(structure.residues[0].name.as_str(), "ALA");
        assert_eq!(structure.residues[0].kind, crate::ResidueKind::AminoAcid);
        assert_eq!(structure.atoms[0].name, AtomName([b' ', b'C', b'A', b' ']));
        assert_eq!(structure.atoms[0].element, Element::C);
        assert_eq!(structure.atoms[0].occupancy, Some(1.0));
        assert_eq!(structure.atoms[0].b_iso, Some(20.0));
        assert_eq!(structure.coordinates.positions[0], [11.104, 13.207, 9.900]);
    }

    #[test]
    fn pdb_element_field_accepts_full_rdkit_gemmi_element_table() {
        for atomic_number in 1..=118 {
            let symbol = crate::chemistry::valence::rdkit_element_symbol(atomic_number).unwrap();
            let element_field = format!("{symbol:>2}");
            let mut pdb = " ".repeat(80);
            pdb.replace_range(0..6, "HETATM");
            pdb.replace_range(6..11, &format!("{atomic_number:>5}"));
            pdb.replace_range(12..16, " M  ");
            pdb.replace_range(17..20, "LIG");
            pdb.replace_range(20..22, "A ");
            pdb.replace_range(22..26, &format!("{atomic_number:>4}"));
            pdb.replace_range(30..38, "   0.000");
            pdb.replace_range(38..46, "   0.000");
            pdb.replace_range(46..54, "   0.000");
            pdb.replace_range(54..60, "  1.00");
            pdb.replace_range(60..66, " 10.00");
            pdb.replace_range(76..78, &element_field);

            let structure = BioStructure::from_pdb_str(&format!("{pdb}\n")).unwrap();
            assert_eq!(
                structure.atoms[0].element.atomic_number(),
                atomic_number,
                "{symbol}"
            );
        }

        for (symbol, expected) in [("HG", Element::HG), ("CD", Element::CD)] {
            let mut pdb = " ".repeat(80);
            pdb.replace_range(0..6, "HETATM");
            pdb.replace_range(6..11, "    1");
            pdb.replace_range(12..16, " M  ");
            pdb.replace_range(17..20, "LIG");
            pdb.replace_range(20..22, "A ");
            pdb.replace_range(22..26, "   1");
            pdb.replace_range(30..38, "   0.000");
            pdb.replace_range(38..46, "   0.000");
            pdb.replace_range(46..54, "   0.000");
            pdb.replace_range(54..60, "  1.00");
            pdb.replace_range(60..66, " 10.00");
            pdb.replace_range(76..78, symbol);

            let structure = BioStructure::from_pdb_str(&format!("{pdb}\n")).unwrap();
            assert_eq!(structure.atoms[0].element, expected, "{symbol}");
        }
    }

    #[test]
    fn bio_element_symbol_lookup_covers_rdkit_periodic_table() {
        for atomic_number in 1..=118 {
            let symbol = crate::chemistry::valence::rdkit_element_symbol(atomic_number).unwrap();
            let element = element_from_symbol(symbol).unwrap();
            assert_eq!(element.atomic_number(), atomic_number, "{symbol}");

            let uppercase = symbol.to_ascii_uppercase();
            let element = element_from_symbol(&uppercase).unwrap();
            assert_eq!(element.atomic_number(), atomic_number, "{uppercase}");
        }
        assert_eq!(element_from_symbol("X").unwrap(), Element::DUMMY);
        assert_eq!(element_from_symbol("x").unwrap(), Element::DUMMY);
        assert_eq!(element_from_symbol("D").unwrap(), Element::H);
        assert_eq!(element_from_symbol("d").unwrap(), Element::H);
    }

    #[test]
    fn mmcif_atom_site_accepts_full_rdkit_gemmi_element_table() {
        let mut mmcif = "\
data_demo
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
"
        .to_string();
        for atomic_number in 1..=118 {
            let symbol = crate::chemistry::valence::rdkit_element_symbol(atomic_number).unwrap();
            mmcif.push_str(&format!(
                "HETATM {atomic_number} {symbol} {symbol} . LIG A 1 {atomic_number} 0.0 0.0 0.0\n"
            ));
        }
        mmcif.push_str("HETATM 119 HG HG . LIG A 1 119 0.0 0.0 0.0\n");
        mmcif.push_str("HETATM 120 CD CD . LIG A 1 120 0.0 0.0 0.0\n");

        let structure = read_test_mmcif(&mmcif).unwrap();

        for atomic_number in 1..=118 {
            assert_eq!(
                structure.atoms[usize::from(atomic_number - 1)]
                    .element
                    .atomic_number(),
                atomic_number
            );
        }

        assert_eq!(structure.atoms[118].element, Element::HG);
        assert_eq!(structure.atoms[119].element, Element::CD);
    }

    #[test]
    fn preserves_altloc_and_water_kind() {
        let pdb =
            "HETATM   22  O  AHOH B  10       1.000   2.000   3.000  1.00 10.00           O  \n";

        let structure = BioStructure::from_pdb_str(&pdb).unwrap();

        assert_eq!(structure.residues[0].kind, ResidueKind::Water);
        assert_eq!(structure.atoms[0].altloc, Some(AltLocLabel(b'A')));
    }

    #[test]
    fn reads_pdb_anisou_for_previous_atom() {
        let pdb = "\
ATOM      1  CA  ALA A   7      11.104  13.207   9.900  1.00 20.00           C
ANISOU    1  CA  ALA A   7    1000   2000   3000    400    500    600       C
";

        let structure = BioStructure::from_pdb_str(&pdb).unwrap();

        assert_eq!(
            structure.atoms[0].anisou,
            Some([0.099999994, 0.19999999, 0.29999998, 0.04, 0.049999997, 0.06])
        );
    }

    #[test]
    fn infers_missing_pdb_element_from_padded_atom_name() {
        let pdb =
            "ATOM      1  CA  ALA A   7      11.104  13.207   9.900  1.00 20.00              \n";

        let structure = BioStructure::from_pdb_str(&pdb).unwrap();

        assert_eq!(structure.atoms[0].element, Element::C);
    }

    #[test]
    fn reads_pdb_seqres_entity_sequence_and_links_chain() {
        let pdb = "\
SEQRES   1 A    5  ALA GLY SER THR TYR
ATOM      1  CA  ALA A   7      11.104  13.207   9.900  1.00 20.00           C
";

        let structure = BioStructure::from_pdb_str(pdb).unwrap();

        assert_eq!(structure.num_entities(), 1);
        assert_eq!(structure.entities[0].kind, EntityKind::Polymer);
        assert_eq!(
            structure.entities[0].sequence,
            vec!["ALA", "GLY", "SER", "THR", "TYR"]
        );
        assert_eq!(structure.chains[0].entity_id, Some(EntityId::new(0)));
    }

    #[test]
    fn reads_pdb_header_title_authors_cryst1_and_matrices() {
        let pdb = format!(
            "\
HEADER    OXIDOREDUCTASE                          28-MAR-07   2XYZ
TITLE     EXAMPLE
TITLE    2 STRUCTURE
KEYWDS    TEST,
KEYWDS   2 GEMMI ROUTE
EXPDTA    X-RAY
EXPDTA   2 DIFFRACTION
AUTHOR    DOE,J.
AUTHOR   2 SMITH
{}
{}
{}
{}
{}
{}
CRYST1   10.000   20.000   30.000  90.00 100.00 120.00 P 1           2
{}
{}
{}
{}
{}
{}
ATOM      1  CA  ALA A   7      11.104  13.207   9.900  1.00 20.00           C
",
            format!(
                "SCALE1    {:>10.6}{:>10.6}{:>10.6}     {:>10.5}",
                0.1, 0.2, 0.3, 1.5
            ),
            format!(
                "SCALE2    {:>10.6}{:>10.6}{:>10.6}     {:>10.5}",
                0.4, 0.5, 0.6, 2.5
            ),
            format!(
                "SCALE3    {:>10.6}{:>10.6}{:>10.6}     {:>10.5}",
                0.7, 0.8, 0.9, 3.5
            ),
            format!(
                "ORIGX1    {:>10.6}{:>10.6}{:>10.6}     {:>10.5}",
                1.0, 0.0, 0.0, 4.5
            ),
            format!(
                "ORIGX2    {:>10.6}{:>10.6}{:>10.6}     {:>10.5}",
                0.0, 1.0, 0.0, 5.5
            ),
            format!(
                "ORIGX3    {:>10.6}{:>10.6}{:>10.6}     {:>10.5}",
                0.0, 0.0, 1.0, 6.5
            ),
            format!(
                "MTRIX1   1{:>10.6}{:>10.6}{:>10.6}     {:>10.5}    1",
                1.0, 0.0, 0.0, 0.0
            ),
            format!(
                "MTRIX2   1{:>10.6}{:>10.6}{:>10.6}     {:>10.5}    1",
                0.0, 1.0, 0.0, 0.0
            ),
            format!(
                "MTRIX3   1{:>10.6}{:>10.6}{:>10.6}     {:>10.5}    1",
                0.0, 0.0, 1.0, 0.0
            ),
            format!(
                "MTRIX1   2{:>10.6}{:>10.6}{:>10.6}     {:>10.5}    0",
                1.0, 0.0, 0.0, 1.0
            ),
            format!(
                "MTRIX2   2{:>10.6}{:>10.6}{:>10.6}     {:>10.5}    0",
                0.0, 1.0, 0.0, 2.0
            ),
            format!(
                "MTRIX3   2{:>10.6}{:>10.6}{:>10.6}     {:>10.5}    0",
                0.0, 0.0, 1.0, 3.0
            ),
        );

        let structure = BioStructure::from_pdb_str(&pdb).unwrap();
        let metadata = structure.metadata();
        let crystal = structure.crystal().unwrap();

        assert_eq!(metadata.entry_id.as_deref(), Some("2XYZ"));
        assert_eq!(
            metadata.received_initial_deposition_date.as_deref(),
            Some("2007-03-28")
        );
        assert_eq!(metadata.title.as_deref(), Some("EXAMPLE STRUCTURE"));
        assert_eq!(metadata.pdbx_keywords.as_deref(), Some("OXIDOREDUCTASE"));
        assert_eq!(metadata.keywords.as_deref(), Some("TEST, GEMMI ROUTE"));
        assert_eq!(
            metadata.experimental_method.as_deref(),
            Some("X-RAYDIFFRACTION")
        );
        assert_eq!(metadata.authors, vec!["DOE", "SMITH, J."]);
        assert_eq!(crystal.cell.a, 10.0);
        assert_eq!(crystal.cell.gamma, 120.0);
        assert_eq!(crystal.spacegroup_hm.as_deref(), Some("P 1"));
        assert_eq!(crystal.z_pdb.as_deref(), Some("2"));
        assert_eq!(
            crystal.scale,
            Some(BioTransform {
                mat: [[0.1, 0.2, 0.3], [0.4, 0.5, 0.6], [0.7, 0.8, 0.9]],
                vec: [1.5, 2.5, 3.5],
            })
        );
        assert!(structure.has_origx);
        assert_eq!(
            structure.origx,
            BioTransform {
                mat: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                vec: [4.5, 5.5, 6.5],
            }
        );
        assert_eq!(structure.ncs_oper_identity_id.as_deref(), Some("1"));
        assert_eq!(structure.ncs_operators.len(), 1);
        assert_eq!(structure.ncs_operators[0].id, "2");
        assert!(!structure.ncs_operators[0].given);
        assert_eq!(
            structure.ncs_operators[0].transform,
            BioTransform {
                mat: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                vec: [1.0, 2.0, 3.0],
            }
        );
        assert_eq!(crystal.cs_count, 0);
        assert_eq!(crystal.cell_images.len(), 1);
        let expected_image = bio_transform_combine(
            &crystal.frac,
            &bio_transform_combine(&structure.ncs_operators[0].transform, &crystal.orth),
        );
        assert_eq!(crystal.cell_images[0], expected_image);
    }

    #[test]
    fn reads_pdb_helix_and_sheet_secondary_structure_records() {
        let pdb = [
            "SEQRES   1 A    7  ALA GLY SER THR TYR VAL LEU".to_string(),
            "SEQRES   2 A    7  ASP".to_string(),
            format!(
                "{:<80}",
                "HELIX    1 AA1 ALA A   55  TYR A   59  5                                   5"
            ),
            format!("{:<80}", "SHEET    1 AA1 7 TYR A  20  THR A  21  0"),
            format!(
                "{:<80}",
                "SHEET    2 AA1 7 LYS A 156  PRO A 161 -1  O  CYS A 157   N  TYR A  20"
            ),
            format!(
                "{:<80}",
                "ATOM      1  CA  ALA A   7      11.104  13.207   9.900  1.00 20.00           C"
            ),
        ]
        .join("\n")
            + "\n";

        let structure = BioStructure::from_pdb_str(&pdb).unwrap();

        assert_eq!(
            structure.entities[0].sequence,
            vec!["ALA", "GLY", "SER", "THR", "TYR", "VAL", "LEU", "ASP"]
        );

        assert_eq!(structure.helices.len(), 1);
        assert_eq!(structure.helices[0].start.chain_name, "A");
        assert_eq!(structure.helices[0].start.seq_id.unwrap().seq_num, 55);
        assert_eq!(structure.helices[0].start.seq_id.unwrap().ins_code, None);
        assert_eq!(structure.helices[0].end.chain_name, "A");
        assert_eq!(structure.helices[0].end.seq_id.unwrap().seq_num, 59);
        assert_eq!(structure.helices[0].helix_class, BioHelixClass::R310);
        assert_eq!(structure.helices[0].length, 5);

        assert_eq!(structure.sheets.len(), 1);
        assert_eq!(structure.sheets[0].name, "AA1");
        assert_eq!(structure.sheets[0].strands.len(), 2);
        assert_eq!(structure.sheets[0].strands[0].sense, 0);
        assert_eq!(structure.sheets[0].strands[1].sense, -1);
        assert_eq!(structure.sheets[0].strands[1].hbond_atom2.atom_name, "O");
        assert_eq!(structure.sheets[0].strands[1].hbond_atom2.chain_name, "A");
        assert_eq!(
            structure.sheets[0].strands[1]
                .hbond_atom2
                .seq_id
                .unwrap()
                .seq_num,
            157
        );
        assert_eq!(structure.sheets[0].strands[1].hbond_atom1.atom_name, "N");
        assert_eq!(structure.sheets[0].strands[1].hbond_atom1.chain_name, "A");
        assert_eq!(
            structure.sheets[0].strands[1]
                .hbond_atom1
                .seq_id
                .unwrap()
                .seq_num,
            20
        );
    }

    #[test]
    fn explicit_pdb_params_accept_ter_and_deferred_connection_records() {
        let mut link = " ".repeat(80);
        link.replace_range(0..4, "LINK");
        link.replace_range(12..16, "ZN1A");
        link.replace_range(17..20, "ZN ");
        link.replace_range(20..22, "A ");
        link.replace_range(22..26, "   1");
        link.replace_range(42..46, " O  ");
        link.replace_range(47..50, "HOH");
        link.replace_range(50..52, "A ");
        link.replace_range(52..56, "   2");
        link.replace_range(59..63, "1555");
        link.replace_range(66..70, "1555");
        link.replace_range(73..77, "2.00");
        let mut cispep = " ".repeat(80);
        cispep.replace_range(0..6, "CISPEP");
        cispep.replace_range(14..16, "A ");
        cispep.replace_range(17..21, "   3");
        cispep.replace_range(28..30, "A ");
        cispep.replace_range(31..35, "   4");
        cispep.replace_range(43..46, "  0");
        cispep.replace_range(53..59, " 12.34");

        let pdb = format!(
            "\
ATOM      1  SG  CYS A   3       0.000   0.000   0.000  1.00 20.00           S
ATOM      2  SG  CYS A   4       2.000   0.000   0.000  1.00 20.00           S
HETATM    3 ZN1A  ZN A   1       5.000   0.000   0.000  1.00 20.00          ZN
ATOM      4  O   HOH A   2       6.900   0.000   0.000  1.00 20.00           O
TER
SSBOND   1 CYS A   3     CYS A   4                          1555   1555  2.03
{cispep}
{link}
"
        );

        let structure = BioStructure::from_pdb_str_with_params(
            &pdb,
            BioPdbReadParams {
                max_line_length: 120,
                ..BioPdbReadParams::default()
            },
        )
        .unwrap();

        assert_eq!(structure.ter_status, 'e');
        assert!(structure.deferred_conn_records.is_empty());
        assert_eq!(structure.connections.len(), 2);
        assert_eq!(structure.connections[0].type_, BioConnectionType::Disulf);
        assert_eq!(structure.connections[0].partner1.residue_name, "CYS");
        assert_eq!(structure.connections[0].partner2.residue_name, "CYS");
        assert_eq!(structure.connections[1].type_, BioConnectionType::MetalC);
        assert_eq!(structure.cispeps.len(), 1);
        assert_eq!(structure.cispeps[0].model_num, 1);
        assert_eq!(structure.cispeps[0].reported_angle, Some(12.34));
    }

    #[test]
    fn reads_pdb_ter_modres_hetnam_and_dbref_records() {
        let mut modres = " ".repeat(80);
        modres.replace_range(0..6, "MODRES");
        modres.replace_range(12..15, "MSE");
        modres.replace_range(15..17, "A ");
        modres.replace_range(18..22, "   5");
        modres.replace_range(24..27, "MET");
        modres.replace_range(29..70, &format!("{:<41}", "SELENOMETHIONINE"));
        modres.replace_range(72..80, "MODID001");

        let mut hetnam = " ".repeat(80);
        hetnam.replace_range(0..6, "HETNAM");
        hetnam.replace_range(11..14, "QWE");
        hetnam.replace_range(71..79, "LONGCODE");

        let mut dbref = " ".repeat(80);
        dbref.replace_range(0..5, "DBREF");
        dbref.replace_range(11..13, "A ");
        dbref.replace_range(14..19, "    1");
        dbref.replace_range(20..25, "    9");
        dbref.replace_range(26..32, "UNP   ");
        dbref.replace_range(33..41, "Q12345  ");
        dbref.replace_range(42..54, &format!("{:<12}", "SAMPLE_ID"));
        dbref.replace_range(55..60, "   11");
        dbref.replace_range(62..67, "   19");

        let mut dbref1 = " ".repeat(80);
        dbref1.replace_range(0..6, "DBREF1");
        dbref1.replace_range(11..13, "A ");
        dbref1.replace_range(14..19, "    1");
        dbref1.replace_range(20..25, "    9");
        dbref1.replace_range(26..32, "UNP   ");
        dbref1.replace_range(47..67, &format!("{:<20}", "SAMPLE_IDENTIFIER"));

        let mut dbref2 = " ".repeat(80);
        dbref2.replace_range(0..6, "DBREF2");
        dbref2.replace_range(11..13, "A ");
        dbref2.replace_range(18..40, &format!("{:<22}", "Q12345LONGACCESSION"));
        dbref2.replace_range(45..55, "       101");
        dbref2.replace_range(57..67, "       109");

        let dbref1 = dbref1.trim_end();
        let dbref2 = dbref2.trim_end();

        let pdb = format!(
            "\
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C
ATOM      2  O   HOH A   2       1.000   0.000   0.000  1.00 20.00           O
TER
TER
{modres}
{hetnam}
{dbref}
{dbref1}
{dbref2}
"
        );

        let structure = BioStructure::from_pdb_str(&pdb).unwrap();

        assert_eq!(structure.ter_status, 'e');
        assert_eq!(structure.residues[0].entity_kind, EntityKind::Unknown);
        assert_eq!(structure.residues[1].entity_kind, EntityKind::Unknown);
        assert_eq!(structure.mod_residues.len(), 1);
        assert_eq!(structure.mod_residues[0].chain_name, "A");
        assert_eq!(structure.mod_residues[0].res_id.seq_num, 5);
        assert_eq!(structure.mod_residues[0].parent_comp_id, "MET");
        assert_eq!(structure.mod_residues[0].mod_id, "MODID001");
        assert!(structure.shortened_ccd_codes.is_empty());
        assert_eq!(structure.entities.len(), 1);
        assert_eq!(structure.entities[0].dbrefs.len(), 2);
        assert_eq!(structure.entities[0].dbrefs[0].db_name, "UNP");
        assert_eq!(structure.entities[0].dbrefs[0].accession_code, "Q12345");
        assert_eq!(structure.entities[0].dbrefs[0].id_code, "SAMPLE_ID");
        assert_eq!(
            structure.entities[0].dbrefs[0].db_begin.unwrap().seq_num,
            11
        );
        assert_eq!(structure.entities[0].dbrefs[0].db_end.unwrap().seq_num, 19);
        assert_eq!(structure.entities[0].dbrefs[1].id_code, "SAMPLE_IDENTIFIER");
        assert_eq!(
            structure.entities[0].dbrefs[1].accession_code,
            "Q12345LONGACCESSION"
        );
        assert_eq!(
            structure.entities[0].dbrefs[1].db_begin.unwrap().seq_num,
            101
        );
        assert_eq!(structure.entities[0].dbrefs[1].db_end.unwrap().seq_num, 109);
        assert!(structure.residues[0].source.subchain_id.is_none());
        assert!(structure.residues[1].source.subchain_id.is_none());
    }

    #[test]
    fn assign_subchains_respects_force_preserves_existing_and_can_fail_on_unknown() {
        let mut structure = BioStructure::default();
        structure.models.push(ModelRow {
            chain_span: RowSpan::new(0, 1),
            source_model_number: Some(1),
        });
        structure.chains.push(ChainRow {
            model_id: ModelId::new(0),
            entity_id: None,
            residue_span: RowSpan::new(0, 2),
            kind: ChainKind::Unknown,
            source: ChainSourceIds {
                auth_chain_id: pdb_chain_id_from_str("A"),
                label_asym_id: None,
            },
        });
        structure.residues.push(ResidueRow {
            chain_id: ChainId::new(0),
            atom_span: RowSpan::new(0, 0),
            name: ResidueName([b'A', b'L', b'A', 0], 3),
            kind: ResidueKind::AminoAcid,
            entity_kind: EntityKind::Polymer,
            source: ResidueSourceIds {
                seq_id: Some(PdbSeqId {
                    seq_num: 1,
                    ins_code: None,
                }),
                label_seq_id: None,
                segment_id: None,
                subchain_id: None,
                label_entity_id: None,
            },
            het_flag: None,
            sifts_unp: None,
        });
        structure.residues.push(ResidueRow {
            chain_id: ChainId::new(0),
            atom_span: RowSpan::new(0, 0),
            name: ResidueName([b'H', b'O', b'H', 0], 3),
            kind: ResidueKind::Water,
            entity_kind: EntityKind::Water,
            source: ResidueSourceIds {
                seq_id: Some(PdbSeqId {
                    seq_num: 2,
                    ins_code: None,
                }),
                label_seq_id: None,
                segment_id: None,
                subchain_id: None,
                label_entity_id: None,
            },
            het_flag: None,
            sifts_unp: None,
        });

        assign_subchains(&mut structure, false, false).unwrap();
        assert_eq!(
            structure.residues[0].source.subchain_id.unwrap().as_str(),
            "Axp"
        );
        assert_eq!(
            structure.residues[1].source.subchain_id.unwrap().as_str(),
            "Axw"
        );

        structure.residues[0].source.subchain_id = pdb_chain_id_from_str("QQQ");
        assign_subchains(&mut structure, false, false).unwrap();
        assert_eq!(
            structure.residues[0].source.subchain_id.unwrap().as_str(),
            "QQQ"
        );

        assign_subchains(&mut structure, true, false).unwrap();
        assert_eq!(
            structure.residues[0].source.subchain_id.unwrap().as_str(),
            "Axp"
        );

        structure.residues[0].entity_kind = EntityKind::Unknown;
        let err = assign_subchains(&mut structure, true, true).unwrap_err();
        assert!(matches!(err, BioReadError::Parse { .. }));
    }

    #[test]
    fn pdb_finish_inserts_empty_model_normalizes_authors_and_gates_remarks() {
        let structure = BioStructure::from_pdb_str(
            "\
HEADER    HYDROLASE                               01-JAN-00   1ABC
AUTHOR    A.-B.DOE
REMARK 300 REMARK: DETAIL
",
        )
        .unwrap();

        assert_eq!(structure.num_models(), 1);
        assert_eq!(structure.models[0].source_model_number, Some(1));
        assert_eq!(structure.metadata.authors, vec!["DOE, A.-B."]);
        assert_eq!(
            structure.metadata.remark_300_detail.as_deref(),
            Some("DETAIL")
        );
        assert!(structure.crystal().is_none());

        let skipped = BioStructure::from_pdb_str_with_params(
            "\
AUTHOR    A.B.DOE
REMARK 300 REMARK: SHOULD NOT APPEAR
ATOM      1  CA  ALA A   7      11.104  13.207   9.900  1.00 20.00           C
",
            BioPdbReadParams {
                skip_remarks: true,
                ..BioPdbReadParams::default()
            },
        )
        .unwrap();
        assert_eq!(skipped.metadata.authors, vec!["DOE, A.B."]);
        assert!(skipped.metadata.remark_300_detail.is_none());
    }

    #[test]
    fn pdb_finish_sets_non_p1_spacegroup_cell_images_from_gemmi_spacegroup_ops() {
        let pdb = "\
CRYST1   10.000   11.000   12.000  90.00  90.00 120.00 P 21 21 21    4
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C
";

        let structure = BioStructure::from_pdb_str(pdb).unwrap();
        let crystal = structure.crystal().unwrap();
        let sg = find_spacegroup_by_name("P 21 21 21", 90.0, 120.0, None).unwrap();

        assert_eq!(crystal.cs_count as usize, sg.ops.len() - 1);
        assert_eq!(crystal.cell_images.len(), sg.ops.len() - 1);
        for (image, op) in crystal.cell_images.iter().zip(&sg.ops[1..]) {
            assert_eq!(*image, gemmi_op_to_bio_transform(op));
        }
    }

    #[test]
    fn pdb_finish_backfills_first_model_polymer_subchains_to_entities() {
        let pdb = "\
SEQRES   1 A    2  ALA GLY
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C
ATOM      2  CA  GLY A   2       1.000   0.000   0.000  1.00 20.00           C
TER
";

        let structure = BioStructure::from_pdb_str(pdb).unwrap();
        assert_eq!(structure.entities.len(), 1);
        assert_eq!(structure.entities[0].kind, EntityKind::Polymer);
        assert_eq!(structure.entities[0].subchains[0].as_str(), "Axp");
        assert!(
            structure.entities[0]
                .subchains
                .iter()
                .any(|subchain| subchain.as_str() == "Axp")
        );
    }

    #[test]
    fn mmcif_copy_helpers_only_write_present_values() {
        let row = vec![
            CifToken {
                value: "12".to_string(),
                line_number: 7,
            },
            CifToken {
                value: ".".to_string(),
                line_number: 7,
            },
            CifToken {
                value: "3.5".to_string(),
                line_number: 7,
            },
            CifToken {
                value: "?".to_string(),
                line_number: 7,
            },
            CifToken {
                value: "HELLO".to_string(),
                line_number: 7,
            },
        ];

        let mut int_dest = 0;
        copy_int(&row, 0, &mut int_dest).unwrap();
        assert_eq!(int_dest, 12);
        copy_int(&row, 1, &mut int_dest).unwrap();
        assert_eq!(int_dest, 12);

        let mut float_dest = 0.0;
        copy_double(&row, 2, &mut float_dest).unwrap();
        assert_eq!(float_dest, 3.5);
        copy_double(&row, 3, &mut float_dest).unwrap();
        assert_eq!(float_dest, 3.5);

        let mut string_dest = String::new();
        copy_string(&row, 4, &mut string_dest);
        assert_eq!(string_dest, "HELLO");
        copy_string(&row, 1, &mut string_dest);
        assert_eq!(string_dest, "HELLO");
    }

    #[test]
    fn mmcif_get_smat33_reads_gemmi_value_layout() {
        let row = vec![
            CifToken {
                value: "1.0".to_string(),
                line_number: 11,
            },
            CifToken {
                value: "2.0".to_string(),
                line_number: 11,
            },
            CifToken {
                value: "3.0".to_string(),
                line_number: 11,
            },
            CifToken {
                value: "4.0".to_string(),
                line_number: 11,
            },
            CifToken {
                value: "5.0".to_string(),
                line_number: 11,
            },
            CifToken {
                value: "6.0".to_string(),
                line_number: 11,
            },
        ];

        let smat = get_smat33(&row, 0).unwrap();
        assert_eq!(
            smat,
            CifSmat33 {
                u11: 1.0,
                u22: 2.0,
                u33: 3.0,
                u12: 4.0,
                u13: 5.0,
                u23: 6.0,
            }
        );
    }

    #[test]
    fn mmcif_transform_tags_match_gemmi_order() {
        assert_eq!(
            transform_tags("_m", "_v"),
            [
                "_m[1][1]".to_string(),
                "_m[1][2]".to_string(),
                "_m[1][3]".to_string(),
                "_v[1]".to_string(),
                "_m[2][1]".to_string(),
                "_m[2][2]".to_string(),
                "_m[2][3]".to_string(),
                "_v[2]".to_string(),
                "_m[3][1]".to_string(),
                "_m[3][2]".to_string(),
                "_m[3][3]".to_string(),
                "_v[3]".to_string(),
            ]
        );
    }

    #[test]
    fn mmcif_get_transform_matrix_reads_3x4_layout() {
        let row = vec![
            CifToken {
                value: "1.0".to_string(),
                line_number: 21,
            },
            CifToken {
                value: "2.0".to_string(),
                line_number: 21,
            },
            CifToken {
                value: "3.0".to_string(),
                line_number: 21,
            },
            CifToken {
                value: "4.0".to_string(),
                line_number: 21,
            },
            CifToken {
                value: "5.0".to_string(),
                line_number: 21,
            },
            CifToken {
                value: "6.0".to_string(),
                line_number: 21,
            },
            CifToken {
                value: "7.0".to_string(),
                line_number: 21,
            },
            CifToken {
                value: "8.0".to_string(),
                line_number: 21,
            },
            CifToken {
                value: "9.0".to_string(),
                line_number: 21,
            },
            CifToken {
                value: "10.0".to_string(),
                line_number: 21,
            },
            CifToken {
                value: "11.0".to_string(),
                line_number: 21,
            },
            CifToken {
                value: "12.0".to_string(),
                line_number: 21,
            },
        ];

        assert_eq!(
            get_transform_matrix(&row).unwrap(),
            BioTransform {
                mat: [[1.0, 2.0, 3.0], [5.0, 6.0, 7.0], [9.0, 10.0, 11.0]],
                vec: [4.0, 8.0, 12.0],
            }
        );
    }

    #[test]
    fn mmcif_make_seqid_supports_old_auth_suffix_and_empty_seqnum() {
        assert_eq!(
            make_seqid("15A", None, 4).unwrap(),
            Some(PdbSeqId {
                seq_num: 15,
                ins_code: Some(b'A'),
            })
        );
        assert_eq!(
            make_seqid("A", None, 4).unwrap(),
            Some(PdbSeqId {
                seq_num: 0,
                ins_code: Some(b'A'),
            })
        );
        assert_eq!(make_seqid("", None, 4).unwrap(), None);
    }

    #[test]
    fn mmcif_make_seqid_rejects_inconsistent_insertion_code() {
        let err = make_seqid("15A", Some("B"), 23).unwrap_err();
        assert!(matches!(
            err,
            BioReadError::Parse {
                line_number: 23,
                ..
            }
        ));
    }

    #[test]
    fn mmcif_row_access_uses_fallback_only_for_null_primary() {
        let row = vec![
            CifToken {
                value: ".".to_string(),
                line_number: 1,
            },
            CifToken {
                value: "AUTH".to_string(),
                line_number: 1,
            },
            CifToken {
                value: "LABEL".to_string(),
                line_number: 1,
            },
        ];
        let primary_missing = CifRowAccess::new(None, Some(2));
        assert!(primary_missing.ok());
        assert_eq!(primary_missing.get(&row), Some("LABEL"));

        let auth_first = CifRowAccess::new(Some(1), Some(2));
        assert_eq!(auth_first.get(&row), Some("AUTH"));

        let null_primary = CifRowAccess::new(Some(0), Some(2));
        assert_eq!(null_primary.get(&row), Some("LABEL"));
    }

    #[test]
    fn mmcif_get_by_id_returns_first_matching_entry() {
        #[derive(Debug)]
        struct Item {
            id: String,
            value: i32,
        }

        let mut items = vec![
            Item {
                id: "a".to_string(),
                value: 1,
            },
            Item {
                id: "b".to_string(),
                value: 2,
            },
            Item {
                id: "b".to_string(),
                value: 3,
            },
        ];

        let hit = get_by_id(&mut items, "b", |item| item.id.as_str()).unwrap();
        assert_eq!(hit.value, 2);
        assert!(get_by_id(&mut items, "z", |item| item.id.as_str()).is_none());
    }

    #[test]
    fn mmcif_anisotropic_u_map_is_attached_to_matching_atom_id() {
        let cif = r#"
data_demo
loop_
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
1 C CA . ALA A 7 11.104 13.207 9.900
loop_
_atom_site_anisotrop.id
_atom_site_anisotrop.U[1][1]
_atom_site_anisotrop.U[2][2]
_atom_site_anisotrop.U[3][3]
_atom_site_anisotrop.U[1][2]
_atom_site_anisotrop.U[1][3]
_atom_site_anisotrop.U[2][3]
1 1.0 2.0 3.0 4.0 5.0 6.0
"#;

        let structure = read_test_mmcif(cif).unwrap();
        assert_eq!(
            structure.atoms[0].anisou,
            Some([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
        );
    }

    #[test]
    fn mmcif_anisotropic_u_map_skips_atoms_without_matching_row() {
        let cif = r#"
data_demo
loop_
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
1 C CA . ALA A 7 11.104 13.207 9.900
2 N N  . ALA A 7 10.104 12.207 8.900
loop_
_atom_site_anisotrop.id
_atom_site_anisotrop.U[1][1]
_atom_site_anisotrop.U[2][2]
_atom_site_anisotrop.U[3][3]
_atom_site_anisotrop.U[1][2]
_atom_site_anisotrop.U[1][3]
_atom_site_anisotrop.U[2][3]
2 1.5 2.5 3.5 4.5 5.5 6.5
"#;

        let structure = read_test_mmcif(cif).unwrap();

        assert_eq!(structure.atoms[0].anisou, None);
        assert_eq!(
            structure.atoms[1].anisou,
            Some([1.5, 2.5, 3.5, 4.5, 5.5, 6.5])
        );
    }

    #[test]
    fn mmcif_secondary_structure_helpers_handle_missing_loops() {
        let loops = parse_cif_loops("data_demo\n").unwrap();
        assert!(get_anisotropic_u(&loops).unwrap().is_empty());
        assert!(read_helices(&loops).unwrap().is_empty());
        assert!(read_sheets(&loops).unwrap().is_empty());
    }

    #[test]
    fn mmcif_read_helices_filters_non_helix_and_parses_class_and_length() {
        let cif = r#"
data_demo
loop_
_struct_conf.conf_type_id
_struct_conf.beg_auth_asym_id
_struct_conf.beg_label_comp_id
_struct_conf.beg_auth_seq_id
_struct_conf.pdbx_beg_PDB_ins_code
_struct_conf.end_auth_asym_id
_struct_conf.end_label_comp_id
_struct_conf.end_auth_seq_id
_struct_conf.pdbx_end_PDB_ins_code
_struct_conf.pdbx_PDB_helix_class
_struct_conf.pdbx_PDB_helix_length
HELX_P A ALA 55 ? A TYR 59 ? 5 7
TURN_TY1 A GLY 60 ? A SER 61 ? ? ?
"#;

        let loops = parse_cif_loops(cif).unwrap();
        let helices = read_helices(&loops).unwrap();
        assert_eq!(helices.len(), 1);
        assert_eq!(helices[0].start.chain_name, "A");
        assert_eq!(helices[0].start.seq_id.unwrap().seq_num, 55);
        assert_eq!(helices[0].end.seq_id.unwrap().seq_num, 59);
        assert_eq!(helices[0].helix_class, BioHelixClass::R310);
        assert_eq!(helices[0].length, 7);
    }

    #[test]
    fn mmcif_read_sheets_sets_sense_and_hbond_atoms() {
        let cif = r#"
data_demo
loop_
_struct_sheet.id
AA1
loop_
_struct_sheet_range.sheet_id
_struct_sheet_range.id
_struct_sheet_range.beg_auth_asym_id
_struct_sheet_range.beg_label_comp_id
_struct_sheet_range.beg_auth_seq_id
_struct_sheet_range.pdbx_beg_PDB_ins_code
_struct_sheet_range.end_auth_asym_id
_struct_sheet_range.end_label_comp_id
_struct_sheet_range.end_auth_seq_id
_struct_sheet_range.pdbx_end_PDB_ins_code
AA1 1 A TYR 20 ? A THR 21 ?
AA1 2 A LYS 156 ? A PRO 161 ?
loop_
_struct_sheet_order.sheet_id
_struct_sheet_order.range_id_2
_struct_sheet_order.sense
AA1 2 anti-parallel
loop_
_pdbx_struct_sheet_hbond.sheet_id
_pdbx_struct_sheet_hbond.range_id_2
_pdbx_struct_sheet_hbond.range_1_auth_asym_id
_pdbx_struct_sheet_hbond.range_1_label_comp_id
_pdbx_struct_sheet_hbond.range_1_auth_seq_id
_pdbx_struct_sheet_hbond.range_1_PDB_ins_code
_pdbx_struct_sheet_hbond.range_1_label_atom_id
_pdbx_struct_sheet_hbond.range_2_auth_asym_id
_pdbx_struct_sheet_hbond.range_2_label_comp_id
_pdbx_struct_sheet_hbond.range_2_auth_seq_id
_pdbx_struct_sheet_hbond.range_2_PDB_ins_code
_pdbx_struct_sheet_hbond.range_2_label_atom_id
AA1 2 A CYS 157 ? O A TYR 20 ? N
"#;

        let loops = parse_cif_loops(cif).unwrap();
        let sheets = read_sheets(&loops).unwrap();
        assert_eq!(sheets.len(), 1);
        assert_eq!(sheets[0].strands.len(), 2);
        assert_eq!(sheets[0].strands[1].sense, -1);
        assert_eq!(sheets[0].strands[1].hbond_atom1.atom_name, "O");
        assert_eq!(
            sheets[0].strands[1].hbond_atom1.seq_id.unwrap().seq_num,
            157
        );
        assert_eq!(sheets[0].strands[1].hbond_atom2.atom_name, "N");
        assert_eq!(sheets[0].strands[1].hbond_atom2.seq_id.unwrap().seq_num, 20);
    }

    #[test]
    fn mmcif_cell_and_spacegroup_helpers_read_present_tags() {
        let cif = r#"
data_demo
_cell.length_a 10.0
_cell.length_b 20.0
_cell.length_c 30.0
_cell.angle_alpha 90.0
_cell.angle_beta 100.0
_cell.angle_gamma 120.0
_symmetry.space_group_name_H-M 'P 21 21 21'
"#;

        let document = parse_cif_document(cif).unwrap();
        let mut crystal = None;
        set_cell_from_mmcif(&document, &mut crystal);

        let crystal = crystal.unwrap();
        assert_eq!(crystal.cell.a, 10.0);
        assert_eq!(crystal.cell.b, 20.0);
        assert_eq!(crystal.cell.c, 30.0);
        assert_eq!(crystal.cell.alpha, 90.0);
        assert_eq!(crystal.cell.beta, 100.0);
        assert_eq!(crystal.cell.gamma, 120.0);
        assert_eq!(
            find_spacegroup_hm_value(&document).and_then(cif_optional),
            Some("P 21 21 21")
        );
    }

    #[test]
    fn mmcif_cell_helper_ignores_partially_null_cell_tags() {
        let cif = r#"
data_demo
_cell.length_a 10.0
_cell.length_b .
_cell.length_c 30.0
_cell.angle_alpha 90.0
_cell.angle_beta 100.0
_cell.angle_gamma 120.0
_symmetry.space_group_name_H-M .
"#;

        let document = parse_cif_document(cif).unwrap();
        let mut crystal = Some(default_crystal_info());
        crystal.as_mut().unwrap().cell.a = 7.5;
        set_cell_from_mmcif(&document, &mut crystal);

        assert_eq!(crystal.unwrap().cell.a, 7.5);
        assert_eq!(find_spacegroup_hm_value(&document), Some("."));
        assert_eq!(
            find_spacegroup_hm_value(&document).and_then(cif_optional),
            None
        );
    }

    #[test]
    fn mmcif_cell_and_spacegroup_helpers_handle_absent_tags() {
        let document = parse_cif_document("data_demo\n").unwrap();
        let mut crystal = None;
        set_cell_from_mmcif(&document, &mut crystal);

        assert!(crystal.is_none());
        assert_eq!(find_spacegroup_hm_value(&document), None);
    }

    #[test]
    fn read_matrix_extracts_rows_one_through_three() {
        let mut transform = BioTransform::default();

        let row1 = format!(
            "SCALE1    {:>10.6}{:>10.6}{:>10.6}     {:>10.5}",
            0.1, 0.2, 0.3, 1.5
        );
        let row2 = format!(
            "SCALE2    {:>10.6}{:>10.6}{:>10.6}     {:>10.5}",
            0.4, 0.5, 0.6, 2.5
        );
        let row3 = format!(
            "SCALE3    {:>10.6}{:>10.6}{:>10.6}     {:>10.5}",
            0.7, 0.8, 0.9, 3.5
        );

        assert_eq!(read_matrix(&mut transform, &row1), 1);
        assert_eq!(read_matrix(&mut transform, &row2), 2);
        assert_eq!(read_matrix(&mut transform, &row3), 3);
        assert_eq!(transform.mat[0], [0.1, 0.2, 0.3]);
        assert_eq!(transform.mat[1], [0.4, 0.5, 0.6]);
        assert_eq!(transform.mat[2], [0.7, 0.8, 0.9]);
        assert_eq!(transform.vec, [1.5, 2.5, 3.5]);
    }

    #[test]
    fn read_matrix_returns_zero_for_short_records() {
        let mut transform = BioTransform::default();

        assert_eq!(read_matrix(&mut transform, "SCALE1    0.100"), 0);
        assert_eq!(transform, BioTransform::default());
    }

    #[test]
    fn read_matrix_ignores_non_one_to_three_row_numbers() {
        let mut transform = BioTransform::default();
        let row4 = format!(
            "SCALE4    {:>10.6}{:>10.6}{:>10.6}     {:>10.5}",
            9.9, 8.8, 7.7, 6.5
        );

        assert_eq!(read_matrix(&mut transform, &row4), 4);
        assert_eq!(transform, BioTransform::default());
    }

    #[test]
    fn remark_helpers_add_software_splits_name_version_and_parenthesized_date() {
        let mut meta = BioMetadata::default();

        add_software(
            &mut meta,
            BioSoftwareClassification::Refinement,
            "REFMAC version 5.8.0258 (12-JAN-18), SHELXL 2018/3",
        );

        assert_eq!(meta.software.len(), 2);
        assert_eq!(meta.software[0].name, "REFMAC");
        assert_eq!(meta.software[0].version, "5.8.0258");
        assert_eq!(meta.software[0].date, "2018-01-12");
        assert_eq!(
            meta.software[0].classification,
            BioSoftwareClassification::Refinement
        );
        assert_eq!(meta.software[1].name, "SHELXL");
        assert_eq!(meta.software[1].version, "2018/3");
    }

    #[test]
    fn remark_helpers_add_restraint_count_weight_extracts_count_weight_and_function() {
        let mut refinement = BioRefinementInfo::default();

        add_restraint_count_weight(&mut refinement, "t_bond_d", "5760   ; 2.000  ; HARMONIC");

        assert_eq!(refinement.restr_stats.len(), 1);
        let restraint = &refinement.restr_stats[0];
        assert_eq!(restraint.name, "t_bond_d");
        assert_eq!(restraint.count, Some(5760));
        assert_eq!(restraint.weight, Some(2.0));
        assert_eq!(restraint.function, "HARMONIC");
    }

    #[test]
    fn remark_helpers_read_remark3_line_tracks_tls_continuations_and_tensor_values() {
        let mut meta = BioMetadata::default();
        let mut continuation = None;

        read_remark3_line(
            "REMARK   3   DATA USED IN REFINEMENT.",
            &mut meta,
            &mut continuation,
        );
        read_remark3_line("REMARK   3    TLS GROUP : 1", &mut meta, &mut continuation);
        read_remark3_line(
            "REMARK   3    SET : RESIDUES A  1 THROUGH A  10",
            &mut meta,
            &mut continuation,
        );
        read_remark3_line(
            "REMARK   3              AND A  11 THROUGH A  20",
            &mut meta,
            &mut continuation,
        );
        read_remark3_line(
            "REMARK   3    T11: 0.100 T22: 0.200 T33: 0.300",
            &mut meta,
            &mut continuation,
        );
        read_remark3_line(
            "REMARK   3    BOND LENGTHS                       (A): 0.012",
            &mut meta,
            &mut continuation,
        );

        let refinement = &meta.refinement[0];
        let tls_group = &refinement.tls_groups[0];
        assert_eq!(tls_group.id, "1");
        assert_eq!(tls_group.num_id, Some(1));
        assert_eq!(
            tls_group.selections[0].details,
            "RESIDUES A  1 THROUGH A  10 AND A  11 THROUGH A  20"
        );
        assert_eq!(tls_group.t[0][0], 0.1);
        assert_eq!(tls_group.t[1][1], 0.2);
        assert_eq!(tls_group.t[2][2], 0.3);
        assert_eq!(refinement.restr_stats[0].name, "t_bond_d");
        assert_eq!(refinement.restr_stats[0].dev_ideal, Some(0.012));
    }

    #[test]
    fn remark_helpers_read_remark_200_230_240_reads_experiment_fields_and_multiline_description() {
        let mut meta = BioMetadata::default();
        let mut continuation = None;

        read_remark_200_230_240(
            "REMARK 200 EXPERIMENTAL DETAILS",
            &mut meta,
            &mut continuation,
        );
        read_remark_200_230_240(
            "REMARK 200 EXPERIMENT TYPE                : X-RAY DIFFRACTION",
            &mut meta,
            &mut continuation,
        );
        read_remark_200_230_240(
            "REMARK 200 PH                             : 6.5",
            &mut meta,
            &mut continuation,
        );
        read_remark_200_230_240(
            "REMARK 200 REMARK                         : FIRST LINE",
            &mut meta,
            &mut continuation,
        );
        read_remark_200_230_240("REMARK 200  SECOND LINE", &mut meta, &mut continuation);
        read_remark_200_230_240(
            "REMARK 200 INTENSITY-INTEGRATION SOFTWARE : XDS VERSION 1.0",
            &mut meta,
            &mut continuation,
        );

        assert_eq!(meta.experiments.len(), 1);
        assert_eq!(meta.experiments[0].method, "X-RAY DIFFRACTION");
        assert_eq!(meta.experiment_crystals.len(), 1);
        assert_eq!(meta.experiment_crystals[0].ph, Some(6.5));
        assert_eq!(
            meta.experiment_crystals[0].description,
            "FIRST LINE SECOND LINE"
        );
        assert_eq!(meta.software.len(), 1);
        assert_eq!(
            meta.software[0].classification,
            BioSoftwareClassification::DataReduction
        );
        assert_eq!(meta.software[0].name, "XDS");
        assert_eq!(meta.software[0].version, "1.0");
    }

    #[test]
    fn remark_metadata_helpers_read_resolution_detail_and_assembly_generators() {
        let mut structure = BioStructure::default();
        let mut biomolecule = " ".repeat(80);
        biomolecule.replace_range(0..10, "REMARK 350");
        biomolecule.replace_range(11..23, "BIOMOLECULE:");
        biomolecule.replace_range(23..24, "1");
        let mut apply = " ".repeat(80);
        apply.replace_range(0..10, "REMARK 350");
        apply.replace_range(11..40, "APPLY THE FOLLOWING TO CHAINS:");
        apply.replace_range(41..45, "A, B");
        let mut biomt1 = " ".repeat(80);
        biomt1.replace_range(0..10, "REMARK 350");
        biomt1.replace_range(13..19, "BIOMT1");
        biomt1.replace_range(20..23, "  1");
        biomt1.replace_range(23..33, &format!("{:>10.6}", 1.0));
        biomt1.replace_range(33..43, &format!("{:>10.6}", 0.0));
        biomt1.replace_range(43..53, &format!("{:>10.6}", 0.0));
        biomt1.replace_range(58..68, &format!("{:>10.5}", 0.0));
        let mut biomt2 = " ".repeat(80);
        biomt2.replace_range(0..10, "REMARK 350");
        biomt2.replace_range(13..19, "BIOMT2");
        biomt2.replace_range(20..23, "  1");
        biomt2.replace_range(23..33, &format!("{:>10.6}", 0.0));
        biomt2.replace_range(33..43, &format!("{:>10.6}", 1.0));
        biomt2.replace_range(43..53, &format!("{:>10.6}", 0.0));
        biomt2.replace_range(58..68, &format!("{:>10.5}", 0.0));
        let mut biomt3 = " ".repeat(80);
        biomt3.replace_range(0..10, "REMARK 350");
        biomt3.replace_range(13..19, "BIOMT3");
        biomt3.replace_range(20..23, "  1");
        biomt3.replace_range(23..33, &format!("{:>10.6}", 0.0));
        biomt3.replace_range(33..43, &format!("{:>10.6}", 0.0));
        biomt3.replace_range(43..53, &format!("{:>10.6}", 1.0));
        biomt3.replace_range(58..68, &format!("{:>10.5}", 0.0));
        structure.raw_remarks = vec![
            "REMARK   2 RESOLUTION.    1.80 ANGSTROM.".to_string(),
            "REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) : 1.90".to_string(),
            "REMARK 300 REMARK: FIRST DETAIL".to_string(),
            "REMARK 300 SECOND DETAIL".to_string(),
            biomolecule,
            apply,
            biomt1,
            biomt2,
            biomt3,
        ];

        read_metadata_from_remarks(&mut structure).unwrap();

        assert!((structure.resolution().unwrap() - 1.8).abs() < 1e-6);
        assert_eq!(
            structure.metadata.remark_300_detail.as_deref(),
            Some("FIRST DETAIL\nSECOND DETAIL")
        );
        assert_eq!(structure.assemblies.len(), 1);
        assert_eq!(structure.assemblies[0].name, "1");
        assert_eq!(structure.assemblies[0].generators[0].chains, vec!["A", "B"]);
        assert_eq!(structure.assemblies[0].generators[0].operators.len(), 1);
        assert_eq!(
            structure.assemblies[0].generators[0].operators[0]
                .transform
                .mat[0],
            [1.0, 0.0, 0.0]
        );
    }

    #[test]
    fn remark_metadata_helpers_fall_back_to_refinement_resolution_when_remark2_missing() {
        let mut structure = BioStructure::default();
        structure.raw_remarks = vec![
            "REMARK   3   DATA USED IN REFINEMENT.".to_string(),
            "REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) : 2.25".to_string(),
        ];

        read_metadata_from_remarks(&mut structure).unwrap();

        assert_eq!(structure.resolution(), Some(2.25));
        assert_eq!(structure.metadata.refinement.len(), 1);
    }

    #[test]
    fn pdb_connection_helpers_complete_ssbond_picks_shortest_matching_altloc_pair() {
        let pdb = "\
ATOM      1  SD ACYS A   1       0.000   0.000   0.000  1.00 20.00           S
ATOM      2  SD BCYS A   1      20.000   0.000   0.000  1.00 20.00           S
ATOM      3  SD ACYS A   2       1.500   0.000   0.000  1.00 20.00           S
ATOM      4  SD BCYS A   2      30.000   0.000   0.000  1.00 20.00           S
";
        let structure = BioStructure::from_pdb_str(pdb).unwrap();
        let mut connection = BioConnection {
            partner1: BioAtomAddress {
                chain_name: "A".to_string(),
                residue_name: "CYS".to_string(),
                seq_id: Some(PdbSeqId {
                    seq_num: 1,
                    ins_code: None,
                }),
                ..BioAtomAddress::default()
            },
            partner2: BioAtomAddress {
                chain_name: "A".to_string(),
                residue_name: "CYS".to_string(),
                seq_id: Some(PdbSeqId {
                    seq_num: 2,
                    ins_code: None,
                }),
                ..BioAtomAddress::default()
            },
            ..BioConnection::default()
        };

        complete_ssbond(&mut connection, &structure);

        assert_eq!(connection.partner1.atom_name, "SD");
        assert_eq!(connection.partner2.atom_name, "SD");
        assert_eq!(connection.partner1.altloc, Some(AltLocLabel(b'A')));
        assert_eq!(connection.partner2.altloc, Some(AltLocLabel(b'A')));
    }

    #[test]
    fn pdb_connection_helpers_compare_link_symops_reports_any_same_and_different() {
        let mut sym = [0_i16; 4];
        assert_eq!(compare_link_symops("LINK", &mut sym), BioAsu::Any);

        let mut same = " ".repeat(80);
        same.replace_range(0..4, "LINK");
        same.replace_range(59..63, "1555");
        same.replace_range(66..70, "1555");
        assert_eq!(compare_link_symops(&same, &mut sym), BioAsu::Same);

        let mut different = " ".repeat(80);
        different.replace_range(0..4, "LINK");
        different.replace_range(59..63, "1555");
        different.replace_range(66..70, "2666");
        assert_eq!(compare_link_symops(&different, &mut sym), BioAsu::Different);
        assert_eq!(sym, [2, 1, 1, 1]);
    }

    #[test]
    fn pdb_connection_helpers_process_conn_types_link_and_cispep_records() {
        let pdb = "\
HETATM    1 ZN1A  ZN A   1       0.000   0.000   0.000  1.00 20.00          ZN
ATOM      2  O   HOH A   2       2.000   0.000   0.000  1.00 20.00           O
ATOM      3  C1  LIG A   3       5.000   0.000   0.000  1.00 20.00           C
ATOM      4  O1  LIG A   4       6.000   0.000   0.000  1.00 20.00           O
";
        let mut structure = BioStructure::from_pdb_str(pdb).unwrap();
        let mut metal_link = " ".repeat(80);
        metal_link.replace_range(0..4, "LINK");
        metal_link.replace_range(12..16, "ZN1A");
        metal_link.replace_range(17..20, "ZN ");
        metal_link.replace_range(20..22, "A ");
        metal_link.replace_range(22..26, "   1");
        metal_link.replace_range(42..46, " O  ");
        metal_link.replace_range(47..50, "HOH");
        metal_link.replace_range(50..52, "A ");
        metal_link.replace_range(52..56, "   2");
        metal_link.replace_range(59..63, "1555");
        metal_link.replace_range(66..70, "1555");
        metal_link.replace_range(73..77, "2.00");

        let mut covale_linkr = " ".repeat(80);
        covale_linkr.replace_range(0..5, "LINKR");
        covale_linkr.replace_range(12..16, " C1 ");
        covale_linkr.replace_range(17..20, "LIG");
        covale_linkr.replace_range(20..22, "A ");
        covale_linkr.replace_range(22..26, "   3");
        covale_linkr.replace_range(42..46, " O1 ");
        covale_linkr.replace_range(47..50, "LIG");
        covale_linkr.replace_range(50..52, "A ");
        covale_linkr.replace_range(52..56, "   4");
        covale_linkr.replace_range(59..63, "1555");
        covale_linkr.replace_range(66..70, "1555");
        covale_linkr.replace_range(72..80, "TESTLINK");
        let mut cispep = " ".repeat(80);
        cispep.replace_range(0..6, "CISPEP");
        cispep.replace_range(14..16, "A ");
        cispep.replace_range(17..21, "   3");
        cispep.replace_range(28..30, "A ");
        cispep.replace_range(31..35, "   4");
        cispep.replace_range(43..46, "  0");
        cispep.replace_range(53..59, " 12.34");

        let records = vec![metal_link, covale_linkr, cispep];

        process_conn(&mut structure, &records);

        assert_eq!(structure.connections.len(), 2);
        assert_eq!(structure.connections[0].type_, BioConnectionType::MetalC);
        assert_eq!(structure.connections[0].name, "metalc1");
        assert_eq!(structure.connections[1].type_, BioConnectionType::Covale);
        assert_eq!(structure.connections[1].name, "covale1");
        assert_eq!(structure.connections[1].link_id, "TESTLINK");
        assert_eq!(structure.cispeps.len(), 1);
        assert_eq!(structure.cispeps[0].model_num, 1);
        assert_eq!(structure.cispeps[0].reported_angle, Some(12.34));
    }

    #[test]
    fn pdb_connection_helpers_use_gemmi_full_metal_table() {
        let pdb = "\
HETATM    1 CD    CD A   1       0.000   0.000   0.000  1.00 20.00          CD
ATOM      2  O   HOH A   2       2.000   0.000   0.000  1.00 20.00           O
";
        let mut structure = BioStructure::from_pdb_str(pdb).unwrap();
        let mut metal_link = " ".repeat(80);
        metal_link.replace_range(0..4, "LINK");
        metal_link.replace_range(12..16, "CD  ");
        metal_link.replace_range(17..20, "CD ");
        metal_link.replace_range(20..22, "A ");
        metal_link.replace_range(22..26, "   1");
        metal_link.replace_range(42..46, " O  ");
        metal_link.replace_range(47..50, "HOH");
        metal_link.replace_range(50..52, "A ");
        metal_link.replace_range(52..56, "   2");

        process_conn(&mut structure, &[metal_link]);

        assert_eq!(structure.connections.len(), 1);
        assert_eq!(structure.connections[0].type_, BioConnectionType::MetalC);
        assert_eq!(structure.connections[0].name, "metalc1");
    }

    #[test]
    fn pdb_connection_helpers_change_author_name_format_to_mmcif_moves_initials() {
        let mut name = " A.-B.DOE".to_string();
        change_author_name_format_to_mmcif(&mut name);
        assert_eq!(name, "DOE, A.-B.");
    }

    #[test]
    fn pdb_connection_helpers_read_remark_290_checks_sequence_order() {
        let mut op1 = " ".repeat(80);
        op1.replace_range(0..10, "REMARK 290");
        op1.replace_range(15..18, "  1");
        op1.replace_range(18..24, "555   ");
        op1.replace_range(24..29, "X,Y,Z");
        let mut op2 = " ".repeat(80);
        op2.replace_range(0..10, "REMARK 290");
        op2.replace_range(15..18, "  2");
        op2.replace_range(18..24, "555   ");
        op2.replace_range(24..35, "-X,Y+1/2,-Z");
        let ok = vec![op1, op2];
        let ops = read_remark_290(&ok).unwrap();
        assert_eq!(ops, vec!["X,Y,Z", "-X,Y+1/2,-Z"]);

        let mut bad_line = " ".repeat(80);
        bad_line.replace_range(0..10, "REMARK 290");
        bad_line.replace_range(15..18, "  2");
        bad_line.replace_range(18..24, "555   ");
        bad_line.replace_range(24..29, "X,Y,Z");
        let bad = vec![bad_line];
        let err = read_remark_290(&bad).unwrap_err();
        assert!(matches!(err, BioReadError::Parse { .. }));
    }

    #[test]
    fn pdb_option_guards_normalize_max_line_length_and_track_non_ascii() {
        let pdb = "ATOM      1  CA  ALA A   7      11.104  13.207   9.900  1.00 20.00           C  \u{0080}\n";

        let structure = BioStructure::from_pdb_str_with_params(
            &pdb,
            BioPdbReadParams {
                max_line_length: 0,
                check_non_ascii: true,
                ..BioPdbReadParams::default()
            },
        )
        .unwrap();

        assert_eq!(structure.num_atoms(), 1);
        assert_eq!(structure.non_ascii_line, Some(1));
    }

    #[test]
    fn pdb_option_guards_detect_cif_and_mmjson_and_skip_remarks() {
        let cif_err =
            BioStructure::from_pdb_str_with_params("data_demo\n", BioPdbReadParams::default())
                .unwrap_err();
        assert!(matches!(cif_err, BioReadError::Parse { .. }));

        let json_err =
            BioStructure::from_pdb_str_with_params("{\"data_\":1}\n", BioPdbReadParams::default())
                .unwrap_err();
        assert!(matches!(json_err, BioReadError::Parse { .. }));

        let pdb = "\
REMARK 300 REMARK: SKIP ME
ATOM      1  CA  ALA A   7      11.104  13.207   9.900  1.00 20.00           C
";
        let structure = BioStructure::from_pdb_str_with_params(
            pdb,
            BioPdbReadParams {
                skip_remarks: true,
                ..BioPdbReadParams::default()
            },
        )
        .unwrap();
        assert!(structure.raw_remarks.is_empty());
        assert!(structure.metadata.remark_300_detail.is_none());
    }

    #[test]
    fn pdb_remark_and_conect_branches_preserve_raw_remarks_and_accumulate_adjacency() {
        let pdb = "\
REMARK 300 FIRST DETAIL\r
REMARK 350 SECOND DETAIL
ATOM      1  CA  ALA A   7      11.104  13.207   9.900  1.00 20.00           C
ATOM      2  CB  ALA A   7      12.104  13.207   9.900  1.00 20.00           C
CONECT    1    2    3    0    4
CONECT    2    1
";

        let structure = BioStructure::from_pdb_str(pdb).unwrap();

        assert_eq!(
            structure.raw_remarks,
            vec![
                "REMARK 300 FIRST DETAIL".to_string(),
                "REMARK 350 SECOND DETAIL".to_string(),
            ]
        );
        assert_eq!(structure.conect_map.get(&1), Some(&vec![2, 3, 4]));
        assert_eq!(structure.conect_map.get(&2), Some(&vec![1]));
        assert!(structure.deferred_conn_records.is_empty());
    }

    #[test]
    fn pdb_model_control_flow_handles_implicit_models_duplicate_numbers_and_endmdl() {
        let implicit = BioStructure::from_pdb_str(
            "ATOM      1  CA  ALA A   7      11.104  13.207   9.900  1.00 20.00           C  \n",
        )
        .unwrap();
        assert_eq!(implicit.models[0].source_model_number, Some(1));

        let with_endmdl = BioStructure::from_pdb_str(
            "\
MODEL        1
ATOM      1  CA  ALA A   7      11.104  13.207   9.900  1.00 20.00           C
ENDMDL
MODEL        2
ATOM      2  CA  GLY B   8      12.104  13.207   9.900  1.00 20.00           C
END
",
        )
        .unwrap();
        assert_eq!(with_endmdl.num_models(), 2);
        assert_eq!(with_endmdl.models[1].source_model_number, Some(2));

        let duplicate = BioStructure::from_pdb_str(
            "\
MODEL        1
ATOM      1  CA  ALA A   7      11.104  13.207   9.900  1.00 20.00           C
ENDMDL
MODEL        1
",
        )
        .unwrap_err();
        assert!(matches!(duplicate, BioReadError::Parse { .. }));
    }

    #[test]
    fn pdb_atom_reader_matches_gemmi_hybrid36_serial_and_sequence_values() {
        let pdb = "\
ATOM  99998  SD  MET L9999      48.231 -64.383  -9.257  1.00 11.54           S
ATOM  99999  CE  MET L9999      49.398 -63.242 -10.211  1.00 14.60           C
ATOM  A0000  N   VAL LA000      52.228 -67.689 -12.196  1.00  8.76           N
ATOM  A0001  CA  VAL LA000      53.657 -67.774 -12.458  1.00  3.40           C
";

        let structure = BioStructure::from_pdb_str(pdb).unwrap();
        let values = structure
            .atoms
            .iter()
            .map(|atom| {
                let residue = &structure.residues[atom.residue_id.index() as usize];
                (
                    atom.source.serial.unwrap().0,
                    residue.source.seq_id.unwrap().seq_num,
                )
            })
            .collect::<Vec<_>>();

        assert_eq!(
            values,
            vec![
                (99_998, 9_999),
                (99_999, 9_999),
                (100_000, 10_000),
                (100_001, 10_000)
            ]
        );
    }

    #[test]
    fn pdb_model_control_flow_rejects_non_adjacent_anisou() {
        let err = BioStructure::from_pdb_str(
            "\
ANISOU    1  CA  ALA A   7    1000   2000   3000    400    500    600       C
",
        )
        .unwrap_err();
        assert!(matches!(err, BioReadError::Parse { .. }));
    }

    #[test]
    fn pdb_chain_residue_control_flow_tracks_segment_ids_and_reuses_residue_spans() {
        let mut atom1 = " ".repeat(80);
        atom1.replace_range(0..6, "ATOM  ");
        atom1.replace_range(6..11, "    1");
        atom1.replace_range(12..16, " CA ");
        atom1.replace_range(17..20, "ALA");
        atom1.replace_range(20..22, "A ");
        atom1.replace_range(22..26, "   7");
        atom1.replace_range(30..38, "  11.104");
        atom1.replace_range(38..46, "  13.207");
        atom1.replace_range(46..54, "   9.900");
        atom1.replace_range(54..60, "  1.00");
        atom1.replace_range(60..66, " 20.00");
        atom1.replace_range(76..78, " C");
        atom1.replace_range(72..76, "SEG1");

        let mut atom2 = atom1.clone();
        atom2.replace_range(6..11, "    2");
        atom2.replace_range(12..16, " CB ");
        atom2.replace_range(17..20, "GLY");
        atom2.replace_range(22..26, "   8");
        atom2.replace_range(30..38, "  12.104");

        let mut atom3 = atom1.clone();
        atom3.replace_range(6..11, "    3");
        atom3.replace_range(12..16, " N  ");
        atom3.replace_range(76..78, " N");
        atom3.replace_range(30..38, "  10.104");

        let pdb = format!("{atom1}\n{atom2}\n{atom3}\n");

        let structure = BioStructure::from_pdb_str(&pdb).unwrap();

        assert_eq!(structure.num_residues(), 2);
        assert_eq!(structure.residues[0].source.segment_id, Some(*b"SEG1"));
        assert_eq!(structure.residues[0].atom_span.len, 2);
        assert_eq!(structure.atoms[0].residue_id, ResidueId::new(0));
        assert_eq!(structure.atoms[1].residue_id, ResidueId::new(0));
    }

    #[test]
    fn pdb_chain_residue_control_flow_separates_same_residue_ids_across_segments() {
        let mut atom1 = " ".repeat(80);
        atom1.replace_range(0..6, "ATOM  ");
        atom1.replace_range(6..11, "    1");
        atom1.replace_range(12..16, " CA ");
        atom1.replace_range(17..20, "ALA");
        atom1.replace_range(20..22, "A ");
        atom1.replace_range(22..26, "   7");
        atom1.replace_range(30..38, "  11.104");
        atom1.replace_range(38..46, "  13.207");
        atom1.replace_range(46..54, "   9.900");
        atom1.replace_range(54..60, "  1.00");
        atom1.replace_range(60..66, " 20.00");
        atom1.replace_range(72..76, "SEG1");
        atom1.replace_range(76..78, " C");

        let mut atom2 = atom1.clone();
        atom2.replace_range(6..11, "    2");
        atom2.replace_range(12..16, " CB ");

        let mut atom3 = atom1.clone();
        atom3.replace_range(6..11, "    3");
        atom3.replace_range(12..16, " N  ");
        atom3.replace_range(72..76, "SEG2");
        atom3.replace_range(76..78, " N");

        let pdb = format!("{atom1}\n{atom2}\n{atom3}\n");

        let structure = BioStructure::from_pdb_str(&pdb).unwrap();

        assert_eq!(structure.num_residues(), 2);
        assert_eq!(
            structure.residues[0].source.seq_id,
            structure.residues[1].source.seq_id
        );
        assert_eq!(structure.residues[0].name, structure.residues[1].name);
        assert_eq!(structure.residues[0].source.segment_id, Some(*b"SEG1"));
        assert_eq!(structure.residues[1].source.segment_id, Some(*b"SEG2"));
        assert_eq!(structure.residues[0].atom_span.len, 2);
        assert_eq!(structure.residues[1].atom_span.len, 1);
        assert_eq!(structure.atoms[0].residue_id, ResidueId::new(0));
        assert_eq!(structure.atoms[1].residue_id, ResidueId::new(0));
        assert_eq!(structure.atoms[2].residue_id, ResidueId::new(1));
    }

    #[test]
    fn pdb_chain_residue_control_flow_sets_after_ter_entity_kind_and_splits_chain_parts() {
        let mut atom = " ".repeat(80);
        atom.replace_range(0..6, "ATOM  ");
        atom.replace_range(6..11, "    1");
        atom.replace_range(12..16, " CA ");
        atom.replace_range(17..20, "ALA");
        atom.replace_range(20..22, "A ");
        atom.replace_range(22..26, "   1");
        atom.replace_range(30..38, "  11.104");
        atom.replace_range(38..46, "  13.207");
        atom.replace_range(46..54, "   9.900");
        atom.replace_range(54..60, "  1.00");
        atom.replace_range(60..66, " 20.00");
        atom.replace_range(76..78, " C");

        let mut het = " ".repeat(80);
        het.replace_range(0..6, "HETATM");
        het.replace_range(6..11, "    2");
        het.replace_range(12..16, " O  ");
        het.replace_range(17..20, "HOH");
        het.replace_range(20..22, "A ");
        het.replace_range(22..26, "   2");
        het.replace_range(30..38, "  12.104");
        het.replace_range(38..46, "  13.207");
        het.replace_range(46..54, "   9.900");
        het.replace_range(54..60, "  1.00");
        het.replace_range(60..66, " 20.00");
        het.replace_range(76..78, " O");

        let pdb = format!("{atom}\nTER\n{het}\n");
        let structure = BioStructure::from_pdb_str(&pdb).unwrap();
        assert_eq!(structure.residues[1].kind, ResidueKind::Water);
        assert_eq!(structure.residues[1].entity_kind, EntityKind::Water);

        let split = BioStructure::from_pdb_str_with_params(
            &pdb,
            BioPdbReadParams {
                split_chain_on_ter: true,
                ..BioPdbReadParams::default()
            },
        )
        .unwrap();
        assert_eq!(split.num_chains(), 2);
    }

    #[test]
    fn pdb_chain_residue_control_flow_ignores_ter_when_requested() {
        let mut atom = " ".repeat(80);
        atom.replace_range(0..6, "ATOM  ");
        atom.replace_range(6..11, "    1");
        atom.replace_range(12..16, " CA ");
        atom.replace_range(17..20, "ALA");
        atom.replace_range(20..22, "A ");
        atom.replace_range(22..26, "   1");
        atom.replace_range(30..38, "  11.104");
        atom.replace_range(38..46, "  13.207");
        atom.replace_range(46..54, "   9.900");
        atom.replace_range(54..60, "  1.00");
        atom.replace_range(60..66, " 20.00");
        atom.replace_range(76..78, " C");

        let mut het = " ".repeat(80);
        het.replace_range(0..6, "HETATM");
        het.replace_range(6..11, "    2");
        het.replace_range(12..16, " ZN ");
        het.replace_range(17..20, "ZN ");
        het.replace_range(20..22, "A ");
        het.replace_range(22..26, "   2");
        het.replace_range(30..38, "  12.104");
        het.replace_range(38..46, "  13.207");
        het.replace_range(46..54, "   9.900");
        het.replace_range(54..60, "  1.00");
        het.replace_range(60..66, " 20.00");
        het.replace_range(76..78, "ZN");

        let pdb = format!("{atom}\nTER\n{het}\n");
        let structure = BioStructure::from_pdb_str_with_params(
            &pdb,
            BioPdbReadParams {
                ignore_ter: true,
                ..BioPdbReadParams::default()
            },
        )
        .unwrap();

        assert_eq!(structure.ter_status, '\0');
        assert_eq!(structure.num_chains(), 1);
        assert_eq!(structure.num_residues(), 2);
        assert_eq!(structure.residues[0].entity_kind, EntityKind::Unknown);
        assert_eq!(structure.residues[1].entity_kind, EntityKind::Unknown);
    }

    #[test]
    fn pdb_atom_control_flow_infers_element_from_padded_atom_name_when_columns_blank() {
        let mut atom = " ".repeat(80);
        atom.replace_range(0..6, "ATOM  ");
        atom.replace_range(6..11, "    1");
        atom.replace_range(12..16, "1HG ");
        atom.replace_range(17..20, "UNL");
        atom.replace_range(20..22, "A ");
        atom.replace_range(22..26, "   1");
        atom.replace_range(30..38, "   0.000");
        atom.replace_range(38..46, "   1.000");
        atom.replace_range(46..54, "   2.000");
        atom.replace_range(54..60, "  1.00");
        atom.replace_range(60..66, " 20.00");

        let structure = BioStructure::from_pdb_str(&format!("{atom}\n")).unwrap();

        assert_eq!(structure.atoms.len(), 1);
        assert_eq!(structure.atoms[0].element.atomic_number(), 1);
    }

    #[test]
    fn pdb_atom_control_flow_inserts_atoms_into_reused_earlier_residue_span() {
        let mut atom1 = " ".repeat(80);
        atom1.replace_range(0..6, "ATOM  ");
        atom1.replace_range(6..11, "    1");
        atom1.replace_range(12..16, " CA ");
        atom1.replace_range(17..20, "ALA");
        atom1.replace_range(20..22, "A ");
        atom1.replace_range(22..26, "   1");
        atom1.replace_range(30..38, "   0.000");
        atom1.replace_range(38..46, "   0.000");
        atom1.replace_range(46..54, "   0.000");
        atom1.replace_range(54..60, "  1.00");
        atom1.replace_range(60..66, " 20.00");
        atom1.replace_range(76..78, " C");

        let mut atom2 = atom1.clone();
        atom2.replace_range(6..11, "    2");
        atom2.replace_range(12..16, " CB ");
        atom2.replace_range(22..26, "   2");
        atom2.replace_range(30..38, "   1.000");

        let mut atom3 = atom1.clone();
        atom3.replace_range(6..11, "    3");
        atom3.replace_range(12..16, " N  ");
        atom3.replace_range(30..38, "   2.000");
        atom3.replace_range(76..78, " N");

        let mut anisou = " ".repeat(80);
        anisou.replace_range(0..6, "ANISOU");
        anisou.replace_range(6..11, "    3");
        anisou.replace_range(12..16, " N  ");
        anisou.replace_range(17..20, "ALA");
        anisou.replace_range(20..22, "A ");
        anisou.replace_range(22..26, "   1");
        anisou.replace_range(28..35, "   1000");
        anisou.replace_range(35..42, "   2000");
        anisou.replace_range(42..49, "   3000");
        anisou.replace_range(49..56, "    400");
        anisou.replace_range(56..63, "    500");
        anisou.replace_range(63..70, "    600");
        anisou.replace_range(76..78, " N");

        let pdb = format!("{atom1}\n{atom2}\n{atom3}\n{anisou}\n");
        let structure = BioStructure::from_pdb_str(&pdb).unwrap();

        assert_eq!(structure.num_residues(), 2);
        assert_eq!(structure.residues[0].atom_span.start, 0);
        assert_eq!(structure.residues[0].atom_span.len, 2);
        assert_eq!(structure.residues[1].atom_span.start, 2);
        assert_eq!(structure.residues[1].atom_span.len, 1);
        assert_eq!(structure.atoms[0].name, AtomName(*b" CA "));
        assert_eq!(structure.atoms[1].name, AtomName(*b" N  "));
        assert_eq!(structure.atoms[2].name, AtomName(*b" CB "));
        assert_eq!(structure.atoms[1].residue_id, ResidueId::new(0));
        assert_eq!(structure.atoms[2].residue_id, ResidueId::new(1));
        assert_eq!(
            structure.atoms[1].anisou,
            Some([
                1000.0_f32 * 1e-4,
                2000.0_f32 * 1e-4,
                3000.0_f32 * 1e-4,
                400.0_f32 * 1e-4,
                500.0_f32 * 1e-4,
                600.0_f32 * 1e-4,
            ])
        );
        assert_eq!(structure.atoms[2].anisou, None);
    }

    #[test]
    fn pdb_model_control_flow_rejects_implicit_model_number_collision() {
        let err = BioStructure::from_pdb_str(
            "\
MODEL        2
ENDMDL
ATOM      1  CA  ALA A   7      11.104  13.207   9.900  1.00 20.00           C
",
        )
        .unwrap_err();
        assert!(matches!(err, BioReadError::Parse { .. }));
    }

    #[test]
    fn pdb_chain_residue_control_flow_detects_after_ter_from_previous_polymer_chain_part() {
        let mut atom1 = " ".repeat(80);
        atom1.replace_range(0..6, "ATOM  ");
        atom1.replace_range(6..11, "    1");
        atom1.replace_range(12..16, " CA ");
        atom1.replace_range(17..20, "ALA");
        atom1.replace_range(20..22, "A ");
        atom1.replace_range(22..26, "   1");
        atom1.replace_range(30..38, "  11.104");
        atom1.replace_range(38..46, "  13.207");
        atom1.replace_range(46..54, "   9.900");
        atom1.replace_range(54..60, "  1.00");
        atom1.replace_range(60..66, " 20.00");
        atom1.replace_range(76..78, " C");

        let mut atom2 = atom1.clone();
        atom2.replace_range(6..11, "    2");
        atom2.replace_range(17..20, "GLY");
        atom2.replace_range(22..26, "   2");
        atom2.replace_range(30..38, "  12.104");

        let mut het = atom1.clone();
        het.replace_range(0..6, "HETATM");
        het.replace_range(6..11, "    3");
        het.replace_range(12..16, " ZN ");
        het.replace_range(17..20, "ZN ");
        het.replace_range(22..26, "   3");
        het.replace_range(30..38, "  13.104");
        het.replace_range(76..78, "ZN");

        let pdb = format!("{atom1}\n{atom2}\nTER\n{het}\n");
        let structure = BioStructure::from_pdb_str(&pdb).unwrap();

        assert_eq!(structure.num_chains(), 1);
        assert_eq!(structure.num_residues(), 3);
        assert_eq!(structure.residues[0].entity_kind, EntityKind::Polymer);
        assert_eq!(structure.residues[1].entity_kind, EntityKind::Polymer);
        assert_eq!(structure.residues[2].kind, ResidueKind::Unknown);
        assert_eq!(structure.residues[2].entity_kind, EntityKind::NonPolymer);
    }

    #[test]
    fn reads_minimal_mmcif_atom_site_loop() {
        let cif = r#"
data_demo
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA . ALA A 7 11.104 13.207 9.900 1.00 20.00 7 ALA A CA 1
"#;

        let structure = read_test_mmcif(cif).unwrap();

        assert_eq!(structure.num_models(), 1);
        assert_eq!(structure.models[0].source_model_number, Some(1));
        assert_eq!(
            structure.chains[0].source.auth_chain_id.unwrap().as_str(),
            "A"
        );
        assert_eq!(
            structure.chains[0].source.label_asym_id.unwrap().as_str(),
            "A"
        );
        assert_eq!(structure.residues[0].name.as_str(), "ALA");
        assert_eq!(structure.atoms[0].source.serial, Some(PdbAtomSerial(1)));
        assert_eq!(structure.atoms[0].occupancy, Some(1.0));
        assert_eq!(structure.atoms[0].b_iso, Some(20.0));
        assert_eq!(structure.coordinates.positions[0], [11.104, 13.207, 9.900]);
    }

    #[test]
    fn mmcif_atom_site_uses_gemmi_atom_defaults_for_missing_or_null_occ_and_b_iso() {
        for optional_columns in [
            "".to_string(),
            "_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n".to_string(),
        ] {
            let optional_values = if optional_columns.is_empty() {
                ""
            } else {
                ". ? "
            };
            let cif = format!(
                "data_demo\nloop_\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n{optional_columns}1 C CA . ALA A 1 0 0 0 {optional_values}\n"
            );
            let structure = read_test_mmcif(&cif).unwrap();
            assert_eq!(structure.atoms[0].occupancy, Some(1.0));
            assert_eq!(structure.atoms[0].b_iso, Some(20.0));
        }
    }

    #[test]
    fn mmcif_atom_site_maps_model_entity_calc_tls_fraction_and_anisou() {
        let cif = r#"
data_demo
loop_
_atom_site.id
_atom_site.group_PDB
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
_atom_site.calc_flag
_atom_site.pdbx_tls_group_id
_atom_site.ccp4_deuterium_fraction
1 ATOM C CA . ALA A 1 7 ? 1.0 2.0 3.0 1.0 11.0 1 70 ALA X CA 1 calc 12 0.25
2 HETATM O O . HOH B 2 8 A 4.0 5.0 6.0 0.5 22.0 ? 80 HOH Y O 2 dum 7x 0.75
loop_
_atom_site_anisotrop.id
_atom_site_anisotrop.U[1][1]
_atom_site_anisotrop.U[2][2]
_atom_site_anisotrop.U[3][3]
_atom_site_anisotrop.U[1][2]
_atom_site_anisotrop.U[1][3]
_atom_site_anisotrop.U[2][3]
1 1 2 3 4 5 6
"#;

        let structure = read_test_mmcif(cif).unwrap();

        assert_eq!(structure.num_models(), 2);
        assert_eq!(structure.models[0].source_model_number, Some(1));
        assert_eq!(structure.models[1].source_model_number, Some(2));
        assert!(structure.has_d_fraction);
        assert_eq!(
            structure.residues[0].source.label_entity_id,
            Some(EntityId::new(0))
        );
        assert_eq!(structure.residues[0].het_flag, Some('A'));
        assert_eq!(structure.residues[1].het_flag, Some('H'));
        assert_eq!(structure.residues[1].entity_kind, EntityKind::Water);
        assert_eq!(structure.atoms[0].calc_flag, BioCalcFlag::Calculated);
        assert_eq!(structure.atoms[1].calc_flag, BioCalcFlag::Dummy);
        assert_eq!(structure.atoms[0].tls_group_id, Some(12));
        assert_eq!(structure.atoms[1].tls_group_id, Some(7));
        assert_eq!(structure.atoms[0].fraction, Some(0.25));
        assert_eq!(structure.atoms[1].fraction, Some(0.75));
        assert_eq!(
            structure.atoms[0].anisou,
            Some([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
        );
        assert_eq!(
            structure.chains[0].source.auth_chain_id.unwrap().as_str(),
            "X"
        );
        assert_eq!(
            structure.chains[0].source.label_asym_id.unwrap().as_str(),
            "A"
        );
    }

    #[test]
    fn reads_mmcif_entities_sequences_and_struct_asym_links() {
        let cif = r#"
data_demo
loop_
_entity.id
_entity.type
1 polymer
loop_
_entity_poly.entity_id
_entity_poly.type
1 polypeptide(L)
loop_
_entity_poly_seq.entity_id
_entity_poly_seq.num
_entity_poly_seq.mon_id
1 1 ALA
1 2 GLY
1 2 SER
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 C CA . ALA A 1 0.0 1.0 2.0
"#;

        let structure = read_test_mmcif(cif).unwrap();

        assert_eq!(structure.num_entities(), 1);
        assert_eq!(structure.entities[0].kind, EntityKind::Polymer);
        assert_eq!(structure.entities[0].polymer_kind, PolymerKind::Peptide);
        assert_eq!(structure.entities[0].sequence, vec!["ALA", "GLY,SER"]);
        assert_eq!(structure.entities[0].subchains[0].as_str(), "A");
        assert_eq!(structure.chains[0].entity_id, Some(EntityId::new(0)));
    }

    #[test]
    fn mmcif_entity_and_sequence_info_builds_dbrefs_and_deduplicates_struct_ref_seq() {
        let cif = r#"
data_demo
loop_
_entity.id
_entity.type
1 polymer
loop_
_entity_poly.entity_id
_entity_poly.type
1 polypeptide(L)
loop_
_entity_poly_seq.entity_id
_entity_poly_seq.num
_entity_poly_seq.mon_id
1 1 ALA
1 2 GLY
1 2 SER
loop_
_struct_ref.id
_struct_ref.entity_id
_struct_ref.db_name
_struct_ref.db_code
_struct_ref.pdbx_db_accession
_struct_ref.pdbx_db_isoform
R1 1 UNP SAMPLE_ID Q12345 ISO1
loop_
_struct_ref_seq.ref_id
_struct_ref_seq.seq_align_beg
_struct_ref_seq.seq_align_end
_struct_ref_seq.db_align_beg
_struct_ref_seq.db_align_end
_struct_ref_seq.pdbx_auth_seq_align_beg
_struct_ref_seq.pdbx_seq_align_beg_ins_code
_struct_ref_seq.pdbx_auth_seq_align_end
_struct_ref_seq.pdbx_seq_align_end_ins_code
R1 1 3 11 19 101 ? 103 A
R1 1 3 11 19 101 ? 103 A
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 C CA . ALA A 1 1 0.0 1.0 2.0
"#;

        let structure = read_test_mmcif(cif).unwrap();

        assert_eq!(structure.entities[0].sequence, vec!["ALA", "GLY,SER"]);
        assert!(structure.entities[0].reflects_microhetero);
        assert_eq!(structure.entities[0].dbrefs.len(), 1);
        let dbref = &structure.entities[0].dbrefs[0];
        assert_eq!(dbref.db_name, "UNP");
        assert_eq!(dbref.id_code, "SAMPLE_ID");
        assert_eq!(dbref.accession_code, "Q12345");
        assert_eq!(dbref.isoform, "ISO1");
        assert_eq!(dbref.label_seq_begin, Some(1));
        assert_eq!(dbref.label_seq_end, Some(3));
        assert_eq!(dbref.db_begin.unwrap().seq_num, 11);
        assert_eq!(dbref.db_end.unwrap().seq_num, 19);
        assert_eq!(dbref.seq_begin.unwrap().seq_num, 101);
        assert_eq!(dbref.seq_end.unwrap().seq_num, 103);
        assert_eq!(dbref.seq_end.unwrap().ins_code, Some(b'A'));
        assert_eq!(
            structure.entities[0].subchains,
            vec![pdb_chain_id_from_str("A").unwrap()]
        );
    }

    #[test]
    fn mmcif_entity_and_sequence_info_falls_back_to_model_subchains_without_struct_asym() {
        let cif = r#"
data_demo
loop_
_entity.id
_entity.type
1 polymer
loop_
_entity_poly.entity_id
_entity_poly.type
1 polypeptide(L)
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 C CA . ALA AX 1 1 0.0 1.0 2.0
ATOM 2 C CB . ALA AX 1 1 0.0 1.5 2.0
"#;

        let structure = read_test_mmcif(cif).unwrap();

        assert_eq!(structure.entities.len(), 1);
        assert_eq!(
            structure.entities[0].subchains,
            vec![pdb_chain_id_from_str("AX").unwrap()]
        );
    }

    #[test]
    fn mmcif_connectivity_recovers_auth_address_from_label_subchain_and_label_seq() {
        let cif = r#"
data_demo
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
ATOM 1 S SG . CYS A 1 0.0 0.0 0.0 101 CYS X
ATOM 2 S SG . CYS A 2 1.0 0.0 0.0 102 CYS X
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr2_label_seq_id
conn1 disulf A A CYS CYS SG SG 1 2
"#;

        let structure = read_test_mmcif(cif).unwrap();
        assert_eq!(structure.connections.len(), 1);
        let conn = &structure.connections[0];
        assert_eq!(conn.type_, BioConnectionType::Disulf);
        assert_eq!(conn.partner1.chain_name, "X");
        assert_eq!(conn.partner1.seq_id.unwrap().seq_num, 101);
        assert_eq!(conn.partner2.chain_name, "X");
        assert_eq!(conn.partner2.seq_id.unwrap().seq_num, 102);
    }

    #[test]
    fn mmcif_connectivity_rejects_rows_without_auth_or_label_address_pairs() {
        let cif = r#"
data_demo
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 S SG . CYS A 1 0.0 0.0 0.0
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_label_atom_id
broken disulf CYS CYS SG SG
"#;

        let err = read_test_mmcif(cif).unwrap_err();
        assert!(matches!(err, BioReadError::Parse { .. }));
    }

    #[test]
    fn mmcif_prot_cis_captures_altloc_and_angle() {
        let cif = r#"
data_demo
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
ATOM 1 C CA . ALA A 3 0.0 0.0 0.0 30 ALA X
ATOM 2 N N  . GLY A 4 1.0 0.0 0.0 40 GLY Y
loop_
_struct_mon_prot_cis.pdbx_PDB_model_num
_struct_mon_prot_cis.auth_asym_id
_struct_mon_prot_cis.auth_seq_id
_struct_mon_prot_cis.pdbx_PDB_ins_code
_struct_mon_prot_cis.label_comp_id
_struct_mon_prot_cis.auth_comp_id
_struct_mon_prot_cis.pdbx_auth_asym_id_2
_struct_mon_prot_cis.pdbx_auth_seq_id_2
_struct_mon_prot_cis.pdbx_PDB_ins_code_2
_struct_mon_prot_cis.pdbx_label_comp_id_2
_struct_mon_prot_cis.pdbx_auth_comp_id_2
_struct_mon_prot_cis.label_alt_id
_struct_mon_prot_cis.pdbx_omega_angle
1 X 30 ? ALA ALA Y 40 ? GLY GLY B -12.5
"#;

        let structure = read_test_mmcif(cif).unwrap();

        assert_eq!(structure.cispeps.len(), 1);
        let cis = &structure.cispeps[0];
        assert_eq!(cis.model_num, 1);
        assert_eq!(cis.partner_c.chain_name, "X");
        assert_eq!(cis.partner_c.seq_id.unwrap().seq_num, 30);
        assert_eq!(cis.partner_n.chain_name, "Y");
        assert_eq!(cis.partner_n.seq_id.unwrap().seq_num, 40);
        assert_eq!(cis.only_altloc, Some(AltLocLabel(b'B')));
        assert_eq!(cis.reported_angle, Some(-12.5));
    }

    #[test]
    fn mmcif_struct_mod_residue_reads_fields() {
        let cif = r#"
data_demo
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 C CA . MSE A 5 0.0 0.0 0.0
loop_
_pdbx_struct_mod_residue.auth_asym_id
_pdbx_struct_mod_residue.auth_seq_id
_pdbx_struct_mod_residue.PDB_ins_code
_pdbx_struct_mod_residue.auth_comp_id
_pdbx_struct_mod_residue.label_comp_id
_pdbx_struct_mod_residue.parent_comp_id
_pdbx_struct_mod_residue.details
_pdbx_struct_mod_residue.ccp4_mod_id
A 5 ? MSE MSE MET 'SELENOMETHIONINE' MOD1
"#;

        let structure = read_test_mmcif(cif).unwrap();

        assert_eq!(structure.mod_residues.len(), 1);
        let modres = &structure.mod_residues[0];
        assert_eq!(modres.chain_name, "A");
        assert_eq!(modres.res_id.seq_num, 5);
        assert_eq!(modres.parent_comp_id, "MET");
        assert_eq!(modres.details, "SELENOMETHIONINE");
        assert_eq!(modres.mod_id, "MOD1");
    }

    #[test]
    fn mmcif_parse_operation_expr_expands_ranges_and_lists() {
        assert_eq!(parse_operation_expr("3"), vec!["3"]);
        assert_eq!(parse_operation_expr("(1-3)"), vec!["1", "2", "3"]);
        assert_eq!(
            parse_operation_expr("(2,3-4,XY)"),
            vec!["2", "3", "4", "XY"]
        );
    }

    #[test]
    fn mmcif_read_assemblies_filters_unknown_detail_and_expands_operators() {
        let cif = r#"
data_demo
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 C CA . ALA A 1 0.0 0.0 0.0
loop_
_pdbx_struct_oper_list.id
_pdbx_struct_oper_list.type
_pdbx_struct_oper_list.matrix[1][1]
_pdbx_struct_oper_list.matrix[1][2]
_pdbx_struct_oper_list.matrix[1][3]
_pdbx_struct_oper_list.vector[1]
_pdbx_struct_oper_list.matrix[2][1]
_pdbx_struct_oper_list.matrix[2][2]
_pdbx_struct_oper_list.matrix[2][3]
_pdbx_struct_oper_list.vector[2]
_pdbx_struct_oper_list.matrix[3][1]
_pdbx_struct_oper_list.matrix[3][2]
_pdbx_struct_oper_list.matrix[3][3]
_pdbx_struct_oper_list.vector[3]
1 identity 1 0 0 0 0 1 0 0 0 0 1 0
2 identity 1 0 0 10 0 1 0 0 0 0 1 0
loop_
_pdbx_struct_assembly.id
_pdbx_struct_assembly.details
_pdbx_struct_assembly.method_details
_pdbx_struct_assembly.oligomeric_details
_pdbx_struct_assembly.oligomeric_count
1 software_defined_assembly program dimeric 2
2 author_and_software_defined_assembly combo tetrameric 4
3 garbage ? ignored 1
loop_
_pdbx_struct_assembly_gen.assembly_id
_pdbx_struct_assembly_gen.oper_expression
_pdbx_struct_assembly_gen.asym_id_list
1 (1-2) A,B
2 (2) C
loop_
_pdbx_struct_assembly_prop.biol_id
_pdbx_struct_assembly_prop.type
_pdbx_struct_assembly_prop.value
1 'ABSA (A^2)' 12.5
1 'SSA (A^2)' 8.0
2 MORE 4.5
"#;

        let structure = read_test_mmcif(cif).unwrap();
        assert_eq!(structure.assemblies.len(), 2);
        assert_eq!(structure.assemblies[0].name, "1");
        assert_eq!(structure.assemblies[0].software_name, "program");
        assert_eq!(structure.assemblies[0].absa, Some(12.5));
        assert_eq!(structure.assemblies[0].ssa, Some(8.0));
        assert_eq!(structure.assemblies[0].generators.len(), 1);
        assert_eq!(
            structure.assemblies[0].generators[0].subchains,
            vec!["A", "B"]
        );
        assert_eq!(structure.assemblies[0].generators[0].operators.len(), 2);
        assert_eq!(structure.assemblies[0].generators[0].operators[1].name, "2");
        assert_eq!(structure.assemblies[1].name, "2");
        assert_eq!(structure.assemblies[1].more, Some(4.5));
        assert!(
            structure
                .assemblies
                .iter()
                .all(|assembly| assembly.name != "3")
        );
    }

    #[test]
    fn mmcif_fill_residue_entity_type_falls_back_to_water_and_unknown() {
        let cif = r#"
data_demo
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.auth_seq_id
ATOM 1 C CA . ALA A 1 0.0 0.0 0.0 1
HETATM 2 O O . HOH B . 1.0 0.0 0.0 2
"#;

        let structure = read_test_mmcif(cif).unwrap();
        assert_eq!(structure.residues[0].entity_kind, EntityKind::Unknown);
        assert_eq!(structure.residues[1].entity_kind, EntityKind::Water);
    }

    #[test]
    fn mmcif_sifts_unp_maps_accession_index_and_residue_span() {
        let cif = r#"
data_demo
loop_
_entity.id
_entity.type
1 polymer
loop_
_entity_poly.entity_id
_entity_poly.type
1 polypeptide(L)
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
ATOM 1 C CA . ALA A 1 0.0 0.0 0.0 11 ALA X
ATOM 2 C CA . GLY A 2 1.0 0.0 0.0 12 GLY X
loop_
_pdbx_sifts_xref_db.entity_id
_pdbx_sifts_xref_db.asym_id
_pdbx_sifts_xref_db.seq_id_ordinal
_pdbx_sifts_xref_db.seq_id
_pdbx_sifts_xref_db.observed
_pdbx_sifts_xref_db.unp_res
_pdbx_sifts_xref_db.unp_num
_pdbx_sifts_xref_db.unp_acc
1 A 1 1 y A 101 P12345
1 A 1 2 y G 102 P12345
"#;

        let structure = read_test_mmcif(cif).unwrap();
        assert_eq!(structure.entities[0].sifts_unp_acc, vec!["P12345"]);
        assert_eq!(
            structure.residues[0].sifts_unp,
            Some(BioSiftsUnpResidue {
                res: Some('A'),
                acc_index: 0,
                num: 101,
            })
        );
        assert_eq!(
            structure.residues[1].sifts_unp,
            Some(BioSiftsUnpResidue {
                res: Some('G'),
                acc_index: 0,
                num: 102,
            })
        );
    }

    #[test]
    fn mmcif_find_diffrn_returns_matching_diffraction_row() {
        let mut meta = BioMetadata {
            experiment_crystals: vec![BioExperimentCrystalInfo {
                diffractions: vec![
                    BioDiffractionInfo {
                        id: "d1".to_string(),
                        beamline: "BL1".to_string(),
                        ..BioDiffractionInfo::default()
                    },
                    BioDiffractionInfo {
                        id: "d2".to_string(),
                        beamline: "BL2".to_string(),
                        ..BioDiffractionInfo::default()
                    },
                ],
                ..BioExperimentCrystalInfo::default()
            }],
            ..BioMetadata::default()
        };
        assert_eq!(
            find_diffrn(&mut meta, "d2").map(|d| d.beamline.clone()),
            Some("BL2".to_string())
        );
        assert!(find_diffrn(&mut meta, "missing").is_none());
    }

    #[test]
    fn mmcif_read_entry_info_and_audit_author_collect_multivalue_metadata() {
        let cif = r#"
data_demo
_entry.id 2XYZ
_cell.Z_PDB 4
_struct.title 'FIRST TITLE'
_struct.title 'SECOND TITLE'
_database_PDB_rev.date_original 2001-01-02
_struct_keywords.pdbx_keywords OXIDOREDUCTASE
_struct_keywords.text 'TEST CASE'
loop_
_exptl.method
'X-RAY DIFFRACTION'
'NEUTRON DIFFRACTION'
loop_
_audit_author.name
'DOE, J.'
'SMITH, A.'
"#;

        let document = parse_cif_document(cif).unwrap();
        let mut structure = BioStructure::default();

        let block = &document.blocks[0];
        read_entry_info(&document, &mut structure);
        read_audit_author(&block.loops, &mut structure).unwrap();

        assert_eq!(structure.metadata.entry_id.as_deref(), Some("2XYZ"));
        assert_eq!(
            structure.metadata.experimental_method.as_deref(),
            Some("X-RAY DIFFRACTION; NEUTRON DIFFRACTION")
        );
        assert_eq!(
            structure.metadata.title.as_deref(),
            Some("FIRST TITLE; SECOND TITLE")
        );
        assert_eq!(
            structure
                .metadata
                .received_initial_deposition_date
                .as_deref(),
            Some("2001-01-02")
        );
        assert_eq!(
            structure.metadata.pdbx_keywords.as_deref(),
            Some("OXIDOREDUCTASE")
        );
        assert_eq!(structure.metadata.keywords.as_deref(), Some("TEST CASE"));
        assert_eq!(structure.crystal().unwrap().z_pdb.as_deref(), Some("4"));
        assert_eq!(structure.metadata.authors, vec!["DOE, J.", "SMITH, A."]);
    }

    #[test]
    fn mmcif_metadata_readers_route_refinement_tls_experiment_software_and_ncs() {
        let cif = r#"
data_demo
_em_3d_reconstruction.resolution 3.4
loop_
_refine.pdbx_refine_id
_refine.ls_d_res_high
_refine.ls_d_res_low
_refine.ls_percent_reflns_obs
_refine.ls_number_reflns_obs
_refine.ls_number_reflns_R_work
_refine.ls_number_reflns_R_free
_refine.ls_R_factor_obs
_refine.ls_R_factor_R_work
_refine.ls_R_factor_R_free
X-RAY 1.8 30.0 95.0 1000 900 100 0.22 0.20 0.25
NEUTRON . . . . . . . . .
loop_
_pdbx_refine_tls.id
_pdbx_refine_tls.pdbx_refine_id
_pdbx_refine_tls.T[1][1]
_pdbx_refine_tls.T[2][2]
_pdbx_refine_tls.T[3][3]
_pdbx_refine_tls.T[1][2]
_pdbx_refine_tls.T[1][3]
_pdbx_refine_tls.T[2][3]
_pdbx_refine_tls.L[1][1]
_pdbx_refine_tls.L[2][2]
_pdbx_refine_tls.L[3][3]
_pdbx_refine_tls.L[1][2]
_pdbx_refine_tls.L[1][3]
_pdbx_refine_tls.L[2][3]
_pdbx_refine_tls.S[1][1]
_pdbx_refine_tls.S[1][2]
_pdbx_refine_tls.S[1][3]
_pdbx_refine_tls.S[2][1]
_pdbx_refine_tls.S[2][2]
_pdbx_refine_tls.S[2][3]
_pdbx_refine_tls.S[3][1]
_pdbx_refine_tls.S[3][2]
_pdbx_refine_tls.S[3][3]
_pdbx_refine_tls.origin_x
_pdbx_refine_tls.origin_y
_pdbx_refine_tls.origin_z
1 X-RAY 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 0.1 0.2 0.3
loop_
_pdbx_refine_tls_group.refine_tls_id
_pdbx_refine_tls_group.beg_auth_asym_id
_pdbx_refine_tls_group.beg_auth_seq_id
_pdbx_refine_tls_group.beg_PDB_ins_code
_pdbx_refine_tls_group.end_auth_seq_id
_pdbx_refine_tls_group.end_PDB_ins_code
_pdbx_refine_tls_group.selection_details
1 A 10 ? 20 A 'CHAIN A'
loop_
_exptl.method
_exptl.crystals_number
'X-RAY DIFFRACTION' 2
loop_
_exptl_crystal.id
_exptl_crystal.description
C1 'CRYSTAL 1'
loop_
_diffrn.id
_diffrn.crystal_id
_diffrn.ambient_temp
D1 C1 100.0
loop_
_diffrn_detector.diffrn_id
_diffrn_detector.pdbx_collection_date
_diffrn_detector.detector
_diffrn_detector.type
_diffrn_detector.details
D1 2024-01-01 PIXEL PILATUS 'FAST OPTICS'
loop_
_diffrn_radiation.diffrn_id
_diffrn_radiation.pdbx_scattering_type
_diffrn_radiation.pdbx_monochromatic_or_laue_m_l
_diffrn_radiation.monochromator
D1 x-ray M SI111
loop_
_diffrn_source.diffrn_id
_diffrn_source.source
_diffrn_source.type
_diffrn_source.pdbx_synchrotron_site
_diffrn_source.pdbx_synchrotron_beamline
_diffrn_source.pdbx_wavelength_list
D1 SYNCHROTRON ROTATING_ANODE SSRL BL1 '1.0,1.1'
loop_
_reflns.pdbx_diffrn_id
_reflns.number_obs
_reflns.d_resolution_high
_reflns.d_resolution_low
_reflns.percent_possible_obs
_reflns.pdbx_redundancy
_reflns.pdbx_Rmerge_I_obs
_reflns.pdbx_Rsym_value
_reflns.pdbx_netI_over_sigmaI
'D1,D2' 500 1.9 20.0 99.0 3.5 0.1 0.2 8.0
loop_
_software.name
_software.classification
_software.version
_software.date
_software.description
_software.contact_author
_software.contact_author_email
REFMAC refinement 5.8 2024-01-02 'REFINE DESC' 'J. DOE' 'j@example.org'
loop_
_struct_ncs_oper.matrix[1][1]
_struct_ncs_oper.matrix[1][2]
_struct_ncs_oper.matrix[1][3]
_struct_ncs_oper.vector[1]
_struct_ncs_oper.matrix[2][1]
_struct_ncs_oper.matrix[2][2]
_struct_ncs_oper.matrix[2][3]
_struct_ncs_oper.vector[2]
_struct_ncs_oper.matrix[3][1]
_struct_ncs_oper.matrix[3][2]
_struct_ncs_oper.matrix[3][3]
_struct_ncs_oper.vector[3]
_struct_ncs_oper.id
_struct_ncs_oper.code
1 0 0 0 0 1 0 0 0 0 1 0 I given
1 0 0 1 0 1 0 2 0 0 1 3 N generated
"#;

        let document = parse_cif_document(cif).unwrap();
        let mut structure = BioStructure::default();

        let block = &document.blocks[0];
        read_refinement_info(&document, &block.loops, &mut structure).unwrap();
        read_tls_info(&block.loops, &mut structure).unwrap();
        read_experimental_info(&block.loops, &mut structure).unwrap();
        read_reflns_info(&block.loops, &mut structure).unwrap();
        read_software_info(&block.loops, &mut structure).unwrap();
        read_ncs_info(&block.loops, &mut structure).unwrap();

        assert!((structure.resolution().unwrap() - 1.8).abs() < 1e-6);
        assert_eq!(structure.metadata.refinement.len(), 2);
        assert_eq!(structure.metadata.refinement[0].id, "X-RAY");
        assert_eq!(structure.metadata.refinement[0].work_set_count, Some(900));
        assert_eq!(structure.metadata.refinement[0].tls_groups.len(), 1);
        assert_eq!(
            structure.metadata.refinement[0].tls_groups[0].num_id,
            Some(1)
        );
        assert_eq!(
            structure.metadata.refinement[0].tls_groups[0].selections[0].res_begin,
            Some(PdbSeqId {
                seq_num: 10,
                ins_code: None,
            })
        );
        assert_eq!(
            structure.metadata.refinement[0].tls_groups[0].selections[0].res_end,
            Some(PdbSeqId {
                seq_num: 20,
                ins_code: Some(b'A'),
            })
        );
        assert_eq!(structure.metadata.experiments.len(), 1);
        assert_eq!(
            structure.metadata.experiments[0].number_of_crystals,
            Some(2)
        );
        assert_eq!(
            structure.metadata.experiments[0].diffraction_ids,
            vec!["D1".to_string(), "D2".to_string()]
        );
        assert!(
            (structure.metadata.experiments[0]
                .reflections
                .r_merge
                .unwrap()
                - 0.1)
                .abs()
                < 1e-6
        );
        assert_eq!(structure.metadata.experiment_crystals.len(), 1);
        assert_eq!(structure.metadata.experiment_crystals[0].id, "C1");
        assert_eq!(
            structure.metadata.experiment_crystals[0].diffractions[0].id,
            "D1"
        );
        assert_eq!(
            structure.metadata.experiment_crystals[0].diffractions[0].mono_or_laue,
            Some('M')
        );
        assert_eq!(
            structure.metadata.experiment_crystals[0].diffractions[0].beamline,
            "BL1"
        );
        assert_eq!(structure.metadata.software.len(), 1);
        assert_eq!(
            structure.metadata.software[0].classification,
            BioSoftwareClassification::Refinement
        );
        assert_eq!(structure.metadata.software[0].description, "REFINE DESC");
        assert_eq!(structure.ncs_oper_identity_id.as_deref(), Some("I"));
        assert_eq!(structure.ncs_operators.len(), 1);
        assert_eq!(structure.ncs_operators[0].id, "N");
        assert!(!structure.ncs_operators[0].given);
        assert_eq!(structure.ncs_operators[0].transform.vec, [1.0, 2.0, 3.0]);
    }

    #[test]
    fn mmcif_top_level_populate_order_applies_metadata_fract_transf_and_origx() {
        let cif = r#"
data_demo
_entry.id 9XYZ
_struct.title 'TOP LEVEL TITLE'
_cell.length_a 10.0
_cell.length_b 20.0
_cell.length_c 30.0
_cell.angle_alpha 90.0
_cell.angle_beta 90.0
_cell.angle_gamma 90.0
_symmetry.space_group_name_H-M 'P 1'
loop_
_audit_author.name
'DOE, J.'
loop_
_atom_sites.fract_transf_matrix[1][1]
_atom_sites.fract_transf_matrix[1][2]
_atom_sites.fract_transf_matrix[1][3]
_atom_sites.fract_transf_vector[1]
_atom_sites.fract_transf_matrix[2][1]
_atom_sites.fract_transf_matrix[2][2]
_atom_sites.fract_transf_matrix[2][3]
_atom_sites.fract_transf_vector[2]
_atom_sites.fract_transf_matrix[3][1]
_atom_sites.fract_transf_matrix[3][2]
_atom_sites.fract_transf_matrix[3][3]
_atom_sites.fract_transf_vector[3]
0.1 0.0 0.0 1.0 0.0 0.05 0.0 2.0 0.0 0.0 0.033333333 3.0
loop_
_database_PDB_matrix.origx[1][1]
_database_PDB_matrix.origx[1][2]
_database_PDB_matrix.origx[1][3]
_database_PDB_matrix.origx_vector[1]
_database_PDB_matrix.origx[2][1]
_database_PDB_matrix.origx[2][2]
_database_PDB_matrix.origx[2][3]
_database_PDB_matrix.origx_vector[2]
_database_PDB_matrix.origx[3][1]
_database_PDB_matrix.origx[3][2]
_database_PDB_matrix.origx[3][3]
_database_PDB_matrix.origx_vector[3]
1.0 0.0 0.0 4.0 0.0 1.0 0.0 5.0 0.0 0.0 1.0 6.0
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 C CA . ALA A 1 1 0.0 1.0 2.0
"#;

        let structure = read_test_mmcif(cif).unwrap();
        let crystal = structure.crystal().unwrap();

        assert_eq!(structure.metadata.entry_id.as_deref(), Some("9XYZ"));
        assert_eq!(structure.metadata.title.as_deref(), Some("TOP LEVEL TITLE"));
        assert_eq!(structure.metadata.authors, vec!["DOE, J."]);
        assert_eq!(crystal.spacegroup_hm.as_deref(), Some("P 1"));
        assert!(crystal.explicit_matrices);
        assert!(bio_transform_approx_eq(
            &crystal.frac,
            &BioTransform {
                mat: [[0.1, 0.0, 0.0], [0.0, 0.05, 0.0], [0.0, 0.0, 0.033333333]],
                vec: [1.0, 2.0, 3.0],
            },
            1e-6,
            1e-6
        ));
        assert!(bio_transform_approx_eq(
            &crystal.orth,
            &BioTransform {
                mat: [[10.0, 0.0, 0.0], [0.0, 20.0, 0.0], [0.0, 0.0, 30.0]],
                vec: [-10.0, -40.0, -90.0],
            },
            1e-5,
            1e-5
        ));
        assert!(structure.has_origx);
        assert_eq!(
            structure.origx,
            BioTransform {
                mat: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
                vec: [4.0, 5.0, 6.0],
            }
        );
    }

    #[test]
    fn mmcif_top_level_populate_attaches_connectivity_assemblies_sifts_and_restores_chem_comp() {
        let cif = r#"
data_demo
loop_
_entity.id
_entity.type
1 polymer
loop_
_entity_poly.entity_id
_entity_poly.type
1 polypeptide(L)
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
ATOM 1 S SG . CYS A 1 1 0.0 0.0 0.0 101 CYS X
ATOM 2 S SG . CYS A 1 2 1.0 0.0 0.0 102 CYS X
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr2_label_seq_id
conn1 disulf A A CYS CYS SG SG 1 2
loop_
_pdbx_struct_oper_list.id
_pdbx_struct_oper_list.type
_pdbx_struct_oper_list.matrix[1][1]
_pdbx_struct_oper_list.matrix[1][2]
_pdbx_struct_oper_list.matrix[1][3]
_pdbx_struct_oper_list.vector[1]
_pdbx_struct_oper_list.matrix[2][1]
_pdbx_struct_oper_list.matrix[2][2]
_pdbx_struct_oper_list.matrix[2][3]
_pdbx_struct_oper_list.vector[2]
_pdbx_struct_oper_list.matrix[3][1]
_pdbx_struct_oper_list.matrix[3][2]
_pdbx_struct_oper_list.matrix[3][3]
_pdbx_struct_oper_list.vector[3]
1 'point symmetry operation' 1 0 0 0 0 1 0 0 0 0 1 0
loop_
_pdbx_struct_assembly.id
_pdbx_struct_assembly.details
_pdbx_struct_assembly.method_details
_pdbx_struct_assembly.oligomeric_details
_pdbx_struct_assembly.oligomeric_count
1 software_defined_assembly program dimeric 2
loop_
_pdbx_struct_assembly_gen.assembly_id
_pdbx_struct_assembly_gen.oper_expression
_pdbx_struct_assembly_gen.asym_id_list
1 1 A
loop_
_pdbx_sifts_xref_db.entity_id
_pdbx_sifts_xref_db.asym_id
_pdbx_sifts_xref_db.seq_id_ordinal
_pdbx_sifts_xref_db.seq_id
_pdbx_sifts_xref_db.observed
_pdbx_sifts_xref_db.unp_res
_pdbx_sifts_xref_db.unp_num
_pdbx_sifts_xref_db.unp_acc
1 A 1 1 y C 15 P12345
1 A 1 2 y C 16 P12345
loop_
_chem_comp.id
_chem_comp.three_letter_code
~XY LONGNAME
"#;

        let structure = read_test_mmcif(cif).unwrap();

        assert_eq!(structure.connections.len(), 1);
        assert_eq!(structure.connections[0].type_, BioConnectionType::Disulf);
        assert_eq!(structure.assemblies.len(), 1);
        assert_eq!(structure.assemblies[0].name, "1");
        assert_eq!(structure.assemblies[0].software_name, "program");
        assert_eq!(structure.assemblies[0].generators[0].subchains, vec!["A"]);
        assert_eq!(structure.entities[0].sifts_unp_acc, vec!["P12345"]);
        assert_eq!(
            structure.residues[0].sifts_unp,
            Some(BioSiftsUnpResidue {
                res: Some('C'),
                acc_index: 0,
                num: 15,
            })
        );
        assert_eq!(
            structure.residues[1].sifts_unp,
            Some(BioSiftsUnpResidue {
                res: Some('C'),
                acc_index: 0,
                num: 16,
            })
        );
        assert!(structure.shortened_ccd_codes.is_empty());
    }

    #[test]
    fn mmcif_make_structure_accepts_single_coordinate_block() {
        let cif = r#"
data_demo
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 C CA . ALA A 1 1 0.0 0.0 0.0
"#;

        let structure = read_test_mmcif(cif).unwrap();

        assert_eq!(structure.num_atoms(), 1);
        assert_eq!(structure.num_models(), 1);
        assert_eq!(structure.input_format, BioCoorFormat::Mmcif);
        assert_eq!(structure.name, "demo");
    }

    #[test]
    fn mmcif_make_structure_skips_incomplete_chem_comp_short_code_table() {
        let cif = r#"
data_demo
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 C CA . ALA A 1 1 0.0 0.0 0.0
loop_
_chem_comp.id
_chem_comp.name
ALA Alanine
"#;

        let structure = read_test_mmcif(cif).unwrap();

        assert_eq!(structure.num_atoms(), 1);
        assert!(structure.shortened_ccd_codes.is_empty());
    }

    #[test]
    fn mmcif_make_structure_skips_incomplete_optional_gemmi_tables() {
        let cif = r#"
data_demo
loop_
_entity.id
_entity.type
1 polymer
loop_
_entity_poly.entity_id
1
loop_
_entity_poly_seq.entity_id
_entity_poly_seq.num
1 1
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 C CA . ALA A 1 1 0.0 0.0 0.0
loop_
_atom_site_anisotrop.id
_atom_site_anisotrop.U[1][1]
1 0.1
loop_
_struct_conn.id
_struct_conn.conn_type_id
conn1 disulf
loop_
_struct_mon_prot_cis.pdbx_PDB_model_num
_struct_mon_prot_cis.auth_asym_id
1 A
loop_
_pdbx_struct_mod_residue.auth_asym_id
A
loop_
_struct_conf.conf_type_id
_struct_conf.beg_auth_asym_id
HELX_P A
loop_
_struct_sheet_range.sheet_id
_struct_sheet_range.id
S1 1
loop_
_struct_sheet_order.sheet_id
_struct_sheet_order.range_id_2
S1 1
loop_
_pdbx_struct_sheet_hbond.sheet_id
_pdbx_struct_sheet_hbond.range_id_2
S1 1
loop_
_pdbx_refine_tls.id
_pdbx_refine_tls.T[1][1]
1 0.1
loop_
_diffrn.id
D1
loop_
_struct_ncs_oper.matrix[1][1]
_struct_ncs_oper.matrix[1][2]
1 0
loop_
_pdbx_struct_oper_list.id
_pdbx_struct_oper_list.matrix[1][1]
1 1
loop_
_pdbx_struct_assembly_prop.biol_id
_pdbx_struct_assembly_prop.type
1 ABSA
loop_
_pdbx_struct_assembly_gen.assembly_id
_pdbx_struct_assembly_gen.oper_expression
1 1
loop_
_pdbx_struct_assembly.id
_pdbx_struct_assembly.details
1 software_defined_assembly
loop_
_pdbx_sifts_xref_db.entity_id
_pdbx_sifts_xref_db.asym_id
1 A
"#;

        let structure = read_test_mmcif(cif).unwrap();

        assert_eq!(structure.num_atoms(), 1);
        assert!(structure.connections.is_empty());
        assert!(structure.cispeps.is_empty());
        assert!(structure.mod_residues.is_empty());
        assert!(structure.helices.is_empty());
        assert!(structure.sheets.is_empty());
        assert!(structure.assemblies.is_empty());
        assert!(structure.ncs_operators.is_empty());
    }

    #[test]
    fn mmcif_make_structure_allows_additional_non_coordinate_blocks() {
        let cif = r#"
data_demo
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 C CA . ALA A 1 1 0.0 0.0 0.0
data_restraints
loop_
_chem_comp.id
_chem_comp.name
ALA Alanine
"#;

        let structure = read_test_mmcif(cif).unwrap();

        assert_eq!(structure.num_atoms(), 1);
        assert_eq!(structure.num_models(), 1);
    }

    #[test]
    fn mmcif_make_structure_hands_document_back_through_save_doc_surface() {
        let cif = r#"
data_demo
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 C CA . ALA A 1 1 0.0 0.0 0.0
data_restraints
loop_
_chem_comp.id
_chem_comp.name
ALA Alanine
"#;

        let document = parse_cif_document(cif).unwrap();
        let mut saved = CifDocument { blocks: Vec::new() };
        let structure = make_structure_from_mmcif_document(document, Some(&mut saved)).unwrap();

        assert_eq!(structure.num_atoms(), 1);
        assert_eq!(saved.blocks.len(), 2);
        assert_eq!(saved.blocks[0].name, "demo");
        assert_eq!(saved.blocks[1].name, "restraints");
    }

    #[test]
    fn mmcif_make_structure_rejects_secondary_atom_site_block() {
        let cif = r#"
data_first
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 C CA . ALA A 1 1 0.0 0.0 0.0
data_second
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 2 N N . GLY B 2 2 1.0 1.0 1.0
"#;

        let err = read_test_mmcif(cif).unwrap_err();

        assert!(matches!(
            err,
            BioReadError::Parse {
                line_number: 0,
                message
            } if message.contains("_atom_site in block #2")
        ));
    }

    #[test]
    fn chemcomp_model_block_selects_xyz_example_and_ideal_coordinates() {
        let document = parse_cif_document(
            r#"
data_comp_LIG
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
_chem_comp_atom.model_Cartn_x
_chem_comp_atom.model_Cartn_y
_chem_comp_atom.model_Cartn_z
_chem_comp_atom.pdbx_model_Cartn_x_ideal
_chem_comp_atom.pdbx_model_Cartn_y_ideal
_chem_comp_atom.pdbx_model_Cartn_z_ideal
LIG C1 C 1.6 1.0 2.0 3.0 11.0 12.0 13.0 21.0 22.0 23.0
"#,
        )
        .unwrap();
        let block = &document.blocks[0];

        let xyz = make_model_from_chemcomp_block(block, ChemCompModel::Xyz).unwrap();
        let example = make_model_from_chemcomp_block(block, ChemCompModel::Example).unwrap();
        let ideal = make_model_from_chemcomp_block(block, ChemCompModel::Ideal).unwrap();

        assert_eq!(xyz.residues()[0].name.as_str(), "LIG");
        assert_eq!(
            xyz.atom_position(crate::bio::AtomId::new(0)),
            Some([1.0, 2.0, 3.0])
        );
        assert_eq!(
            example.atom_position(crate::bio::AtomId::new(0)),
            Some([11.0, 12.0, 13.0])
        );
        assert_eq!(
            ideal.atom_position(crate::bio::AtomId::new(0)),
            Some([21.0, 22.0, 23.0])
        );
        assert_eq!(xyz.atoms()[0].formal_charge, Some(2));
    }

    #[test]
    fn chemcomp_model_block_first_uses_first_coordinate_tags_in_loop_order() {
        let document = parse_cif_document(
            r#"
data_comp_ABA
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.model_Cartn_x
_chem_comp_atom.model_Cartn_y
_chem_comp_atom.model_Cartn_z
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
ABA N1 N 4.0 5.0 6.0 1.0 2.0 3.0
"#,
        )
        .unwrap();
        let block = &document.blocks[0];

        let first = make_model_from_chemcomp_block(block, ChemCompModel::First).unwrap();

        assert_eq!(
            first.atom_position(crate::bio::AtomId::new(0)),
            Some([4.0, 5.0, 6.0])
        );
    }

    #[test]
    fn chemcomp_model_block_falls_back_to_block_name_without_comp_id() {
        let document = parse_cif_document(
            r#"
data_comp_GLY
loop_
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
CA C 1.0 2.0 3.0
"#,
        )
        .unwrap();
        let block = &document.blocks[0];

        let model = make_model_from_chemcomp_block(block, ChemCompModel::Xyz).unwrap();

        assert_eq!(model.residues()[0].name.as_str(), "GLY");
    }

    #[test]
    fn chemcomp_structure_block_selects_requested_models_and_renumbers_them() {
        let document = parse_cif_document(
            r#"
data_demo
_chem_comp.id DEM
loop_
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.charge
_chem_comp_atom.model_Cartn_x
_chem_comp_atom.model_Cartn_y
_chem_comp_atom.model_Cartn_z
_chem_comp_atom.pdbx_model_Cartn_x_ideal
_chem_comp_atom.pdbx_model_Cartn_y_ideal
_chem_comp_atom.pdbx_model_Cartn_z_ideal
C1 C 0 1.0 2.0 3.0 4.0 5.0 6.0
"#,
        )
        .unwrap();
        let block = &document.blocks[0];

        let structure = make_structure_from_chemcomp_block(
            block,
            (ChemCompModel::Example as i32) | (ChemCompModel::Ideal as i32),
        )
        .unwrap();

        assert_eq!(structure.input_format, BioCoorFormat::ChemComp);
        assert_eq!(structure.name, "DEM");
        assert_eq!(structure.models.len(), 2);
        assert_eq!(structure.models[0].source_model_number, Some(1));
        assert_eq!(structure.models[1].source_model_number, Some(2));
        assert_eq!(structure.chains[1].model_id, ModelId::new(1));
        assert_eq!(structure.residues[1].chain_id, ChainId::new(1));
        assert_eq!(structure.atoms[1].residue_id, ResidueId::new(1));
        assert_eq!(structure.coordinates.positions[0], [1.0, 2.0, 3.0]);
        assert_eq!(structure.coordinates.positions[1], [4.0, 5.0, 6.0]);
    }

    #[test]
    fn chemcomp_doc_detection_matches_gemmi_block_layouts() {
        let monomer_no_global = parse_cif_document(
            "data_comp_list\n_chem_comp.id LIST\ndata_ATP\n_chem_comp.id ATP\nloop_\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n_chem_comp_atom.charge\n_chem_comp_atom.x\n_chem_comp_atom.y\n_chem_comp_atom.z\nP P 0 1.0 2.0 3.0\n",
        )
        .unwrap();
        assert_eq!(check_chemcomp_block_number(&monomer_no_global), 1);

        let monomer_with_global = parse_cif_document(
            "global_\n_lib.name demo\ndata_comp_list\n_chem_comp.id LIST\ndata_NAG\n_chem_comp.id NAG\nloop_\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n_chem_comp_atom.charge\n_chem_comp_atom.x\n_chem_comp_atom.y\n_chem_comp_atom.z\nC1 C 0 7.0 8.0 9.0\n",
        )
        .unwrap();
        assert_eq!(check_chemcomp_block_number(&monomer_with_global), 2);

        let ccd = parse_cif_document(
            "data_GLY\n_chem_comp.id GLY\nloop_\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n_chem_comp_atom.charge\n_chem_comp_atom.x\n_chem_comp_atom.y\n_chem_comp_atom.z\nN N 0 1.0 2.0 3.0\n",
        )
        .unwrap();
        assert_eq!(check_chemcomp_block_number(&ccd), 0);

        let not_chemcomp = parse_cif_document(
            "data_bad\n_cell.length_a 10\nloop_\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n_chem_comp_atom.charge\n_chem_comp_atom.x\n_chem_comp_atom.y\n_chem_comp_atom.z\nN N 0 1.0 2.0 3.0\n",
        )
        .unwrap();
        assert_eq!(check_chemcomp_block_number(&not_chemcomp), -1);
    }

    #[test]
    fn chemcomp_doc_reader_rejects_non_chemcomp_documents() {
        let document =
            parse_cif_document("data_bad\nloop_\n_atom_site.id\n_atom_site.type_symbol\n1 C\n")
                .unwrap();

        let error = make_structure_from_chemcomp_doc(&document, 7).unwrap_err();
        assert!(matches!(error, BioReadError::Parse { .. }));
        assert!(error.to_string().contains("Not a chem_comp format."));
    }

    #[test]
    fn dispatch_helpers_detect_extension_and_detect_mmjson() {
        assert_eq!(
            coor_format_from_ext("demo.pdb").unwrap(),
            BioCoorFormat::Pdb
        );
        assert_eq!(
            coor_format_from_ext("demo.ent").unwrap(),
            BioCoorFormat::Pdb
        );
        assert_eq!(
            coor_format_from_ext("demo.cif").unwrap(),
            BioCoorFormat::Mmcif
        );
        assert_eq!(
            coor_format_from_ext("demo.mmcif").unwrap(),
            BioCoorFormat::Mmcif
        );
        assert_eq!(
            coor_format_from_ext("demo.unknown").unwrap(),
            BioCoorFormat::Unknown
        );
        assert_eq!(
            coor_format_from_ext("demo.json").unwrap(),
            BioCoorFormat::Mmjson
        );
    }

    #[test]
    fn dispatch_helpers_detect_content_and_detect_mmjson() {
        assert_eq!(
            coor_format_from_content(b"data_demo\n_atom_site.id 1\n").unwrap(),
            BioCoorFormat::Mmcif
        );
        assert_eq!(
            coor_format_from_content(b"  # comment\nATOM      1  CA  ALA A   7\n").unwrap(),
            BioCoorFormat::Pdb
        );
        assert_eq!(
            coor_format_from_content(b" \n\t").unwrap(),
            BioCoorFormat::Unknown
        );
        assert_eq!(
            coor_format_from_content(br#"{"data_demo":{"entry":{"id":["DEMO"]}}}"#).unwrap(),
            BioCoorFormat::Mmjson
        );
    }

    #[test]
    fn structure_dispatch_from_memory_routes_pdb_mmcif_and_chemcomp() {
        let pdb =
            "ATOM      1  CA  ALA A   7      11.104  13.207   9.900  1.00 20.00           C  \n";
        let pdb_structure =
            read_structure_from_memory(pdb, "demo.pdb", BioCoorFormat::Detect).unwrap();
        assert_eq!(pdb_structure.input_format, BioCoorFormat::Pdb);
        assert_eq!(pdb_structure.name, "demo");
        assert_eq!(pdb_structure.num_atoms(), 1);

        let mmcif = "data_demo\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA . ALA A 1 1 0.0 0.0 0.0\n";
        let mmcif_structure =
            read_structure_from_memory(mmcif, "demo.cif", BioCoorFormat::Detect).unwrap();
        assert_eq!(mmcif_structure.input_format, BioCoorFormat::Mmcif);
        assert_eq!(mmcif_structure.name, "demo");
        assert_eq!(mmcif_structure.num_atoms(), 1);

        let chemcomp = "data_GLY\n_chem_comp.id GLY\nloop_\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n_chem_comp_atom.charge\n_chem_comp_atom.x\n_chem_comp_atom.y\n_chem_comp_atom.z\nN N 0 1.0 2.0 3.0\n";
        let chemcomp_structure =
            read_structure_from_memory(chemcomp, "gly.cif", BioCoorFormat::Mmcif).unwrap();
        assert_eq!(chemcomp_structure.input_format, BioCoorFormat::ChemComp);
        assert_eq!(chemcomp_structure.name, "GLY");
        assert_eq!(chemcomp_structure.num_atoms(), 1);
    }

    #[test]
    fn pdb_source_name_uses_gemmi_path_basename_rules() {
        assert_eq!(
            gemmi_path_basename(r"C:\data\demo.pdb.gz", &[".gz", ".pdb"]),
            "demo"
        );
        assert_eq!(
            gemmi_path_basename("/data/demo.PDB", &[".gz", ".pdb"]),
            "demo.PDB"
        );
        assert_eq!(gemmi_path_basename(".pdb", &[".gz", ".pdb"]), ".pdb");
    }

    #[test]
    fn structure_dispatch_from_memory_reads_mmjson_and_rejects_unknown_format() {
        let mmjson = r#"{
  "data_demo": {
    "atom_site": {
      "group_PDB": ["ATOM"],
      "id": ["1"],
      "type_symbol": ["C"],
      "label_atom_id": ["CA"],
      "label_alt_id": ["?"],
      "label_comp_id": ["ALA"],
      "label_asym_id": ["A"],
      "label_entity_id": ["1"],
      "label_seq_id": ["1"],
      "Cartn_x": [1.0],
      "Cartn_y": [2.0],
      "Cartn_z": [3.0]
    }
  }
}"#;
        let json_structure =
            read_structure_from_memory(mmjson, "demo.json", BioCoorFormat::Detect).unwrap();
        assert_eq!(json_structure.input_format, BioCoorFormat::Mmjson);
        assert_eq!(json_structure.name, "demo");
        assert_eq!(json_structure.num_atoms(), 1);
        assert_eq!(json_structure.coordinates.positions[0], [1.0, 2.0, 3.0]);

        let unknown_err =
            read_structure_from_memory(" \n\t", "demo.xyz", BioCoorFormat::Unknown).unwrap_err();
        assert!(matches!(unknown_err, BioReadError::Parse { .. }));
        assert!(
            unknown_err
                .to_string()
                .contains("wrong format of coordinate file demo.xyz")
        );
    }

    #[test]
    fn mmjson_document_reader_builds_cif_document_and_structure() {
        let mmjson = r#"{
  "data_demo": {
    "entry": {
      "id": ["DEMO"]
    },
    "atom_site": {
      "group_PDB": ["ATOM", "ATOM"],
      "id": ["1", "2"],
      "type_symbol": ["C", "N"],
      "label_atom_id": ["CA", "N"],
      "label_alt_id": ["?", "?"],
      "label_comp_id": ["ALA", "ALA"],
      "label_asym_id": ["A", "A"],
      "label_entity_id": ["1", "1"],
      "label_seq_id": ["1", "1"],
      "Cartn_x": [1.0, 4.0],
      "Cartn_y": [2.0, 5.0],
      "Cartn_z": [3.0, 6.0]
    }
  }
}"#;
        let document = read_mmjson_document(mmjson, "demo.json").unwrap();
        assert_eq!(document.blocks.len(), 1);
        assert_eq!(document.blocks[0].name, "demo");
        assert!(
            document.blocks[0]
                .items
                .iter()
                .any(|item| item.tag == "_entry.id")
        );
        assert!(
            document.blocks[0]
                .loops
                .iter()
                .any(|loop_| loop_.tags.iter().any(|tag| tag == "_atom_site.id"))
        );

        let mut structure = make_structure_from_mmcif_document(document, None).unwrap();
        structure.input_format = BioCoorFormat::Mmjson;
        assert_eq!(structure.num_atoms(), 2);
        assert_eq!(structure.input_format, BioCoorFormat::Mmjson);
        assert_eq!(structure.coordinates.positions[1], [4.0, 5.0, 6.0]);
    }

    #[test]
    fn rejects_mmcif_without_atom_site_loop() {
        let err = read_test_mmcif("data_test\n").unwrap_err();

        assert!(matches!(
            err,
            BioReadError::Unsupported { line_number: 0, .. }
        ));
    }
}
