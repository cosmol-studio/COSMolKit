//! PDB/mmCIF readers for `BioStructure`.
//!
//! Source-derived parser work in this module follows the source-reproduction
//! protocol documented in `io::sdf`: copied Gemmi source lines stay adjacent to
//! the Rust code that ports them, and every unsupported source branch remains
//! visible with a two-axis marker.
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
    AltLocLabel, AtomName, AtomRow, AtomSourceIds, BioStructure, ChainId, ChainKind, ChainRow,
    ChainSourceIds, ModelId, ModelRow, PdbAtomSerial, PdbChainId, PdbSeqId, ResidueId, ResidueKind,
    ResidueName, ResidueRow, ResidueSourceIds, RowSpan, classify_residue_name,
};
use crate::bio_invariants::enforce_bio_structure_invariants;
use crate::support::{
    BIO_MMCIF_ATOM_SITE_SUBSET_READ_FEATURE, BIO_PDB_COORDINATE_SUBSET_READ_FEATURE, FeatureSpec,
};

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
    /// Reject records whose Gemmi behavior has not been ported into this
    /// coordinate-subset reader.
    pub reject_unported_records: bool,
}

impl Default for BioPdbReadParams {
    fn default() -> Self {
        Self {
            reject_unported_records: false,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct ResidueKey {
    chain_key: PdbChainId,
    seq_id: PdbSeqId,
    name: ResidueName,
}

#[derive(Debug, Default)]
struct PdbBioBuilder {
    structure: BioStructure,
    current_model: Option<ModelId>,
    current_chain: Option<(ChainId, PdbChainId)>,
    current_residue: Option<(ResidueId, ResidueKey)>,
}

impl PdbBioBuilder {
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
        self.ensure_model(source_model_number);
    }

    fn end_model(&mut self) {
        self.current_model = None;
        self.current_chain = None;
        self.current_residue = None;
    }

    fn ensure_chain(
        &mut self,
        model_id: ModelId,
        chain_key: PdbChainId,
        source: ChainSourceIds,
    ) -> ChainId {
        if let Some((chain_id, current_source_id)) = self.current_chain
            && current_source_id == chain_key
        {
            return chain_id;
        }

        let chain_id = ChainId::new(self.structure.chains.len() as u32);
        let residue_start = self.structure.residues.len() as u32;
        self.structure.chains.push(ChainRow {
            model_id,
            entity_id: None,
            residue_span: RowSpan::new(residue_start, 0),
            kind: ChainKind::Unknown,
            source,
        });
        self.structure.models[model_id.index() as usize]
            .chain_span
            .len += 1;
        self.current_chain = Some((chain_id, chain_key));
        self.current_residue = None;
        chain_id
    }

    fn ensure_residue(
        &mut self,
        chain_id: ChainId,
        key: ResidueKey,
        kind: ResidueKind,
    ) -> ResidueId {
        if let Some((residue_id, current_key)) = self.current_residue
            && current_key == key
        {
            return residue_id;
        }

        let residue_id = ResidueId::new(self.structure.residues.len() as u32);
        let atom_start = self.structure.atoms.len() as u32;
        self.structure.residues.push(ResidueRow {
            chain_id,
            atom_span: RowSpan::new(atom_start, 0),
            name: key.name,
            kind,
            source: ResidueSourceIds {
                seq_id: Some(key.seq_id),
            },
        });
        self.structure.chains[chain_id.index() as usize]
            .residue_span
            .len += 1;
        self.current_residue = Some((residue_id, key));
        residue_id
    }

    fn push_atom(&mut self, record: PdbAtomRecord) {
        let model_id = self.ensure_model(None);
        let chain_id = self.ensure_chain(model_id, record.chain_key, record.chain_source);
        let residue_id = self.ensure_residue(chain_id, record.residue_key(), record.residue_kind);

        self.structure.atoms.push(AtomRow {
            residue_id,
            name: record.atom_name,
            element: record.element,
            altloc: record.altloc,
            occupancy: record.occupancy,
            b_iso: record.b_iso,
            formal_charge: record.formal_charge,
            anisou: None,
            source: AtomSourceIds {
                serial: record.serial.map(PdbAtomSerial),
            },
        });
        self.structure.coordinates.positions.push(record.position);
        self.structure.residues[residue_id.index() as usize]
            .atom_span
            .len += 1;
    }

    fn set_last_atom_anisou(
        &mut self,
        line_number: usize,
        anisou: [f32; 6],
    ) -> Result<(), BioReadError> {
        let atom = self
            .structure
            .atoms
            .last_mut()
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

    fn finish(self) -> Result<BioStructure, BioReadError> {
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
    position: [f32; 3],
    occupancy: Option<f32>,
    b_iso: Option<f32>,
    formal_charge: Option<i8>,
    element: Element,
}

impl PdbAtomRecord {
    fn residue_key(self) -> ResidueKey {
        ResidueKey {
            chain_key: self.chain_key,
            seq_id: self.seq_id,
            name: self.residue_name,
        }
    }
}

/// Reads the currently ported PDB coordinate subset into a `BioStructure`.
///
/// This is intentionally not named as a complete PDB reader. It models the
/// source-aligned coordinate-bearing subset documented by
/// `BIO_PDB_COORDINATE_SUBSET_READ_FEATURE`; unported PDB records remain
/// marked in this module's Gemmi source frame.
pub fn read_pdb_coordinate_subset_from_str(text: &str) -> Result<BioStructure, BioReadError> {
    read_pdb_coordinate_subset_from_str_with_params(text, BioPdbReadParams::default())
}

/// Reads the currently ported PDB coordinate subset with explicit subset policy.
pub fn read_pdb_coordinate_subset_from_str_with_params(
    text: &str,
    params: BioPdbReadParams,
) -> Result<BioStructure, BioReadError> {
    // BEGIN GEMMI CPP FUNCTION gemmi::populate_structure_from_pdb_stream
    // Gemmi✔️✔️: if (is_record_type4(line, "ATOM") || is_record_type4(line, "HETATM")) {
    // Gemmi✔️✔️:   if (len < 55)
    // Gemmi✔️✔️:     wrong("The line is too short to be correct:\n" + std::string(line));
    // Gemmi✔️✔️:   std::string chain_name = read_string(line+20, 2);
    // Gemmi✔️✔️:   ResidueId rid = read_res_id(line+22, line+17);
    // Gemmi❗✔️:   if (!chain || chain_name != chain->name) {
    // Gemmi❗✔️:     if (!model) {
    // Gemmi❗✔️:       int num = (int) st.models.size() + 1;
    // Gemmi❌❌:       if (st.find_model(num))
    // Gemmi❌❌:         wrong("ATOM/HETATM between models");
    // Gemmi❗✔️:       st.models.emplace_back(num);
    // Gemmi❗✔️:       model = &st.models.back();
    // Gemmi❌❌:     const Chain* prev_part = model->find_chain(chain_name);
    // Gemmi❌❌:     after_ter = prev_part &&
    // Gemmi❌❌:                   prev_part->residues[0].entity_type == EntityType::Polymer;
    // Gemmi❗✔️:     model->chains.emplace_back(chain_name);
    // Gemmi❗✔️:     chain = &model->chains.back();
    // Gemmi❗✔️:     resmap.clear();
    // Gemmi❗✔️:     resi = nullptr;
    // Gemmi❌❌:   if (len > 72)
    // Gemmi❌❌:     rid.segment = read_string(line+72, 4);
    // Gemmi❗✔️:   if (!resi || !resi->matches(rid)) {
    // Gemmi❗✔️:     auto it = resmap.find(rid);
    // Gemmi❗✔️:     if (it == resmap.end()) {
    // Gemmi❗✔️:       resmap.emplace(rid, (int) chain->residues.size());
    // Gemmi❗✔️:       chain->residues.emplace_back(rid);
    // Gemmi❗✔️:       resi = &chain->residues.back();
    // Gemmi❗✔️:       resi->het_flag = line[0] & ~0x20;
    // Gemmi❌❌:       if (after_ter)
    // Gemmi❌❌:         resi->entity_type = resi->is_water() ? EntityType::Water
    // Gemmi❌❌:                                                  : EntityType::NonPolymer;
    // Gemmi❗✔️:     } else {
    // Gemmi❗✔️:       resi = &chain->residues[it->second];
    // Gemmi❗✔️:   Atom atom;
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
    // Gemmi❗✔️:   else
    // Gemmi❗✔️:     atom.element = infer_element_from_padded_name(line+12);
    // Gemmi✔️✔️:   atom.charge = (len > 78 ? read_charge(line[78], line[79]) : 0);
    // Gemmi❗✔️:   resi->atoms.emplace_back(atom);
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
    // Gemmi❌❌: } else if (is_record_type4(line, "REMARK")) {
    // Gemmi❌❌: } else if (is_record_type4(line, "CONECT")) {
    // Gemmi❌❌: } else if (is_record_type4(line, "SEQRES")) {
    // Gemmi❌❌: } else if (is_record_type4(line, "HELIX")) {
    // Gemmi❌❌: } else if (is_record_type4(line, "SHEET")) {
    // Gemmi❌❌: } else if (is_record_type3(line, "TER") && !options.ignore_ter) {
    // END GEMMI CPP FUNCTION

    let mut builder = PdbBioBuilder::default();

    for (line_index, line) in text.lines().enumerate() {
        let line_number = line_index + 1;
        let record = record_type(line);
        match record {
            "ATOM" | "HETA" => {
                let atom = parse_pdb_atom_record(line, line_number)?;
                builder.push_atom(atom);
            }
            "MODE" if starts_record(line, "MODEL") => {
                let source_model_number = parse_optional_i32(field(line, 10, 14), line_number)?;
                builder.begin_model(source_model_number);
            }
            "ENDM" if starts_record(line, "ENDMDL") => builder.end_model(),
            "ANIS" if starts_record(line, "ANISOU") => {
                let anisou = parse_pdb_anisou_record(line, line_number)?;
                builder.set_last_atom_anisou(line_number, anisou)?;
            }
            "CONE" if starts_record(line, "CONECT") && params.reject_unported_records => {
                return Err(unsupported(
                    &BIO_PDB_COORDINATE_SUBSET_READ_FEATURE,
                    line_number,
                    "PDB CONECT bond semantics are not ported",
                ));
            }
            "SEQR" | "HELI" | "SHEE" | "SSBO" | "LINK" | "CISP" | "MODR" | "HETN" | "DBRE"
                if params.reject_unported_records =>
            {
                return Err(unsupported(
                    &BIO_PDB_COORDINATE_SUBSET_READ_FEATURE,
                    line_number,
                    "PDB metadata, sequence, secondary-structure, and connection records are not ported",
                ));
            }
            "TER " | "TER" if starts_record3(line, "TER") && params.reject_unported_records => {
                return Err(unsupported(
                    &BIO_PDB_COORDINATE_SUBSET_READ_FEATURE,
                    line_number,
                    "PDB TER chain/entity semantics are not ported",
                ));
            }
            _ => {}
        }
    }

    builder.finish()
}

/// Reads the currently ported mmCIF `_atom_site` subset into a `BioStructure`.
///
/// This is intentionally not named as a complete mmCIF reader. It only models
/// the atom-site coordinate subset documented by
/// `BIO_MMCIF_ATOM_SITE_SUBSET_READ_FEATURE`; broader mmCIF structure semantics
/// remain marked in this module's Gemmi source frame.
pub fn read_mmcif_atom_site_subset_from_str(text: &str) -> Result<BioStructure, BioReadError> {
    // BEGIN GEMMI CPP FUNCTION gemmi::make_structure
    // Gemmi❌❌: for (size_t i = 1; i < doc.blocks.size(); ++i)
    // Gemmi❌❌:   if (doc.blocks[i].has_tag("_atom_site.id"))
    // Gemmi❌❌:     fail("2+ blocks are ok if only the first one has coordinates;\n"
    // Gemmi❌❌:          "_atom_site in block #" + std::to_string(i+1) + ": " + doc.source);
    // Gemmi❌❌: Structure st = make_structure_from_block(doc.blocks.at(0));
    // Gemmi❌❌: if (save_doc)
    // Gemmi❌❌:   *save_doc = std::move(doc);
    // Gemmi❌❌: return st;
    // END GEMMI CPP FUNCTION

    // BEGIN GEMMI CPP FUNCTION gemmi::read_atom_sites
    // Gemmi❌❌: auto aniso_map = get_anisotropic_u(block);
    // Gemmi✔️✔️: cif::Table atom_table = block.find("_atom_site.",
    // Gemmi✔️✔️:                                {"id",
    // Gemmi✔️✔️:                                 "?group_PDB",
    // Gemmi✔️✔️:                                 "type_symbol",
    // Gemmi✔️✔️:                                 "?label_atom_id",
    // Gemmi✔️✔️:                                 "label_alt_id",
    // Gemmi✔️✔️:                                 "?label_comp_id",
    // Gemmi✔️✔️:                                 "label_asym_id",
    // Gemmi❌❌:                                 "?label_entity_id",
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
    // Gemmi❌❌:                                 "?calc_flag",
    // Gemmi❌❌:                                 "?pdbx_tls_group_id",
    // Gemmi❌❌:                                 "?ccp4_deuterium_fraction",
    // Gemmi✔️✔️: if (atom_table.length() != 0) {
    // Gemmi✔️✔️:   RowAccess asym_id(atom_table, kAuthAsymId, kLabelAsymId);
    // Gemmi✔️✔️:   RowAccess comp_id(atom_table, kAuthCompId, kLabelCompId);
    // Gemmi✔️✔️:   RowAccess atom_id(atom_table, kAuthAtomId, kLabelAtomId);
    // Gemmi✔️✔️:   RowAccess seq_id(atom_table, kAuthSeqId, kLabelSeqId);
    // Gemmi✔️✔️:   if (!asym_id.ok()) fail("Neither _atom_site.label_asym_id nor auth_asym_id found");
    // Gemmi✔️✔️:   if (!comp_id.ok()) fail("Neither _atom_site.label_comp_id nor auth_comp_id found");
    // Gemmi✔️✔️:   if (!atom_id.ok()) fail("Neither _atom_site.label_atom_id nor auth_atom_id found");
    // Gemmi✔️✔️:   if (!seq_id.ok()) fail("Neither _atom_site.label_seq_id nor auth_seq_id found");
    // Gemmi❗✔️:   for (auto row : atom_table) {
    // Gemmi✔️✔️:     atom.name = cif::as_string(atom_id.get(gap));
    // Gemmi✔️✔️:     atom.altloc = cif::as_char(row[kAltId], '\0');
    // Gemmi✔️✔️:     atom.charge = row.has2(kCharge) ? cif::as_int(row[kCharge]) : 0;
    // Gemmi✔️✔️:     atom.element = gemmi::Element(cif::as_string(row[kSymbol]));
    // Gemmi✔️✔️:     atom.serial = string_to_int(row[kId], false);
    // Gemmi❌❌:     if (st.has_d_fraction)
    // Gemmi❌❌:       atom.fraction = (float) cif::as_number(row[kDeuterium], 0.);
    // Gemmi❌❌:     if (row.has2(kCalcFlag)) { ... }
    // Gemmi❌❌:     if (row.has2(kTlsGroupId)) { ... }
    // Gemmi✔️✔️:     atom.pos.x = cif::as_number(row[kX]);
    // Gemmi✔️✔️:     atom.pos.y = cif::as_number(row[kY]);
    // Gemmi✔️✔️:     atom.pos.z = cif::as_number(row[kZ]);
    // Gemmi✔️✔️:     if (row.has2(kOcc))
    // Gemmi✔️✔️:       atom.occ = (float) cif::as_number(row[kOcc]);
    // Gemmi✔️✔️:     if (row.has2(kBiso))
    // Gemmi✔️✔️:       atom.b_iso = (float) cif::as_number(row[kBiso]);
    // Gemmi❌❌:     if (!aniso_map.empty()) { ... }
    // Gemmi❗✔️:     resi->atoms.emplace_back(atom);
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️: }
    // END GEMMI CPP FUNCTION

    let loops = parse_cif_loops(text)?;
    let atom_site = loops
        .iter()
        .find(|loop_| loop_.tags.iter().any(|tag| tag == "_atom_site.id"))
        .ok_or_else(|| {
            unsupported(
                &BIO_MMCIF_ATOM_SITE_SUBSET_READ_FEATURE,
                0,
                "mmCIF _atom_site loop is missing",
            )
        })?;
    read_mmcif_atom_site(atom_site)
}

fn parse_pdb_atom_record(line: &str, line_number: usize) -> Result<PdbAtomRecord, BioReadError> {
    if line.len() < 55 {
        return Err(BioReadError::Parse {
            line_number,
            message: "ATOM/HETATM line is too short to be correct".to_string(),
        });
    }

    let serial = parse_decimal_i32(
        &BIO_PDB_COORDINATE_SUBSET_READ_FEATURE,
        field(line, 6, 11),
        line_number,
        "atom serial",
    )?;
    let atom_name = atom_name_from_field(field_raw(line, 12, 16));
    let altloc = parse_altloc(byte_at(line, 16));
    let residue_name = residue_name_from_field(field(line, 17, 20));
    let residue_kind = classify_residue_name(residue_name);
    let chain_id = pdb_chain_id_from_field(field(line, 20, 22));
    let seq_num = parse_decimal_i32(
        &BIO_PDB_COORDINATE_SUBSET_READ_FEATURE,
        field(line, 22, 26),
        line_number,
        "residue sequence number",
    )?;
    let ins_code = match byte_at(line, 26) {
        b' ' | 0 => None,
        value => Some(value),
    };
    let x = parse_f32(field(line, 30, 38), line_number, "x coordinate")?;
    let y = parse_f32(field(line, 38, 46), line_number, "y coordinate")?;
    let z = parse_f32(field(line, 46, 54), line_number, "z coordinate")?;
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
        seq_id: PdbSeqId { seq_num, ins_code },
        position: [x, y, z],
        occupancy,
        b_iso,
        formal_charge,
        element,
    })
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
            &BIO_PDB_COORDINATE_SUBSET_READ_FEATURE,
            field(line, 28, 35),
            line_number,
            "ANISOU u11",
        )? as f32
            * 1e-4,
        parse_decimal_i32(
            &BIO_PDB_COORDINATE_SUBSET_READ_FEATURE,
            field(line, 35, 42),
            line_number,
            "ANISOU u22",
        )? as f32
            * 1e-4,
        parse_decimal_i32(
            &BIO_PDB_COORDINATE_SUBSET_READ_FEATURE,
            field(line, 42, 49),
            line_number,
            "ANISOU u33",
        )? as f32
            * 1e-4,
        parse_decimal_i32(
            &BIO_PDB_COORDINATE_SUBSET_READ_FEATURE,
            field(line, 49, 56),
            line_number,
            "ANISOU u12",
        )? as f32
            * 1e-4,
        parse_decimal_i32(
            &BIO_PDB_COORDINATE_SUBSET_READ_FEATURE,
            field(line, 56, 63),
            line_number,
            "ANISOU u13",
        )? as f32
            * 1e-4,
        parse_decimal_i32(
            &BIO_PDB_COORDINATE_SUBSET_READ_FEATURE,
            field(line, 63, 70),
            line_number,
            "ANISOU u23",
        )? as f32
            * 1e-4,
    ])
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct CifToken {
    value: String,
    line_number: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct CifLoop {
    tags: Vec<String>,
    values: Vec<CifToken>,
}

fn read_mmcif_atom_site(atom_site: &CifLoop) -> Result<BioStructure, BioReadError> {
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
    let mut builder = PdbBioBuilder::default();

    for row in atom_site.values.chunks(width) {
        let line_number = row.first().map_or(0, |token| token.line_number);
        let model_number = columns
            .model_num
            .and_then(|idx| cif_optional(row[idx].value.as_str()))
            .map(|value| {
                parse_decimal_i32(
                    &BIO_MMCIF_ATOM_SITE_SUBSET_READ_FEATURE,
                    value,
                    line_number,
                    "model number",
                )
            })
            .transpose()?;
        if builder.current_model.is_none()
            || model_number
                != builder.current_model.and_then(|id| {
                    builder.structure.models[id.index() as usize].source_model_number
                })
        {
            builder.begin_model(model_number.or(Some(1)));
        }

        let label_chain = cif_optional(row[columns.label_asym_id].value.as_str())
            .ok_or_else(|| missing_cif_value(line_number, "_atom_site.label_asym_id"))?;
        let auth_chain = columns
            .auth_asym_id
            .and_then(|idx| cif_optional(row[idx].value.as_str()));
        let chain_key = pdb_chain_id_from_cif(auth_chain.unwrap_or(label_chain), line_number)?;
        let label_chain_id = pdb_chain_id_from_cif(label_chain, line_number)?;
        let auth_chain_id = auth_chain
            .map(|value| pdb_chain_id_from_cif(value, line_number))
            .transpose()?;

        let atom_name = columns
            .auth_atom_id
            .and_then(|idx| cif_optional(row[idx].value.as_str()))
            .or_else(|| cif_optional(row[columns.label_atom_id].value.as_str()))
            .ok_or_else(|| {
                missing_cif_value(line_number, "_atom_site.label_atom_id/auth_atom_id")
            })?;
        let residue_name = columns
            .auth_comp_id
            .and_then(|idx| cif_optional(row[idx].value.as_str()))
            .or_else(|| cif_optional(row[columns.label_comp_id].value.as_str()))
            .ok_or_else(|| {
                missing_cif_value(line_number, "_atom_site.label_comp_id/auth_comp_id")
            })?;
        let seq_text = columns
            .auth_seq_id
            .and_then(|idx| cif_optional(row[idx].value.as_str()))
            .or_else(|| cif_optional(row[columns.label_seq_id].value.as_str()))
            .ok_or_else(|| missing_cif_value(line_number, "_atom_site.label_seq_id/auth_seq_id"))?;

        let serial = cif_optional(row[columns.id].value.as_str())
            .and_then(|value| value.parse::<i32>().ok());
        let seq_num = parse_decimal_i32(
            &BIO_MMCIF_ATOM_SITE_SUBSET_READ_FEATURE,
            seq_text,
            line_number,
            "mmCIF residue sequence number",
        )?;
        let ins_code = columns
            .ins_code
            .and_then(|idx| cif_optional(row[idx].value.as_str()))
            .and_then(|value| value.as_bytes().first().copied());
        let altloc = cif_optional(row[columns.alt_id].value.as_str())
            .and_then(|value| value.as_bytes().first().copied())
            .and_then(parse_altloc);
        let element = cif_optional(row[columns.type_symbol].value.as_str())
            .and_then(element_from_symbol)
            .ok_or_else(|| missing_cif_value(line_number, "_atom_site.type_symbol"))?;

        let x = parse_f32(row[columns.x].value.as_str(), line_number, "Cartn_x")?;
        let y = parse_f32(row[columns.y].value.as_str(), line_number, "Cartn_y")?;
        let z = parse_f32(row[columns.z].value.as_str(), line_number, "Cartn_z")?;
        let occupancy = columns
            .occupancy
            .and_then(|idx| cif_optional(row[idx].value.as_str()))
            .map(|value| parse_f32(value, line_number, "occupancy"))
            .transpose()?;
        let b_iso = columns
            .b_iso
            .and_then(|idx| cif_optional(row[idx].value.as_str()))
            .map(|value| parse_f32(value, line_number, "B_iso_or_equiv"))
            .transpose()?;
        let formal_charge = columns
            .formal_charge
            .and_then(|idx| cif_optional(row[idx].value.as_str()))
            .map(|value| parse_i8(value, line_number, "formal charge"))
            .transpose()?;
        let residue_name = residue_name_from_field(residue_name);

        builder.push_atom(PdbAtomRecord {
            serial,
            atom_name: atom_name_from_cif(atom_name),
            altloc,
            residue_name,
            residue_kind: classify_residue_name(residue_name),
            chain_key,
            chain_source: ChainSourceIds {
                auth_chain_id,
                label_asym_id: Some(label_chain_id),
            },
            seq_id: PdbSeqId { seq_num, ins_code },
            position: [x, y, z],
            occupancy,
            b_iso,
            formal_charge,
            element,
        });
    }

    builder.finish()
}

#[derive(Debug, Clone, Copy)]
struct AtomSiteColumns {
    id: usize,
    type_symbol: usize,
    label_atom_id: usize,
    alt_id: usize,
    label_comp_id: usize,
    label_asym_id: usize,
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
            type_symbol: required("_atom_site.type_symbol")?,
            label_atom_id: required("_atom_site.label_atom_id")?,
            alt_id: required("_atom_site.label_alt_id")?,
            label_comp_id: required("_atom_site.label_comp_id")?,
            label_asym_id: required("_atom_site.label_asym_id")?,
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
        })
    }
}

fn parse_cif_loops(text: &str) -> Result<Vec<CifLoop>, BioReadError> {
    let tokens = tokenize_cif(text)?;
    let mut loops = Vec::new();
    let mut idx = 0;
    while idx < tokens.len() {
        if !tokens[idx].value.eq_ignore_ascii_case("loop_") {
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
        loops.push(CifLoop { tags, values });
    }
    Ok(loops)
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
    field_raw(line, 0, 4)
}

fn starts_record(line: &str, record: &str) -> bool {
    line.get(0..record.len())
        .is_some_and(|prefix| prefix.eq_ignore_ascii_case(record))
}

fn starts_record3(line: &str, record: &str) -> bool {
    let prefix = field_raw(line, 0, 4);
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

fn byte_at(line: &str, index: usize) -> u8 {
    line.as_bytes().get(index).copied().unwrap_or(b' ')
}

fn parse_decimal_i32(
    feature: &'static FeatureSpec,
    value: &str,
    line_number: usize,
    field_name: &'static str,
) -> Result<i32, BioReadError> {
    let trimmed = value.trim();
    if trimmed.bytes().any(|byte| byte.is_ascii_alphabetic()) {
        return Err(unsupported(
            feature,
            line_number,
            "Gemmi hybrid-36 atom serial/residue sequence decoding is not ported",
        ));
    }
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
    parse_decimal_i32(
        &BIO_PDB_COORDINATE_SUBSET_READ_FEATURE,
        trimmed,
        line_number,
        "model number",
    )
    .map(Some)
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
            &BIO_MMCIF_ATOM_SITE_SUBSET_READ_FEATURE,
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
    let atomic_number = match symbol.trim().to_ascii_uppercase().as_str() {
        "H" => 1,
        "HE" => 2,
        "LI" => 3,
        "BE" => 4,
        "B" => 5,
        "C" => 6,
        "N" => 7,
        "O" => 8,
        "F" => 9,
        "NE" => 10,
        "NA" => 11,
        "MG" => 12,
        "AL" => 13,
        "SI" => 14,
        "P" => 15,
        "S" => 16,
        "CL" => 17,
        "AR" => 18,
        "K" => 19,
        "CA" => 20,
        "MN" => 25,
        "FE" => 26,
        "CO" => 27,
        "NI" => 28,
        "CU" => 29,
        "ZN" => 30,
        "SE" => 34,
        "BR" => 35,
        "I" => 53,
        _ => return None,
    };
    Element::from_atomic_number(atomic_number)
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
    fn reads_single_pdb_atom_record() {
        let pdb =
            "ATOM      1  CA  ALA A   7      11.104  13.207   9.900  1.00 20.00           C  \n";

        let structure = read_pdb_coordinate_subset_from_str(pdb).unwrap();

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
    fn preserves_altloc_and_water_kind() {
        let pdb =
            "HETATM   22  O  AHOH B  10       1.000   2.000   3.000  1.00 10.00           O  \n";

        let structure = read_pdb_coordinate_subset_from_str(pdb).unwrap();

        assert_eq!(structure.residues[0].kind, ResidueKind::Water);
        assert_eq!(structure.atoms[0].altloc, Some(AltLocLabel(b'A')));
    }

    #[test]
    fn reads_pdb_anisou_for_previous_atom() {
        let pdb = "\
ATOM      1  CA  ALA A   7      11.104  13.207   9.900  1.00 20.00           C  
ANISOU    1  CA  ALA A   7    1000   2000   3000    400    500    600       C  
";

        let structure = read_pdb_coordinate_subset_from_str(pdb).unwrap();

        assert_eq!(
            structure.atoms[0].anisou,
            Some([0.099999994, 0.19999999, 0.29999998, 0.04, 0.049999997, 0.06])
        );
    }

    #[test]
    fn infers_missing_pdb_element_from_padded_atom_name() {
        let pdb =
            "ATOM      1  CA  ALA A   7      11.104  13.207   9.900  1.00 20.00              \n";

        let structure = read_pdb_coordinate_subset_from_str(pdb).unwrap();

        assert_eq!(structure.atoms[0].element, Element::C);
    }

    #[test]
    fn strict_pdb_mode_rejects_unported_records() {
        let pdb = "\
ATOM      1  CA  ALA A   7      11.104  13.207   9.900  1.00 20.00           C  
CONECT    1    2
";

        let err = read_pdb_coordinate_subset_from_str_with_params(
            pdb,
            BioPdbReadParams {
                reject_unported_records: true,
            },
        )
        .unwrap_err();

        assert!(matches!(err, BioReadError::Unsupported { .. }));
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

        let structure = read_mmcif_atom_site_subset_from_str(cif).unwrap();

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
    fn rejects_mmcif_without_atom_site_loop() {
        let err = read_mmcif_atom_site_subset_from_str("data_test\n").unwrap_err();

        assert!(matches!(
            err,
            BioReadError::Unsupported { line_number: 0, .. }
        ));
    }
}
