//! `_atom_site` column discovery for the coordinate writer.

use cosmolkit_bio::{BioCalcFlag, BioStructureData};

use super::BioMmcifWriteParams;
use super::value::{
    atom_name_text, coordinate_text, float_field_text, pdbx_icode, qchain, residue_name_text,
    subchain_or_dot,
};
use crate::cif::{CifBlock, quote_cif_value};

/// The discovered `_atom_site` tag list and the scanned atom-row count.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct AtomSiteSchema {
    tags: Vec<String>,
    atom_site_count: usize,
}

impl AtomSiteSchema {
    /// Full `_atom_site.*` tags in source emission order.
    pub(crate) fn tags(&self) -> &[String] {
        &self.tags
    }

    /// Number of atom rows the single discovery scan counted.
    pub(crate) fn atom_site_count(&self) -> usize {
        self.atom_site_count
    }
}

/// Gemmi `add_cif_atoms` first phase: build the exact `_atom_site` tag list
/// and atom-row count with exactly one hierarchy scan.
pub(crate) fn discover_atom_site_schema(
    data: &BioStructureData,
    params: &BioMmcifWriteParams,
) -> AtomSiteSchema {
    // Gemmi✔️✔️: cif::Loop& atom_loop = block.init_mmcif_loop("_atom_site.", {
    // Gemmi✔️✔️:     "id",
    // Gemmi✔️✔️:     "type_symbol",
    // Gemmi✔️✔️:     "label_atom_id",
    // Gemmi✔️✔️:     "label_alt_id",
    // Gemmi✔️✔️:     "label_comp_id",
    // Gemmi✔️✔️:     "label_asym_id",
    // Gemmi✔️✔️:     "label_entity_id",
    // Gemmi✔️✔️:     "label_seq_id",
    // Gemmi✔️✔️:     "pdbx_PDB_ins_code",
    // Gemmi✔️✔️:     "Cartn_x",
    // Gemmi✔️✔️:     "Cartn_y",
    // Gemmi✔️✔️:     "Cartn_z",
    // Gemmi✔️✔️:     "occupancy",
    // Gemmi✔️✔️:     "B_iso_or_equiv",
    // Gemmi✔️✔️:     "pdbx_formal_charge",
    // Gemmi✔️✔️:     "auth_atom_id",  // optional (tags[15] is removed if !auth_all)
    // Gemmi✔️✔️:     "auth_comp_id",  // optional (tags[16] is removed if !auth_all)
    // Gemmi✔️✔️:     "auth_seq_id",
    // Gemmi✔️✔️:     "auth_asym_id",
    // Gemmi✔️✔️:     "pdbx_PDB_model_num"});
    // Gemmi✔️✔️: if (!auth_all)
    // Gemmi✔️✔️:   atom_loop.tags.erase(atom_loop.tags.begin() + 15, atom_loop.tags.begin() + 17);
    // Gemmi✔️✔️: if (use_group_pdb)
    // Gemmi✔️✔️:   atom_loop.tags.emplace(atom_loop.tags.begin(), "_atom_site.group_PDB");
    // Gemmi✔️✔️: bool has_calc_flag = false;
    // Gemmi✔️✔️: bool has_tls_group_id = false;
    // Gemmi✔️✔️: size_t atom_site_count = 0;
    // Gemmi✔️✔️: for (const Model& model : st.models)
    // Gemmi✔️✔️:   for (const Chain& chain : model.chains)
    // Gemmi✔️✔️:     for (const Residue& res : chain.residues)
    // Gemmi✔️✔️:       for (const Atom& atom : res.atoms) {
    // Gemmi✔️✔️:         ++atom_site_count;
    // Gemmi✔️✔️:         if (atom.calc_flag != CalcFlag::NotSet &&
    // Gemmi✔️✔️:             atom.calc_flag != CalcFlag::NoHydrogen)
    // Gemmi✔️✔️:           has_calc_flag = true;
    // Gemmi✔️✔️:         if (atom.tls_group_id >= 0)
    // Gemmi✔️✔️:           has_tls_group_id = true;
    // Gemmi✔️✔️:       }
    // Gemmi✔️✔️: if (has_calc_flag)
    // Gemmi✔️✔️:   atom_loop.tags.emplace_back("_atom_site.calc_flag");
    // Gemmi✔️✔️: if (has_tls_group_id)
    // Gemmi✔️✔️:   atom_loop.tags.emplace_back("_atom_site.pdbx_tls_group_id");
    // Gemmi✔️✔️: if (st.has_d_fraction)
    // Gemmi✔️✔️:   atom_loop.tags.emplace_back("_atom_site.ccp4_deuterium_fraction");
    // Behavior: the fixed base order is authoritative; the two auth columns
    // occupy positions 15/16 and vanish together when !auth_all; group_PDB
    // is prepended (position 0) when requested; calc_flag appears when any
    // atom carries a set flag other than NoHydrogen; the TLS column appears
    // when any atom has a non-negative (non--1) group id; the deuterium
    // fraction column follows the structure-level source_state flag exactly
    // as Gemmi's st.has_d_fraction. Appended columns keep this fixed tail
    // order. The hierarchy is scanned exactly once here; row emission and
    // per-atom schema questions never rescan.
    // Complexity: O(atoms) one pass plus O(tags) list edits; no per-atom
    // schema recomputation or repeated hierarchy walks.
    let mut tags: Vec<String> = [
        "id",
        "type_symbol",
        "label_atom_id",
        "label_alt_id",
        "label_comp_id",
        "label_asym_id",
        "label_entity_id",
        "label_seq_id",
        "pdbx_PDB_ins_code",
        "Cartn_x",
        "Cartn_y",
        "Cartn_z",
        "occupancy",
        "B_iso_or_equiv",
        "pdbx_formal_charge",
        "auth_atom_id",
        "auth_comp_id",
        "auth_seq_id",
        "auth_asym_id",
        "pdbx_PDB_model_num",
    ]
    .iter()
    .map(|suffix| format!("_atom_site.{suffix}"))
    .collect();
    if !params.auth_all {
        tags.drain(15..17);
    }
    if params.group_pdb {
        tags.insert(0, "_atom_site.group_PDB".to_string());
    }

    let mut has_calc_flag = false;
    let mut has_tls_group_id = false;
    let mut atom_site_count = 0usize;
    for atom in data.atoms() {
        atom_site_count += 1;
        let calc = atom.calc_flag();
        if !matches!(calc, BioCalcFlag::NotSet | BioCalcFlag::NoHydrogen) {
            has_calc_flag = true;
        }
        if atom.tls_group_id() >= 0 {
            has_tls_group_id = true;
        }
    }
    if has_calc_flag {
        tags.push("_atom_site.calc_flag".to_string());
    }
    if has_tls_group_id {
        tags.push("_atom_site.pdbx_tls_group_id".to_string());
    }
    if data.source_state().has_d_fraction {
        tags.push("_atom_site.ccp4_deuterium_fraction".to_string());
    }

    AtomSiteSchema {
        tags,
        atom_site_count,
    }
}

fn chain_text_string(chain: &cosmolkit_bio::BioChainRow) -> String {
    // Gemmi qchain(chain.name): the stored four-byte chain id spelled out.
    chain
        .source()
        .auth_chain_id()
        .map_or_else(String::new, |id| id.as_str().to_string())
}

/// Gemmi `pdbx_one_letter_code`-independent per-anisotropic-row record kept in
/// atom emission order; A03 consumes it for the `_atom_site_anisotrop` tail.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct AnisoRowRef {
    pub(crate) serial: i32,
    pub(crate) model_num: i32,
    pub(crate) atom_index: usize,
}

/// Gemmi `add_cif_atoms` row phase: emit `_atom_site` values aligned to the
/// A01 schema, plus the ordered anisotropic selection.
pub(crate) fn emit_atom_site_rows(
    data: &BioStructureData,
    params: &BioMmcifWriteParams,
    schema: &AtomSiteSchema,
) -> Result<(Vec<String>, Vec<AnisoRowRef>), super::BioMmcifWriteError> {
    // Gemmi✔️✔️:   for (const Model& model : st.models) {
    // Gemmi✔️✔️:     for (const Chain& chain : model.chains) {
    // Gemmi✔️✔️:       for (const Residue& res : chain.residues) {
    // Gemmi✔️✔️:         bool as_het = use_hetatm(res);
    // Gemmi✔️✔️:         std::string label_seq_id = res.label_seq.str('.');
    // Gemmi✔️✔️:         std::string auth_seq_id = res.seqid.num.str();
    // Gemmi✔️✔️:         std::string entity_id;
    // Gemmi✔️✔️:         if (const Entity* ent = gemmi::find_entity_of_subchain(res.subchain, st.entities))
    // Gemmi✔️✔️:           entity_id = cif::quote(ent->name);
    // Gemmi✔️✔️:         else
    // Gemmi✔️✔️:           entity_id = string_or_dot(res.entity_id);
    // Gemmi✔️✔️:         for (const Atom& atom : res.atoms) {
    // Gemmi✔️✔️:           if (use_group_pdb)
    // Gemmi✔️✔️:             vv.emplace_back(as_het ? "HETATM" : "ATOM");
    // Gemmi✔️✔️:           vv.emplace_back(std::to_string(++serial));
    // Gemmi✔️✔️:           vv.emplace_back(atom.element.uname());
    // Gemmi✔️✔️:           vv.emplace_back(cif::quote(atom.name));
    // Gemmi✔️✔️:           vv.emplace_back(1, atom.altloc_or('.'));
    // Gemmi✔️✔️:           vv.emplace_back(cif::quote(res.name));
    // Gemmi✔️✔️:           vv.emplace_back(subchain_or_dot(res));
    // Gemmi✔️✔️:           vv.emplace_back(entity_id);
    // Gemmi✔️✔️:           vv.emplace_back(label_seq_id);
    // Gemmi✔️✔️:           vv.emplace_back(pdbx_icode(res));
    // Gemmi✔️✔️:           vv.emplace_back(to_str(atom.pos.x));
    // Gemmi✔️✔️:           vv.emplace_back(to_str(atom.pos.y));
    // Gemmi✔️✔️:           vv.emplace_back(to_str(atom.pos.z));
    // Gemmi✔️✔️:           vv.emplace_back(to_str(atom.occ));
    // Gemmi✔️✔️:           vv.emplace_back(to_str(atom.b_iso));
    // Gemmi✔️✔️:           vv.emplace_back(atom.charge == 0 ? "?" : std::to_string(atom.charge));
    // Gemmi✔️✔️:           if (auth_all) {
    // Gemmi✔️✔️:             size_t atom_name_idx = vv.size() - 13;
    // Gemmi✔️✔️:             vv.emplace_back(vv[atom_name_idx]);  // auth_atom_id = label_atom_id
    // Gemmi✔️✔️:             vv.emplace_back(vv[atom_name_idx + 2]);  // auth_comp_id = label_comp_id
    // Gemmi✔️✔️:           }
    // Gemmi✔️✔️:           vv.emplace_back(auth_seq_id);
    // Gemmi✔️✔️:           vv.emplace_back(qchain(chain.name));
    // Gemmi✔️✔️:           vv.emplace_back(std::to_string(model.num));
    // Gemmi✔️✔️:           if (has_calc_flag)
    // Gemmi✔️✔️:             vv.emplace_back(&".\0.\0d\0c\0dum"[2 * (int) atom.calc_flag]);
    // Gemmi✔️✔️:           if (has_tls_group_id)
    // Gemmi✔️✔️:             vv.emplace_back(int_or_qmark(atom.tls_group_id));
    // Gemmi✔️✔️:           if (st.has_d_fraction)
    // Gemmi✔️✔️:             vv.emplace_back(to_str(atom.fraction));
    // Gemmi✔️✔️:           if (atom.aniso.nonzero())
    // Gemmi✔️✔️:             aniso.emplace_back(serial, model.num, &atom);
    // Gemmi✔️✔️:         }
    // Gemmi✔️✔️:       }
    // Gemmi✔️✔️:     }
    // Gemmi✔️✔️:   }
    // Behavior: values align one-to-one with schema.tags(); serial is a
    // single 1-based counter across the whole structure; label/auth sequence
    // ids follow Gemmi OptionalNum str rules (label '.', auth '?', CK PdbSeqId
    // stores the INT_MIN sentinel); altloc-or-'.' and charge-0-'?' are
    // verbatim; coordinates use the double profile and occ/B/fraction/aniso
    // the float profile; auth_atom_id/auth_comp_id duplicate the exact
    // label strings including quoting; the conditional tail mirrors the
    // schema flag order; anisotropic rows are collected in emission order
    // keyed by serial, model number and atom index. Gemmi Model.num defaults
    // to 0; CK readers always populate source_model_number, and a None only
    // arises for programmatically built parts, where unwrap_or_default keeps
    // the Gemmi default rather than inventing a number. Element uname() is
    // the uppercase canonical symbol, matching Element::symbol().
    // Complexity: one hierarchy walk, one output String per column per atom
    // and constant per-row work; no repeated schema or entity rescans (the
    // per-residue entity id is looked up once per residue, like the source).
    let mut values: Vec<String> =
        Vec::with_capacity(schema.tags().len() * schema.atom_site_count());
    let mut aniso: Vec<AnisoRowRef> = Vec::new();
    let has_calc_flag = schema
        .tags()
        .iter()
        .any(|tag| tag == "_atom_site.calc_flag");
    let has_tls_group_id = schema
        .tags()
        .iter()
        .any(|tag| tag == "_atom_site.pdbx_tls_group_id");
    let has_d_fraction = schema
        .tags()
        .iter()
        .any(|tag| tag == "_atom_site.ccp4_deuterium_fraction");
    let input_format = data.input_format();
    let positions = data.coordinates().positions();
    let mut serial: i32 = 0;
    for model in data.models() {
        let model_num = model.source_model_number().unwrap_or_default();
        let model_num_text = model_num.to_string();
        for chain_index in model.chain_span().start() as usize..model.chain_span().end() as usize {
            let chain = &data.chains()[chain_index];
            let chain_text = qchain(&chain_text_string(chain));
            for residue_index in
                chain.residue_span().start() as usize..chain.residue_span().end() as usize
            {
                let residue = &data.residues()[residue_index];
                let as_het = super::tags::use_hetatm(residue);
                let label_seq_id = residue
                    .source()
                    .label_seq_id()
                    .map_or_else(|| ".".to_string(), |value| value.to_string());
                let auth_seq_id = residue.source().seq_id().map_or_else(
                    || "?".to_string(),
                    // Gemmi✔️✔️: bool has_value() const { return value != None; }
                    // Gemmi✔️✔️: std::string str(char null='?') const {
                    // Gemmi✔️✔️:   return has_value() ? std::to_string(value) : std::string(1, null);
                    // Gemmi✔️✔️: }
                    |seq| {
                        if seq.seq_num() == i32::MIN {
                            "?".to_string()
                        } else {
                            seq.seq_num().to_string()
                        }
                    },
                );
                let entity_id = super::tags::entity_id_for_residue(residue, data.entities())?;
                let residue_name =
                    quote_cif_value(residue_name_text(&residue.name(), input_format).to_string());
                let subchain = subchain_or_dot(residue.source());
                let icode = pdbx_icode(residue.source().seq_id());
                for atom_index in
                    residue.atom_span().start() as usize..residue.atom_span().end() as usize
                {
                    let atom = &data.atoms()[atom_index];
                    serial = serial.checked_add(1).ok_or_else(|| {
                        super::BioMmcifWriteError::Structure(
                            cosmolkit_bio::BioStructureError::RowIndexTooLarge {
                                value: usize::MAX,
                            },
                        )
                    })?;
                    if params.group_pdb {
                        values.push(if as_het { "HETATM" } else { "ATOM" }.to_string());
                    }
                    values.push(serial.to_string());
                    values.push(atom.element().symbol().to_string());
                    let label_atom = quote_cif_value(atom_name_text(atom.name(), input_format));
                    values.push(label_atom.clone());
                    values.push(
                        (atom.altloc().map_or(b'.', |label| label.value()) as char).to_string(),
                    );
                    values.push(residue_name.clone());
                    values.push(subchain.clone());
                    values.push(entity_id.clone());
                    values.push(label_seq_id.clone());
                    values.push(icode.clone());
                    let position = positions
                        .get(atom_index)
                        .map_or([f64::NAN, f64::NAN, f64::NAN], |row| *row);
                    values.push(coordinate_text(position[0]));
                    values.push(coordinate_text(position[1]));
                    values.push(coordinate_text(position[2]));
                    values.push(float_field_text(atom.occupancy()));
                    values.push(float_field_text(atom.b_iso()));
                    values.push(match atom.formal_charge() {
                        0 => "?".to_string(),
                        charge => charge.to_string(),
                    });
                    if params.auth_all {
                        values.push(label_atom);
                        values.push(residue_name.clone());
                    }
                    values.push(auth_seq_id.clone());
                    values.push(chain_text.clone());
                    values.push(model_num_text.clone());
                    if has_calc_flag {
                        values.push(super::tags::calc_flag_text(atom.calc_flag()).to_string());
                    }
                    if has_tls_group_id {
                        values.push(super::value::int_or_qmark(match atom.tls_group_id() {
                            value if value >= 0 => Some(value as i32),
                            _ => None,
                        }));
                    }
                    if has_d_fraction {
                        values.push(float_field_text(atom.fraction()));
                    }
                    if atom.anisou().iter().any(|value| *value != 0.0) {
                        aniso.push(AnisoRowRef {
                            serial,
                            model_num,
                            atom_index,
                        });
                    }
                }
            }
        }
    }
    Ok((values, aniso))
}

/// Gemmi `add_cif_atoms` (complete): write the `_atom_site` loop from the
/// A01 schema and A02 row emission, then the conditional
/// `_atom_site_anisotrop` tail, into a mutable CIF block.
///
/// The two phases stay separate functions because the source splits tag
/// discovery/emission only conceptually; composition here follows
/// to_mmcif.cpp:164-187 verbatim for the tail.
pub(crate) fn add_cif_atoms(
    block: &mut CifBlock,
    data: &BioStructureData,
    params: &BioMmcifWriteParams,
) -> Result<(), super::BioMmcifWriteError> {
    // Gemmi✔️✔️:   if (aniso.empty()) {
    // Gemmi✔️✔️:     block.find_mmcif_category("_atom_site_anisotrop.").erase();
    // Gemmi✔️✔️:   } else {
    // Gemmi✔️✔️:     cif::Loop& aniso_loop = block.init_mmcif_loop("_atom_site_anisotrop.", {
    // Gemmi✔️✔️:                                   "id", "type_symbol", "U[1][1]", "U[2][2]",
    // Gemmi✔️✔️:                                   "U[3][3]", "U[1][2]", "U[1][3]", "U[2][3]"});
    // Gemmi✔️✔️:     if (st.models.size() > 1)
    // Gemmi✔️✔️:       aniso_loop.tags.push_back("_atom_site_anisotrop.pdbx_PDB_model_num");
    // Gemmi✔️✔️:     std::vector<std::string>& aniso_val = aniso_loop.values;
    // Gemmi✔️✔️:     aniso_val.reserve(aniso_loop.tags.size() * aniso.size());
    // Gemmi✔️✔️:     for (const auto& a : aniso) {
    // Gemmi✔️✔️:       aniso_val.emplace_back(std::to_string(std::get<0>(a)));
    // Gemmi✔️✔️:       const Atom* atom = std::get<2>(a);
    // Gemmi✔️✔️:       aniso_val.emplace_back(atom->element.uname());
    // Gemmi✔️✔️:       aniso_val.emplace_back(to_str(atom->aniso.u11));
    // Gemmi✔️✔️:       aniso_val.emplace_back(to_str(atom->aniso.u22));
    // Gemmi✔️✔️:       aniso_val.emplace_back(to_str(atom->aniso.u33));
    // Gemmi✔️✔️:       aniso_val.emplace_back(to_str(atom->aniso.u12));
    // Gemmi✔️✔️:       aniso_val.emplace_back(to_str(atom->aniso.u13));
    // Gemmi✔️✔️:       aniso_val.emplace_back(to_str(atom->aniso.u23));
    // Gemmi✔️✔️:       if (st.models.size() > 1)
    // Gemmi✔️✔️:         aniso_loop.values.push_back(std::to_string(std::get<1>(a)));
    // Gemmi✔️✔️:     }
    // Gemmi✔️✔️:   }
    // Behavior: the `_atom_site` loop is initialized with the A01 schema
    // suffixes and filled with the A02 flat value vector in one bulk move
    // (the source writes directly into `atom_loop.values`). The anisotropic
    // tail is erased when the A02 selection is empty; otherwise it carries
    // the fixed eight base columns in source order plus
    // `pdbx_PDB_model_num` only when the structure has more than one model.
    // Each anisotropic row links back to the emitting atom through the A02
    // serial (id) and, for multimodel structures, the Gemmi model number;
    // the type symbol and the six SMat33 components (u11, u22, u33, u12,
    // u13, u23 — the same order the mmCIF reader stores) use the float
    // profile because Gemmi's aniso is SMat33<float>.
    // Complexity: one schema discovery, one row emission and one linear tail
    // pass; the atom values move once with no per-row copies, and the empty
    // tail is a single category erase like the source.
    let schema = discover_atom_site_schema(data, params);
    let (values, aniso) = emit_atom_site_rows(data, params, &schema)?;
    let suffixes: Vec<&str> = schema
        .tags()
        .iter()
        .map(|tag| tag.strip_prefix("_atom_site.").unwrap_or(tag.as_str()))
        .collect();
    block
        .init_mmcif_loop("_atom_site.", &suffixes)?
        .set_string_values(values)?;

    if aniso.is_empty() {
        block.erase_mmcif_category("_atom_site_anisotrop.");
    } else {
        let multiple_models = data.models().len() > 1;
        let aniso_loop = block.init_mmcif_loop(
            "_atom_site_anisotrop.",
            &[
                "id",
                "type_symbol",
                "U[1][1]",
                "U[2][2]",
                "U[3][3]",
                "U[1][2]",
                "U[1][3]",
                "U[2][3]",
            ],
        )?;
        if multiple_models {
            aniso_loop.push_tag("_atom_site_anisotrop.pdbx_PDB_model_num".to_string());
        }
        let width = aniso_loop.tags().len();
        let mut aniso_values = Vec::with_capacity(width * aniso.len());
        for row in &aniso {
            let atom = &data.atoms()[row.atom_index];
            aniso_values.push(row.serial.to_string());
            aniso_values.push(atom.element().symbol().to_string());
            aniso_values.push(float_field_text(atom.anisou()[0]));
            aniso_values.push(float_field_text(atom.anisou()[1]));
            aniso_values.push(float_field_text(atom.anisou()[2]));
            aniso_values.push(float_field_text(atom.anisou()[3]));
            aniso_values.push(float_field_text(atom.anisou()[4]));
            aniso_values.push(float_field_text(atom.anisou()[5]));
            if multiple_models {
                aniso_values.push(row.model_num.to_string());
            }
        }
        aniso_loop.set_string_values(aniso_values)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{AnisoRowRef, AtomSiteSchema, discover_atom_site_schema, emit_atom_site_rows};
    use crate::bio_write::BioMmcifWriteParams;
    use cosmolkit_bio::{
        AtomName, AtomSourceIds, BioAtomRow, BioCalcFlag, BioChainId, BioChainRow,
        BioCoordinateBlock, BioCoordinateFormat, BioModelId, BioModelRow, BioResidueId,
        BioResidueRow, BioRowSpan, BioSiftsUnpResidue, BioStructureData, BioStructureParts,
        BioStructureSourceState, ChainKind, ChainSourceIds, EntityKind, ResidueInfoKind,
        ResidueName, ResidueSourceIds,
    };
    use cosmolkit_types::Element;

    fn atom(calc_flag: BioCalcFlag, tls_group_id: i16) -> BioAtomRow {
        atom_in_residue(0, calc_flag, tls_group_id)
    }

    fn atom_in_residue(residue: u32, calc_flag: BioCalcFlag, tls_group_id: i16) -> BioAtomRow {
        flex_atom(residue, calc_flag, tls_group_id, None, 0, false)
    }

    #[allow(clippy::too_many_arguments)]
    fn flex_atom(
        residue: u32,
        calc_flag: BioCalcFlag,
        tls_group_id: i16,
        altloc: Option<u8>,
        charge: i8,
        anisou: bool,
    ) -> BioAtomRow {
        BioAtomRow::new(
            BioResidueId::new(residue),
            AtomName::from_ascii(b"CA").unwrap(),
            Element::C,
            None,
            altloc.map(cosmolkit_bio::AltLocLabel::new),
            charge,
            calc_flag,
            1.0,
            20.0,
            if anisou {
                [0.05, 0.0, 0.0, 0.0, 0.0, 0.0]
            } else {
                [0.0; 6]
            },
            tls_group_id,
            0.0,
            AtomSourceIds::default(),
        )
    }

    fn entity(name: &str, kind: EntityKind, subchains: &[&str]) -> cosmolkit_bio::BioEntityRow {
        cosmolkit_bio::BioEntityRow::new(
            kind,
            cosmolkit_bio::PolymerKind::Unknown,
            false,
            Vec::new(),
            Vec::new(),
            Vec::new(),
            subchains.iter().map(|value| value.to_string()).collect(),
            cosmolkit_bio::EntitySourceIds::new(name.to_string()),
        )
    }

    fn empty_parts() -> BioStructureParts {
        BioStructureParts {
            input_format: BioCoordinateFormat::Unknown,
            models: vec![],
            chains: vec![],
            residues: vec![],
            atoms: vec![],
            entities: vec![],
            connections: vec![],
            cispeps: vec![],
            mod_residues: vec![],
            helices: vec![],
            sheets: vec![],
            metadata: Default::default(),
            source_state: Default::default(),
            coordinates: BioCoordinateBlock::default(),
            crystal: None,
            ncs_operators: vec![],
            assemblies: vec![],
        }
    }

    fn structure(atoms: Vec<BioAtomRow>, has_d_fraction: bool) -> BioStructureData {
        let atom_count = atoms.len() as u32;
        let mut parts = empty_parts();
        parts.models = vec![BioModelRow::new(BioRowSpan::new(0, 1).unwrap(), Some(1))];
        parts.chains = vec![BioChainRow::new(
            BioModelId::new(0),
            None,
            BioRowSpan::new(0, 1).unwrap(),
            ChainKind::Protein,
            ChainSourceIds::default(),
        )];
        parts.residues = vec![BioResidueRow::new(
            BioChainId::new(0),
            BioRowSpan::new(0, atom_count).unwrap(),
            ResidueName::from_ascii(b"ALA").unwrap(),
            ResidueInfoKind::Aa,
            EntityKind::Polymer,
            None,
            None,
            ResidueSourceIds::default(),
            BioSiftsUnpResidue::default(),
        )];
        parts.atoms = atoms;
        parts.coordinates = BioCoordinateBlock::new(vec![[1.0, 2.0, 3.0]; parts.atoms.len()]);
        parts.source_state = BioStructureSourceState {
            has_d_fraction,
            ..BioStructureSourceState::default()
        };
        BioStructureData::from_parts(parts).unwrap()
    }

    fn suffixes(schema: &AtomSiteSchema) -> Vec<&str> {
        schema
            .tags()
            .iter()
            .map(|tag| tag.strip_prefix("_atom_site.").unwrap())
            .collect()
    }

    #[test]
    fn bio_pdbscope_a02_rows_match_schema_order_serial_and_profiles() {
        // Two-chain/two-model fixture with author/label differences, altloc,
        // charge, anisou and a non-trivial entity: every column order and
        // value spelling below is derived from to_mmcif.cpp:115-160, not
        // from CK output.
        let mut parts = empty_parts();
        parts.input_format = BioCoordinateFormat::Mmcif;
        parts.models = vec![
            BioModelRow::new(BioRowSpan::new(0, 1).unwrap(), Some(1)),
            BioModelRow::new(BioRowSpan::new(1, 1).unwrap(), Some(3)),
        ];
        parts.chains = vec![
            BioChainRow::new(
                BioModelId::new(0),
                None,
                BioRowSpan::new(0, 1).unwrap(),
                ChainKind::Protein,
                ChainSourceIds::new(
                    Some(cosmolkit_bio::PdbChainId::from_ascii(b"A").unwrap()),
                    Some("subA".to_string()),
                ),
            ),
            BioChainRow::new(
                BioModelId::new(1),
                None,
                BioRowSpan::new(1, 1).unwrap(),
                ChainKind::Protein,
                ChainSourceIds::new(
                    Some(cosmolkit_bio::PdbChainId::from_ascii(b"B").unwrap()),
                    Some("subB".to_string()),
                ),
            ),
        ];
        parts.residues = vec![
            BioResidueRow::new(
                BioChainId::new(0),
                BioRowSpan::new(0, 2).unwrap(),
                ResidueName::from_ascii(b"ALA").unwrap(),
                ResidueInfoKind::Aa,
                EntityKind::Polymer,
                None,
                None,
                ResidueSourceIds::new(
                    Some(cosmolkit_bio::PdbSeqId::new(12, Some(b'A'))),
                    Some(7),
                    None,
                    Some("subA".to_string()),
                    None,
                )
                .unwrap(),
                BioSiftsUnpResidue::default(),
            ),
            BioResidueRow::new(
                BioChainId::new(1),
                BioRowSpan::new(2, 1).unwrap(),
                ResidueName::from_ascii(b"HOH").unwrap(),
                ResidueInfoKind::Paa,
                EntityKind::Water,
                None,
                Some(b'H'),
                ResidueSourceIds::new(
                    Some(cosmolkit_bio::PdbSeqId::new(i32::MIN, None)),
                    None,
                    None,
                    Some("subB".to_string()),
                    None,
                )
                .unwrap(),
                BioSiftsUnpResidue::default(),
            ),
        ];
        parts.entities = vec![entity("entB", EntityKind::Water, &["subB"])];
        parts.atoms = vec![
            flex_atom(0, BioCalcFlag::Determined, 2, Some(b'B'), -1, true),
            flex_atom(0, BioCalcFlag::NotSet, -1, None, 0, false),
            flex_atom(1, BioCalcFlag::NotSet, -1, None, 0, false),
        ];
        parts.coordinates = BioCoordinateBlock::new(vec![
            [11.111, 2.222, 3.333],
            [-0.0, 5.0, 6.0],
            [9.5, 8.5, 7.5],
        ]);
        parts.source_state.has_d_fraction = true;
        let data = BioStructureData::from_parts(parts).unwrap();

        for (group_pdb, auth_all) in [(false, false), (true, true)] {
            let params = BioMmcifWriteParams {
                group_pdb,
                auth_all,
                ..BioMmcifWriteParams::default()
            };
            let schema = discover_atom_site_schema(&data, &params);
            let (values, aniso) = emit_atom_site_rows(&data, &params, &schema).unwrap();
            assert_eq!(values.len(), schema.tags().len() * schema.atom_site_count());
            let width = schema.tags().len();
            let row = |index: usize| -> Vec<&str> {
                values[index * width..(index + 1) * width]
                    .iter()
                    .map(String::as_str)
                    .collect()
            };
            // Row 0: ALA chain A, altloc B, charge -1, calc d, tls 2, fraction 0
            let first = row(0);
            assert_eq!(first[0], if group_pdb { "ATOM" } else { "1" });
            let base = if group_pdb { 1 } else { 0 };
            assert_eq!(first[base], "1", "serial is a global 1-based counter");
            assert_eq!(first[base + 1], "C");
            assert_eq!(first[base + 2], "CA");
            assert_eq!(first[base + 3], "B");
            assert_eq!(first[base + 4], "ALA");
            assert_eq!(first[base + 5], "subA");
            assert_eq!(first[base + 6], ".", "no entity owns subA");
            assert_eq!(first[base + 7], "7", "label_seq_id prints the number");
            assert_eq!(first[base + 8], "A", "insertion code");
            assert_eq!(first[base + 9], "11.111");
            assert_eq!(first[base + 10], "2.222");
            assert_eq!(first[base + 11], "3.333");
            assert_eq!(first[base + 12], "1");
            assert_eq!(first[base + 13], "20");
            assert_eq!(first[base + 14], "-1");
            if auth_all {
                assert_eq!(first[base + 15], "CA", "auth_atom_id duplicates label");
                assert_eq!(first[base + 16], "ALA", "auth_comp_id duplicates label");
            }
            let auth_offset = base + if auth_all { 17 } else { 15 };
            assert_eq!(first[auth_offset], "12", "auth seq prints the number");
            assert_eq!(first[auth_offset + 1], "A");
            assert_eq!(first[auth_offset + 2], "1", "model 1");
            assert_eq!(first[auth_offset + 3], "d", "Determined -> d");
            assert_eq!(first[auth_offset + 4], "2", "tls_group_id 2");
            assert_eq!(first[auth_offset + 5], "0", "fraction float profile");
            // Row 1: same residue, defaults
            let second = row(1);
            assert_eq!(second[base], "2");
            assert_eq!(second[base + 3], ".", "no altloc prints '.'");
            assert_eq!(second[base + 14], "?", "charge 0 prints '?'");
            assert_eq!(second[auth_offset + 3], ".", "NotSet -> '.'");
            assert_eq!(second[auth_offset + 4], "?", "tls -1 -> '?'");
            // Row 2: model 3, HETATM via water entity kind + het_flag H,
            // INT_MIN auth seq and missing label seq
            let third = row(2);
            assert_eq!(third[base], "3");
            if group_pdb {
                assert_eq!(third[0], "HETATM");
            }
            assert_eq!(third[base + 4], "HOH");
            assert_eq!(third[base + 5], "subB");
            assert_eq!(third[base + 6], "entB", "entity name is quoted");
            assert_eq!(third[base + 7], ".", "absent label_seq_id prints '.'");
            assert_eq!(third[base + 8], "?", "no insertion code");
            assert_eq!(third[auth_offset], "?", "INT_MIN auth seq prints '?'");
            assert_eq!(third[auth_offset + 2], "3", "model 3");
            assert_eq!(
                aniso,
                vec![AnisoRowRef {
                    serial: 1,
                    model_num: 1,
                    atom_index: 0
                }]
            );
        }
    }

    #[test]
    fn bio_pdbscope_a01_four_option_combinations_and_base_order() {
        let plain_atoms = structure(vec![atom(BioCalcFlag::NotSet, -1)], false);
        for (group_pdb, auth_all) in [(false, false), (false, true), (true, false), (true, true)] {
            let params = BioMmcifWriteParams {
                group_pdb,
                auth_all,
                ..BioMmcifWriteParams::default()
            };
            let schema = discover_atom_site_schema(&plain_atoms, &params);
            let mut expected: Vec<&str> = [
                "id",
                "type_symbol",
                "label_atom_id",
                "label_alt_id",
                "label_comp_id",
                "label_asym_id",
                "label_entity_id",
                "label_seq_id",
                "pdbx_PDB_ins_code",
                "Cartn_x",
                "Cartn_y",
                "Cartn_z",
                "occupancy",
                "B_iso_or_equiv",
                "pdbx_formal_charge",
            ]
            .to_vec();
            if auth_all {
                expected.push("auth_atom_id");
                expected.push("auth_comp_id");
            }
            expected.extend(["auth_seq_id", "auth_asym_id", "pdbx_PDB_model_num"]);
            if group_pdb {
                expected.insert(0, "group_PDB");
            }
            assert_eq!(suffixes(&schema), expected);
            assert_eq!(schema.atom_site_count(), 1);
        }
    }

    #[test]
    fn bio_pdbscope_a01_optional_columns_empty_and_multiple_models() {
        let none = BioStructureData::from_parts(empty_parts()).unwrap();
        let schema = discover_atom_site_schema(&none, &BioMmcifWriteParams::default());
        assert_eq!(schema.atom_site_count(), 0);
        // Defaults carry group_pdb = true, so the 18-column no-auth base gains
        // the prepended group_PDB column (19 total).
        assert_eq!(suffixes(&schema).len(), 19);
        assert!(!schema.tags().iter().any(|tag| tag.ends_with("calc_flag")));
        assert!(
            !schema
                .tags()
                .iter()
                .any(|tag| tag.ends_with("pdbx_tls_group_id"))
        );
        assert!(
            !schema
                .tags()
                .iter()
                .any(|tag| tag.ends_with("ccp4_deuterium_fraction"))
        );

        // NoHydrogen alone does not add the calc column; Determined does.
        let no_hydrogen = structure(vec![atom(BioCalcFlag::NoHydrogen, -1)], false);
        let schema = discover_atom_site_schema(&no_hydrogen, &BioMmcifWriteParams::default());
        assert!(!schema.tags().iter().any(|tag| tag.ends_with("calc_flag")));
        let determined = structure(vec![atom(BioCalcFlag::Determined, -1)], false);
        let schema = discover_atom_site_schema(&determined, &BioMmcifWriteParams::default());
        assert!(schema.tags().iter().any(|tag| tag.ends_with("calc_flag")));
        assert!(
            !schema
                .tags()
                .iter()
                .any(|tag| tag.ends_with("pdbx_tls_group_id"))
        );

        // TLS >= 0 adds its column; the deuterium column follows the
        // structure-level source_state flag alone.
        let tls = structure(vec![atom(BioCalcFlag::NotSet, 2)], false);
        let schema = discover_atom_site_schema(&tls, &BioMmcifWriteParams::default());
        assert!(
            schema
                .tags()
                .iter()
                .any(|tag| tag.ends_with("pdbx_tls_group_id"))
        );
        assert!(
            !schema
                .tags()
                .iter()
                .any(|tag| tag.ends_with("ccp4_deuterium_fraction"))
        );
        let deuterium = structure(vec![atom(BioCalcFlag::NotSet, -1)], true);
        let schema = discover_atom_site_schema(&deuterium, &BioMmcifWriteParams::default());
        assert_eq!(
            *schema.tags().last().unwrap(),
            "_atom_site.ccp4_deuterium_fraction"
        );

        // Multiple atoms across models are counted once each; any set flag
        // in any atom adds the column.
        let mut multi_parts = empty_parts();
        multi_parts.models = vec![
            BioModelRow::new(BioRowSpan::new(0, 1).unwrap(), Some(1)),
            BioModelRow::new(BioRowSpan::new(1, 1).unwrap(), Some(2)),
        ];
        multi_parts.chains = vec![
            BioChainRow::new(
                BioModelId::new(0),
                None,
                BioRowSpan::new(0, 1).unwrap(),
                ChainKind::Protein,
                ChainSourceIds::default(),
            ),
            BioChainRow::new(
                BioModelId::new(1),
                None,
                BioRowSpan::new(1, 1).unwrap(),
                ChainKind::Protein,
                ChainSourceIds::default(),
            ),
        ];
        multi_parts.residues = vec![
            BioResidueRow::new(
                BioChainId::new(0),
                BioRowSpan::new(0, 1).unwrap(),
                ResidueName::from_ascii(b"ALA").unwrap(),
                ResidueInfoKind::Aa,
                EntityKind::Polymer,
                None,
                None,
                ResidueSourceIds::default(),
                BioSiftsUnpResidue::default(),
            ),
            BioResidueRow::new(
                BioChainId::new(1),
                BioRowSpan::new(1, 1).unwrap(),
                ResidueName::from_ascii(b"GLY").unwrap(),
                ResidueInfoKind::Aa,
                EntityKind::Polymer,
                None,
                None,
                ResidueSourceIds::default(),
                BioSiftsUnpResidue::default(),
            ),
        ];
        multi_parts.atoms = vec![
            atom_in_residue(0, BioCalcFlag::NotSet, -1),
            atom_in_residue(1, BioCalcFlag::Calculated, 0),
        ];
        multi_parts.coordinates = BioCoordinateBlock::new(vec![[0.0, 0.0, 0.0]; 2]);
        let multi = BioStructureData::from_parts(multi_parts).unwrap();
        let schema = discover_atom_site_schema(&multi, &BioMmcifWriteParams::default());
        assert_eq!(schema.atom_site_count(), 2);
        assert!(schema.tags().iter().any(|tag| tag.ends_with("calc_flag")));
        assert!(
            schema
                .tags()
                .iter()
                .any(|tag| tag.ends_with("pdbx_tls_group_id"))
        );
    }

    fn anisou_atom(residue: u32, anisou: [f64; 6]) -> BioAtomRow {
        BioAtomRow::new(
            BioResidueId::new(residue),
            AtomName::from_ascii(b"CA").unwrap(),
            Element::C,
            None,
            None,
            0,
            BioCalcFlag::NotSet,
            1.0,
            20.0,
            anisou,
            -1,
            0.0,
            AtomSourceIds::default(),
        )
    }

    /// Two models (numbers 1 and 3), one residue per model: model 1 has a
    /// nonzero-tensor atom followed by a zero-tensor atom, model 2 has one
    /// nonzero-tensor atom. Expected columns/values derive from
    /// to_mmcif.cpp:164-187.
    fn a03_structure(model_nums: Option<(i32, i32)>, tensors: [[f64; 6]; 3]) -> BioStructureData {
        let mut parts = empty_parts();
        parts.input_format = BioCoordinateFormat::Mmcif;
        let (models, chains, residues) = match model_nums {
            Some((first, second)) => (
                vec![
                    BioModelRow::new(BioRowSpan::new(0, 1).unwrap(), Some(first)),
                    BioModelRow::new(BioRowSpan::new(1, 1).unwrap(), Some(second)),
                ],
                vec![
                    BioChainRow::new(
                        BioModelId::new(0),
                        None,
                        BioRowSpan::new(0, 1).unwrap(),
                        ChainKind::Protein,
                        ChainSourceIds::default(),
                    ),
                    BioChainRow::new(
                        BioModelId::new(1),
                        None,
                        BioRowSpan::new(1, 1).unwrap(),
                        ChainKind::Protein,
                        ChainSourceIds::default(),
                    ),
                ],
                vec![
                    BioResidueRow::new(
                        BioChainId::new(0),
                        BioRowSpan::new(0, 2).unwrap(),
                        ResidueName::from_ascii(b"ALA").unwrap(),
                        ResidueInfoKind::Aa,
                        EntityKind::Polymer,
                        None,
                        None,
                        ResidueSourceIds::default(),
                        BioSiftsUnpResidue::default(),
                    ),
                    BioResidueRow::new(
                        BioChainId::new(1),
                        BioRowSpan::new(2, 1).unwrap(),
                        ResidueName::from_ascii(b"GLY").unwrap(),
                        ResidueInfoKind::Aa,
                        EntityKind::Polymer,
                        None,
                        None,
                        ResidueSourceIds::default(),
                        BioSiftsUnpResidue::default(),
                    ),
                ],
            ),
            None => (
                vec![BioModelRow::new(BioRowSpan::new(0, 1).unwrap(), Some(1))],
                vec![BioChainRow::new(
                    BioModelId::new(0),
                    None,
                    BioRowSpan::new(0, 1).unwrap(),
                    ChainKind::Protein,
                    ChainSourceIds::default(),
                )],
                vec![BioResidueRow::new(
                    BioChainId::new(0),
                    BioRowSpan::new(0, 3).unwrap(),
                    ResidueName::from_ascii(b"ALA").unwrap(),
                    ResidueInfoKind::Aa,
                    EntityKind::Polymer,
                    None,
                    None,
                    ResidueSourceIds::default(),
                    BioSiftsUnpResidue::default(),
                )],
            ),
        };
        parts.models = models;
        parts.chains = chains;
        parts.residues = residues;
        parts.atoms = vec![
            anisou_atom(0, tensors[0]),
            anisou_atom(0, tensors[1]),
            anisou_atom(if model_nums.is_some() { 1 } else { 0 }, tensors[2]),
        ];
        parts.coordinates = BioCoordinateBlock::new(vec![[1.0, 2.0, 3.0]; 3]);
        BioStructureData::from_parts(parts).unwrap()
    }

    fn anisotrop_loop(block: &crate::cif::CifBlock) -> &crate::cif::CifLoop {
        block
            .items()
            .iter()
            .find_map(|item| match item {
                crate::cif::CifItem::Loop(loop_)
                    if loop_
                        .tags()
                        .first()
                        .is_some_and(|tag| tag.starts_with("_atom_site_anisotrop.")) =>
                {
                    Some(loop_)
                }
                _ => None,
            })
            .expect("_atom_site_anisotrop loop present")
    }

    fn loop_values(loop_: &crate::cif::CifLoop) -> Vec<String> {
        loop_
            .values()
            .iter()
            .map(|value| value.raw().to_string())
            .collect()
    }

    #[test]
    fn bio_pdbscope_a03_multimodel_tail_order_ids_and_tensor_profiles() {
        let data = a03_structure(
            Some((1, 3)),
            [
                [0.05, 0.01, 0.02, 0.03, 0.04, 0.06],
                [0.0; 6],
                [0.11, 0.12, 0.13, 0.14, 0.15, 0.16],
            ],
        );
        let mut document =
            crate::cif::read_cif_document("data_a03\n", "a03", crate::cif::CifCheckLevel::Syntax)
                .unwrap();
        super::add_cif_atoms(
            document.sole_block_mut().unwrap(),
            &data,
            &BioMmcifWriteParams::default(),
        )
        .unwrap();
        let block = document.sole_block().unwrap();
        assert!(
            block
                .find_mmcif_category("_atom_site.")
                .unwrap()
                .is_present()
        );
        let aniso_loop = anisotrop_loop(block);
        // Multimodel: the eight base columns plus pdbx_PDB_model_num, last.
        assert_eq!(
            aniso_loop
                .tags()
                .iter()
                .map(String::as_str)
                .collect::<Vec<_>>(),
            vec![
                "_atom_site_anisotrop.id",
                "_atom_site_anisotrop.type_symbol",
                "_atom_site_anisotrop.U[1][1]",
                "_atom_site_anisotrop.U[2][2]",
                "_atom_site_anisotrop.U[3][3]",
                "_atom_site_anisotrop.U[1][2]",
                "_atom_site_anisotrop.U[1][3]",
                "_atom_site_anisotrop.U[2][3]",
                "_atom_site_anisotrop.pdbx_PDB_model_num",
            ]
        );
        // Rows follow emission order: serial 1 (model 1) then serial 3
        // (model 3); the zero-tensor serial 2 is suppressed. Tensor values
        // use the float profile (SMat33<float> to_str).
        assert_eq!(
            loop_values(aniso_loop),
            vec![
                "1", "C", "0.05", "0.01", "0.02", "0.03", "0.04", "0.06", "1", "3", "C", "0.11",
                "0.12", "0.13", "0.14", "0.15", "0.16", "3",
            ]
        );
    }

    #[test]
    fn bio_pdbscope_a03_single_model_omits_model_num_column() {
        let data = a03_structure(
            None,
            [
                [0.05, 0.01, 0.02, 0.03, 0.04, 0.06],
                [0.0; 6],
                [0.07, 0.08, 0.09, 0.1, 0.11, 0.12],
            ],
        );
        let mut document =
            crate::cif::read_cif_document("data_a03\n", "a03", crate::cif::CifCheckLevel::Syntax)
                .unwrap();
        super::add_cif_atoms(
            document.sole_block_mut().unwrap(),
            &data,
            &BioMmcifWriteParams::default(),
        )
        .unwrap();
        let block = document.sole_block().unwrap();
        let aniso_loop = anisotrop_loop(block);
        assert_eq!(aniso_loop.tags().len(), 8);
        assert!(
            aniso_loop
                .tags()
                .iter()
                .all(|tag| !tag.ends_with("pdbx_PDB_model_num"))
        );
        assert_eq!(
            loop_values(aniso_loop),
            vec![
                "1", "C", "0.05", "0.01", "0.02", "0.03", "0.04", "0.06", "3", "C", "0.07", "0.08",
                "0.09", "0.1", "0.11", "0.12",
            ]
        );
    }

    #[test]
    fn bio_pdbscope_a03_all_zero_tensors_erase_existing_category() {
        // A pre-existing _atom_site_anisotrop pair proves the source's erase
        // branch actually removes the category from the parsed block.
        let data = a03_structure(Some((1, 3)), [[0.0; 6], [0.0; 6], [0.0; 6]]);
        let mut document = crate::cif::read_cif_document(
            "data_a03\n_atom_site_anisotrop.id 1\n",
            "a03",
            crate::cif::CifCheckLevel::Syntax,
        )
        .unwrap();
        assert!(
            document
                .sole_block()
                .unwrap()
                .find_mmcif_category("_atom_site_anisotrop.")
                .unwrap()
                .is_present()
        );
        super::add_cif_atoms(
            document.sole_block_mut().unwrap(),
            &data,
            &BioMmcifWriteParams::default(),
        )
        .unwrap();
        let block = document.sole_block().unwrap();
        assert!(
            !block
                .find_mmcif_category("_atom_site_anisotrop.")
                .unwrap()
                .is_present(),
            "empty aniso selection must erase the category"
        );
        assert!(
            block
                .find_mmcif_category("_atom_site.")
                .unwrap()
                .is_present(),
            "the atom loop is still written"
        );
    }
}
