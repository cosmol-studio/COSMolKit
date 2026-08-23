//! Gemmi-aligned `BioStructure` mmCIF writer.

use super::cif::{CifBlock, CifDocument, CifEntry, CifToken, CifWriteOptions, quote_cif_value};
use super::sprintf::{format_fixed, format_general};
use crate::bio::{
    AltLocLabel, AtomName, BioAssemblySpecialKind, BioAsu, BioAtomAddress, BioCalcFlag,
    BioConnectionType, BioStructure, BioTransform, CrystalCell, EntityKind, PdbChainId,
    PolymerKind, ResidueRow,
};
use std::collections::BTreeSet;

#[derive(Debug, thiserror::Error)]
pub enum BioWriteError {
    #[error("BioStructure invariant violation while writing: {0}")]
    Invariant(&'static str),
    #[error("BioStructure contains non-UTF-8 bytes in {field}")]
    InvalidText { field: &'static str },
    #[error(transparent)]
    Io(#[from] std::io::Error),
}

// BEGIN GEMMI CPP TYPE gemmi::MmcifOutputGroups
// Gemmi✔️✔️: bool atoms:1;
// Gemmi✔️✔️: bool block_name:1;
// Gemmi✔️✔️: bool entry:1;
// Gemmi✔️✔️: bool database_status:1;
// Gemmi✔️✔️: bool author:1;
// Gemmi✔️✔️: bool cell:1;
// Gemmi✔️✔️: bool symmetry:1;
// Gemmi✔️✔️: bool entity:1;
// Gemmi✔️✔️: bool entity_poly:1;
// Gemmi✔️✔️: bool struct_ref:1;
// Gemmi✔️✔️: bool chem_comp:1;
// Gemmi✔️✔️: bool exptl:1;
// Gemmi✔️✔️: bool diffrn:1;
// Gemmi✔️✔️: bool reflns:1;
// Gemmi✔️✔️: bool refine:1;
// Gemmi✔️✔️: bool title_keywords:1;
// Gemmi✔️✔️: bool ncs:1;
// Gemmi✔️✔️: bool struct_asym:1;
// Gemmi✔️✔️: bool origx:1;
// Gemmi✔️✔️: bool struct_conf:1;
// Gemmi✔️✔️: bool struct_sheet:1;
// Gemmi✔️✔️: bool struct_biol:1;
// Gemmi✔️✔️: bool assembly:1;
// Gemmi✔️✔️: bool conn:1;
// Gemmi✔️✔️: bool cis:1;
// Gemmi✔️✔️: bool modres:1;
// Gemmi✔️✔️: bool scale:1;
// Gemmi✔️✔️: bool atom_type:1;
// Gemmi✔️✔️: bool entity_poly_seq:1;
// Gemmi✔️✔️: bool tls:1;
// Gemmi✔️✔️: bool software:1;
// Gemmi✔️✔️: bool group_pdb:1;
// Gemmi✔️✔️: bool auth_all:1;
// END GEMMI CPP TYPE
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MmcifOutputGroups {
    pub atoms: bool,
    pub block_name: bool,
    pub entry: bool,
    pub database_status: bool,
    pub author: bool,
    pub cell: bool,
    pub symmetry: bool,
    pub entity: bool,
    pub entity_poly: bool,
    pub struct_ref: bool,
    pub chem_comp: bool,
    pub exptl: bool,
    pub diffrn: bool,
    pub reflns: bool,
    pub refine: bool,
    pub title_keywords: bool,
    pub ncs: bool,
    pub struct_asym: bool,
    pub origx: bool,
    pub struct_conf: bool,
    pub struct_sheet: bool,
    pub struct_biol: bool,
    pub assembly: bool,
    pub conn: bool,
    pub cis: bool,
    pub modres: bool,
    pub scale: bool,
    pub atom_type: bool,
    pub entity_poly_seq: bool,
    pub tls: bool,
    pub software: bool,
    pub group_pdb: bool,
    pub auth_all: bool,
}

impl MmcifOutputGroups {
    // BEGIN GEMMI CPP FUNCTION gemmi::MmcifOutputGroups::MmcifOutputGroups
    // Gemmi✔️✔️: explicit MmcifOutputGroups(bool all)
    // Gemmi✔️✔️:   : atoms(all), block_name(all), entry(all), database_status(all),
    // Gemmi✔️✔️:     author(all), cell(all), symmetry(all), entity(all), entity_poly(all),
    // Gemmi✔️✔️:     struct_ref(all), chem_comp(all), exptl(all), diffrn(all),
    // Gemmi✔️✔️:     reflns(all), refine(all), title_keywords(all), ncs(all),
    // Gemmi✔️✔️:     struct_asym(all), origx(all), struct_conf(all), struct_sheet(all),
    // Gemmi✔️✔️:     struct_biol(all), assembly(all), conn(all), cis(all), modres(all),
    // Gemmi✔️✔️:     scale(all), atom_type(all), entity_poly_seq(all), tls(all),
    // Gemmi✔️✔️:     software(all), group_pdb(all), auth_all(false) {}
    // END GEMMI CPP FUNCTION
    #[must_use]
    pub const fn all(all: bool) -> Self {
        Self {
            atoms: all,
            block_name: all,
            entry: all,
            database_status: all,
            author: all,
            cell: all,
            symmetry: all,
            entity: all,
            entity_poly: all,
            struct_ref: all,
            chem_comp: all,
            exptl: all,
            diffrn: all,
            reflns: all,
            refine: all,
            title_keywords: all,
            ncs: all,
            struct_asym: all,
            origx: all,
            struct_conf: all,
            struct_sheet: all,
            struct_biol: all,
            assembly: all,
            conn: all,
            cis: all,
            modres: all,
            scale: all,
            atom_type: all,
            entity_poly_seq: all,
            tls: all,
            software: all,
            group_pdb: all,
            auth_all: false,
        }
    }
}

impl Default for MmcifOutputGroups {
    fn default() -> Self {
        Self::all(true)
    }
}

/// Options for canonical Gemmi-aligned mmCIF output.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct MmcifWriteOptions {
    pub groups: MmcifOutputGroups,
    pub prefer_pairs: bool,
    pub compact: bool,
    pub misuse_hash: bool,
    pub align_pairs: u16,
    pub align_loops: u16,
}

impl MmcifWriteOptions {
    pub(super) const fn cif_options(self) -> CifWriteOptions {
        CifWriteOptions {
            prefer_pairs: self.prefer_pairs,
            compact: self.compact,
            misuse_hash: self.misuse_hash,
            align_pairs: self.align_pairs,
            align_loops: self.align_loops,
        }
    }
}

// BEGIN GEMMI CPP FUNCTION gemmi::to_str(float)
// Gemmi✔️✔️: char buf[16];
// Gemmi✔️✔️: int len = sprintf_z(buf, "%.6g", d);
// Gemmi✔️✔️: return std::string(buf, len > 0 ? len : 0);
// END GEMMI CPP FUNCTION
pub(super) fn gemmi_to_string_f32(value: f32) -> String {
    format_general(value as f64, 6)
}

// BEGIN GEMMI CPP FUNCTION gemmi::to_str(double)
// Gemmi✔️✔️: char buf[24];
// Gemmi✔️✔️: int len = sprintf_z(buf, "%.9g", d);
// Gemmi✔️✔️: return std::string(buf, len > 0 ? len : 0);
// END GEMMI CPP FUNCTION
pub(super) fn gemmi_to_string_f64(value: f64) -> String {
    format_general(value, 9)
}

// BEGIN GEMMI CPP FUNCTION gemmi::number_or_dot / gemmi::number_or_qmark
// Gemmi✔️✔️: return std::isnan(d) ? "." : to_str(d);
// Gemmi✔️✔️: return std::isnan(d) ? "?" : to_str(d);
// END GEMMI CPP FUNCTION
pub(super) fn number_or_dot(value: Option<f64>) -> String {
    value
        .filter(|value| !value.is_nan())
        .map_or_else(|| ".".to_string(), gemmi_to_string_f64)
}

pub(super) fn number_or_qmark(value: Option<f64>) -> String {
    value
        .filter(|value| !value.is_nan())
        .map_or_else(|| "?".to_string(), gemmi_to_string_f64)
}

// BEGIN GEMMI CPP FUNCTION gemmi::int_or_dot / gemmi::int_or_qmark
// Gemmi✔️✔️: return n == -1 ? "." : std::to_string(n);
// Gemmi✔️✔️: return n == -1 ? "?" : std::to_string(n);
// END GEMMI CPP FUNCTION
pub(super) fn int_or_dot(value: Option<i32>) -> String {
    value.map_or_else(|| ".".to_string(), |value| value.to_string())
}

pub(super) fn int_or_qmark(value: Option<i32>) -> String {
    value.map_or_else(|| "?".to_string(), |value| value.to_string())
}

// BEGIN GEMMI CPP FUNCTION gemmi::is_valid_block_name
// Gemmi✔️✔️: return !name.empty() &&
// Gemmi✔️✔️:        std::all_of(name.begin(), name.end(), [](char c){ return c >= '!' && c <= '~'; });
// END GEMMI CPP FUNCTION
pub(super) fn is_valid_block_name(name: &str) -> bool {
    !name.is_empty() && name.bytes().all(|byte| (b'!'..=b'~').contains(&byte))
}

// BEGIN GEMMI CPP FUNCTION gemmi::write_cell_parameters
// Gemmi✔️✔️: span.set_pair("_cell.length_a",    to_str(cell.a));
// Gemmi✔️✔️: span.set_pair("_cell.length_b",    to_str(cell.b));
// Gemmi✔️✔️: span.set_pair("_cell.length_c",    to_str(cell.c));
// Gemmi✔️✔️: span.set_pair("_cell.angle_alpha", to_str(cell.alpha));
// Gemmi✔️✔️: span.set_pair("_cell.angle_beta",  to_str(cell.beta));
// Gemmi✔️✔️: span.set_pair("_cell.angle_gamma", to_str(cell.gamma));
// END GEMMI CPP FUNCTION
pub(super) fn write_cell_parameters(cell: CrystalCell, block: &mut CifBlock) {
    for (tag, value) in [
        ("_cell.length_a", cell.a),
        ("_cell.length_b", cell.b),
        ("_cell.length_c", cell.c),
        ("_cell.angle_alpha", cell.alpha),
        ("_cell.angle_beta", cell.beta),
        ("_cell.angle_gamma", cell.gamma),
    ] {
        block.set_pair_in_category("_cell.", tag, gemmi_to_string_f64(value));
    }
}

fn atom_name_text(name: AtomName) -> Result<String, BioWriteError> {
    let bytes = name.0;
    let start = bytes
        .iter()
        .position(|byte| *byte != b' ')
        .unwrap_or(bytes.len());
    let end = bytes
        .iter()
        .rposition(|byte| *byte != b' ')
        .map_or(start, |index| index + 1);
    std::str::from_utf8(&bytes[start..end])
        .map(str::to_string)
        .map_err(|_| BioWriteError::InvalidText { field: "atom name" })
}

fn pdb_chain_text<'a>(id: &'a PdbChainId, field: &'static str) -> Result<&'a str, BioWriteError> {
    let length = usize::from(id.1);
    if length > id.0.len() {
        return Err(BioWriteError::Invariant(
            "PDB chain identifier length exceeds its storage",
        ));
    }
    std::str::from_utf8(&id.0[..length]).map_err(|_| BioWriteError::InvalidText { field })
}

fn residue_name_text(residue: &ResidueRow) -> Result<&str, BioWriteError> {
    let length = usize::from(residue.name.1);
    if length > residue.name.0.len() {
        return Err(BioWriteError::Invariant(
            "residue name length exceeds its storage",
        ));
    }
    std::str::from_utf8(&residue.name.0[..length]).map_err(|_| BioWriteError::InvalidText {
        field: "residue name",
    })
}

// BEGIN GEMMI CPP FUNCTION gemmi::use_hetatm
// Gemmi✔️✔️: bool use_hetatm(const Residue& res) {
// Gemmi✔️✔️:   if (res.het_flag == 'H')
// Gemmi✔️✔️:     return true;
// Gemmi✔️✔️:   if (res.het_flag == 'A')
// Gemmi✔️✔️:     return false;
// Gemmi✔️✔️:   if (res.entity_type == EntityType::Branched ||
// Gemmi✔️✔️:       res.entity_type == EntityType::NonPolymer ||
// Gemmi✔️✔️:       res.entity_type == EntityType::Water)
// Gemmi✔️✔️:     return true;
// Gemmi✔️✔️:   return !find_tabulated_residue(res.name).is_standard();
// Gemmi✔️✔️: }
// END GEMMI CPP FUNCTION
fn use_hetatm(residue: &ResidueRow) -> Result<bool, BioWriteError> {
    if residue.het_flag == Some('H') {
        return Ok(true);
    }
    if residue.het_flag == Some('A') {
        return Ok(false);
    }
    if matches!(
        residue.entity_kind,
        EntityKind::Branched | EntityKind::NonPolymer | EntityKind::Water
    ) {
        return Ok(true);
    }
    Ok(!crate::bio::resinfo::find_tabulated_residue(residue_name_text(residue)?).is_standard())
}

fn push_cif_value(loop_: &mut super::cif::CifLoop, value: String) {
    loop_.values.push(CifToken {
        value,
        line_number: 0,
    });
}

fn calc_flag_text(flag: BioCalcFlag) -> &'static str {
    match flag {
        BioCalcFlag::NotSet | BioCalcFlag::NoHydrogen => ".",
        BioCalcFlag::Determined => "d",
        BioCalcFlag::Calculated => "c",
        BioCalcFlag::Dummy => "dum",
    }
}

fn entity_id_for_residue(structure: &BioStructure, residue: &ResidueRow) -> String {
    // Gemmi✔️✔️: std::string entity_id;
    // Gemmi✔️✔️: if (const Entity* ent = gemmi::find_entity_of_subchain(res.subchain, st.entities))
    // Gemmi✔️✔️:   entity_id = cif::quote(ent->name);
    // Gemmi✔️✔️: else
    // Gemmi✔️✔️:   entity_id = string_or_dot(res.entity_id);
    if let Some(subchain) = residue.source.subchain_id
        && let Some(entity) = structure
            .entities
            .iter()
            .find(|entity| entity.subchains.contains(&subchain))
    {
        return quote_cif_value(&entity.source.source_entity_id);
    }
    residue
        .source
        .label_entity_id
        .and_then(|entity_id| structure.entities.get(entity_id.index() as usize))
        .map_or_else(
            || ".".to_string(),
            |entity| quote_cif_value(&entity.source.source_entity_id),
        )
}

// BEGIN GEMMI CPP FUNCTION gemmi::add_cif_atoms
// Gemmi✔️✔️: void add_cif_atoms(const Structure& st, cif::Block& block,
// Gemmi✔️✔️:                    bool use_group_pdb, bool auth_all) {
// Gemmi✔️✔️:   // atom list
// Gemmi✔️✔️:   cif::Loop& atom_loop = block.init_mmcif_loop("_atom_site.", {
// Gemmi✔️✔️:       "id",
// Gemmi✔️✔️:       "type_symbol",
// Gemmi✔️✔️:       "label_atom_id",
// Gemmi✔️✔️:       "label_alt_id",
// Gemmi✔️✔️:       "label_comp_id",
// Gemmi✔️✔️:       "label_asym_id",
// Gemmi✔️✔️:       "label_entity_id",
// Gemmi✔️✔️:       "label_seq_id",
// Gemmi✔️✔️:       "pdbx_PDB_ins_code",
// Gemmi✔️✔️:       "Cartn_x",
// Gemmi✔️✔️:       "Cartn_y",
// Gemmi✔️✔️:       "Cartn_z",
// Gemmi✔️✔️:       "occupancy",
// Gemmi✔️✔️:       "B_iso_or_equiv",
// Gemmi✔️✔️:       "pdbx_formal_charge",
// Gemmi✔️✔️:       "auth_atom_id",  // optional (tags[15] is removed if !auth_all)
// Gemmi✔️✔️:       "auth_comp_id",  // optional (tags[16] is removed if !auth_all)
// Gemmi✔️✔️:       "auth_seq_id",
// Gemmi✔️✔️:       "auth_asym_id",
// Gemmi✔️✔️:       "pdbx_PDB_model_num"});
// Gemmi✔️✔️:   if (!auth_all)
// Gemmi✔️✔️:     atom_loop.tags.erase(atom_loop.tags.begin() + 15, atom_loop.tags.begin() + 17);
// Gemmi✔️✔️:   if (use_group_pdb)
// Gemmi✔️✔️:     atom_loop.tags.emplace(atom_loop.tags.begin(), "_atom_site.group_PDB");
// Gemmi✔️✔️:   bool has_calc_flag = false;
// Gemmi✔️✔️:   bool has_tls_group_id = false;
// Gemmi✔️✔️:   size_t atom_site_count = 0;
// Gemmi✔️✔️:   for (const Model& model : st.models)
// Gemmi✔️✔️:     for (const Chain& chain : model.chains)
// Gemmi✔️✔️:       for (const Residue& res : chain.residues)
// Gemmi✔️✔️:         for (const Atom& atom : res.atoms) {
// Gemmi✔️✔️:           ++atom_site_count;
// Gemmi✔️✔️:           if (atom.calc_flag != CalcFlag::NotSet &&
// Gemmi✔️✔️:               atom.calc_flag != CalcFlag::NoHydrogen)
// Gemmi✔️✔️:             has_calc_flag = true;
// Gemmi✔️✔️:           if (atom.tls_group_id >= 0)
// Gemmi✔️✔️:             has_tls_group_id = true;
// Gemmi✔️✔️:         }
// Gemmi✔️✔️:   if (has_calc_flag)
// Gemmi✔️✔️:     atom_loop.tags.emplace_back("_atom_site.calc_flag");
// Gemmi✔️✔️:   if (has_tls_group_id)
// Gemmi✔️✔️:     atom_loop.tags.emplace_back("_atom_site.pdbx_tls_group_id");
// Gemmi✔️✔️:   if (st.has_d_fraction)
// Gemmi✔️✔️:     atom_loop.tags.emplace_back("_atom_site.ccp4_deuterium_fraction");
// Gemmi✔️✔️:
// Gemmi✔️✔️:   std::vector<std::string>& vv = atom_loop.values;
// Gemmi✔️✔️:   vv.reserve(atom_site_count * atom_loop.tags.size());
// Gemmi✔️✔️:   std::vector<std::tuple<int, int, const Atom*>> aniso;
// Gemmi✔️✔️:   int serial = 0;
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
// Gemmi✔️✔️: }
// END GEMMI CPP FUNCTION
pub(super) fn add_cif_atoms(
    structure: &BioStructure,
    block: &mut CifBlock,
    use_group_pdb: bool,
    auth_all: bool,
) -> Result<(), BioWriteError> {
    let atom_site_suffixes = [
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
    ];
    let has_calc_flag = structure.atoms.iter().any(|atom| {
        !matches!(
            atom.calc_flag,
            BioCalcFlag::NotSet | BioCalcFlag::NoHydrogen
        )
    });
    let has_tls_group_id = structure
        .atoms
        .iter()
        .any(|atom| atom.tls_group_id.is_some_and(|group_id| group_id >= 0));
    let atom_site_count = structure.atoms.len();

    let mut anisotropic_atoms = Vec::<(i32, i32, usize)>::new();
    let atom_loop = block.init_mmcif_loop("_atom_site.", &atom_site_suffixes);
    if !auth_all {
        atom_loop.tags.drain(15..17);
    }
    if use_group_pdb {
        atom_loop.tags.insert(0, "_atom_site.group_PDB".to_string());
    }
    if has_calc_flag {
        atom_loop.tags.push("_atom_site.calc_flag".to_string());
    }
    if has_tls_group_id {
        atom_loop
            .tags
            .push("_atom_site.pdbx_tls_group_id".to_string());
    }
    if structure.has_d_fraction {
        atom_loop
            .tags
            .push("_atom_site.ccp4_deuterium_fraction".to_string());
    }
    let value_capacity = atom_site_count
        .checked_mul(atom_loop.tags.len())
        .ok_or(BioWriteError::Invariant("atom-site value count overflow"))?;
    atom_loop.values.reserve(value_capacity);

    let mut serial = 0_i32;
    for model in &structure.models {
        let model_number = model.source_model_number.unwrap_or(0);
        let chain_start = model.chain_span.start as usize;
        let chain_end = model.chain_span.end() as usize;
        let chains =
            structure
                .chains
                .get(chain_start..chain_end)
                .ok_or(BioWriteError::Invariant(
                    "model chain span is out of bounds",
                ))?;
        for chain in chains {
            let auth_chain = match chain.source.auth_chain_id {
                Some(ref id) => pdb_chain_text(id, "author chain identifier")?.to_string(),
                None => String::new(),
            };
            let residue_start = chain.residue_span.start as usize;
            let residue_end = chain.residue_span.end() as usize;
            let residues = structure.residues.get(residue_start..residue_end).ok_or(
                BioWriteError::Invariant("chain residue span is out of bounds"),
            )?;
            for residue in residues {
                let as_het = use_hetatm(residue)?;
                let residue_name = residue_name_text(residue)?.to_string();
                let label_seq_id = residue
                    .source
                    .label_seq_id
                    .map_or_else(|| ".".to_string(), |number| number.to_string());
                let auth_seq_id = residue
                    .source
                    .seq_id
                    .map_or_else(|| "?".to_string(), |seq_id| seq_id.seq_num.to_string());
                let insertion_code = residue
                    .source
                    .seq_id
                    .and_then(|seq_id| seq_id.ins_code)
                    .map_or_else(|| "?".to_string(), |code| char::from(code).to_string());
                let label_asym_id = match residue.source.subchain_id {
                    Some(ref id) => quote_cif_value(pdb_chain_text(id, "label asym identifier")?),
                    None => ".".to_string(),
                };
                let entity_id = entity_id_for_residue(structure, residue);
                let atom_start = residue.atom_span.start as usize;
                let atom_end = residue.atom_span.end() as usize;
                let atoms =
                    structure
                        .atoms
                        .get(atom_start..atom_end)
                        .ok_or(BioWriteError::Invariant(
                            "residue atom span is out of bounds",
                        ))?;
                for (atom_offset, atom) in atoms.iter().enumerate() {
                    let atom_index = atom_start + atom_offset;
                    let position = structure.coordinates.positions.get(atom_index).ok_or(
                        BioWriteError::Invariant("atom coordinate index is out of bounds"),
                    )?;
                    serial = serial
                        .checked_add(1)
                        .ok_or(BioWriteError::Invariant("atom-site serial overflow"))?;
                    if use_group_pdb {
                        push_cif_value(
                            atom_loop,
                            if as_het { "HETATM" } else { "ATOM" }.to_string(),
                        );
                    }
                    push_cif_value(atom_loop, serial.to_string());
                    push_cif_value(atom_loop, atom.element.symbol().to_ascii_uppercase());
                    let quoted_atom_name = quote_cif_value(&atom_name_text(atom.name)?);
                    push_cif_value(atom_loop, quoted_atom_name.clone());
                    push_cif_value(
                        atom_loop,
                        atom.altloc
                            .map_or_else(|| ".".to_string(), |alt| char::from(alt.0).to_string()),
                    );
                    let quoted_residue_name = quote_cif_value(&residue_name);
                    push_cif_value(atom_loop, quoted_residue_name.clone());
                    push_cif_value(atom_loop, label_asym_id.clone());
                    push_cif_value(atom_loop, entity_id.clone());
                    push_cif_value(atom_loop, label_seq_id.clone());
                    push_cif_value(atom_loop, insertion_code.clone());
                    push_cif_value(atom_loop, gemmi_to_string_f64(position[0]));
                    push_cif_value(atom_loop, gemmi_to_string_f64(position[1]));
                    push_cif_value(atom_loop, gemmi_to_string_f64(position[2]));
                    push_cif_value(
                        atom_loop,
                        gemmi_to_string_f32(atom.occupancy.unwrap_or(f32::NAN)),
                    );
                    push_cif_value(
                        atom_loop,
                        gemmi_to_string_f32(atom.b_iso.unwrap_or(f32::NAN)),
                    );
                    push_cif_value(
                        atom_loop,
                        atom.formal_charge
                            .filter(|charge| *charge != 0)
                            .map_or_else(|| "?".to_string(), |charge| charge.to_string()),
                    );
                    if auth_all {
                        push_cif_value(atom_loop, quoted_atom_name);
                        push_cif_value(atom_loop, quoted_residue_name);
                    }
                    push_cif_value(atom_loop, auth_seq_id.clone());
                    push_cif_value(atom_loop, quote_cif_value(&auth_chain));
                    push_cif_value(atom_loop, model_number.to_string());
                    if has_calc_flag {
                        push_cif_value(atom_loop, calc_flag_text(atom.calc_flag).to_string());
                    }
                    if has_tls_group_id {
                        push_cif_value(
                            atom_loop,
                            atom.tls_group_id
                                .filter(|group_id| *group_id != -1)
                                .map_or_else(|| "?".to_string(), |group_id| group_id.to_string()),
                        );
                    }
                    if structure.has_d_fraction {
                        push_cif_value(
                            atom_loop,
                            gemmi_to_string_f32(atom.fraction.unwrap_or(0.0)),
                        );
                    }
                    if atom
                        .anisou
                        .is_some_and(|aniso| aniso[0] + aniso[1] + aniso[2] != 0.0)
                    {
                        anisotropic_atoms.push((serial, model_number, atom_index));
                    }
                }
            }
        }
    }

    if anisotropic_atoms.is_empty() {
        block.erase_mmcif_category("_atom_site_anisotrop.");
    } else {
        let aniso_suffixes = [
            "id",
            "type_symbol",
            "U[1][1]",
            "U[2][2]",
            "U[3][3]",
            "U[1][2]",
            "U[1][3]",
            "U[2][3]",
        ];
        let aniso_loop = block.init_mmcif_loop("_atom_site_anisotrop.", &aniso_suffixes);
        if structure.models.len() > 1 {
            aniso_loop
                .tags
                .push("_atom_site_anisotrop.pdbx_PDB_model_num".to_string());
        }
        let aniso_capacity = aniso_loop
            .tags
            .len()
            .checked_mul(anisotropic_atoms.len())
            .ok_or(BioWriteError::Invariant(
                "atom-site anisotropic value count overflow",
            ))?;
        aniso_loop.values.reserve(aniso_capacity);
        for (atom_serial, model_number, atom_index) in anisotropic_atoms {
            let atom = &structure.atoms[atom_index];
            let aniso = atom
                .anisou
                .expect("anisotropic index only records Some values");
            push_cif_value(aniso_loop, atom_serial.to_string());
            push_cif_value(aniso_loop, atom.element.symbol().to_ascii_uppercase());
            for value in aniso {
                push_cif_value(aniso_loop, gemmi_to_string_f32(value));
            }
            if structure.models.len() > 1 {
                push_cif_value(aniso_loop, model_number.to_string());
            }
        }
    }
    Ok(())
}

// BEGIN GEMMI CPP FUNCTION gemmi::xmeric_to_number
// Gemmi✔️✔️: // the names are: monomeric, dimeric, ...meric, 21-meric, 22-meric, ...
// Gemmi✔️✔️: int xmeric_to_number(const std::string& oligomeric) {
// Gemmi✔️✔️:   static const char names[20][10] = {
// Gemmi✔️✔️:     "mono", "di", "tri", "tetra", "penta",
// Gemmi✔️✔️:     "hexa", "hepta", "octa", "nona", "deca",
// Gemmi✔️✔️:     "undeca", "dodeca", "trideca", "tetradeca", "pentadeca",
// Gemmi✔️✔️:     "hexadeca", "heptadeca", "octadeca", "nonadeca", "eicosa"
// Gemmi✔️✔️:   };
// Gemmi✔️✔️:   size_t len = oligomeric.length();
// Gemmi✔️✔️:   const char* p = oligomeric.c_str();
// Gemmi✔️✔️:   for (int i = 0; i != 20; ++i)
// Gemmi✔️✔️:     if (len == std::strlen(names[i]) + 5 && strncmp(p, names[i], len-5) == 0)
// Gemmi✔️✔️:       return i + 1;
// Gemmi✔️✔️:   return no_sign_atoi(p);
// Gemmi✔️✔️: }
// END GEMMI CPP FUNCTION
fn xmeric_to_number(oligomeric: &str) -> i32 {
    const NAMES: [&str; 20] = [
        "mono",
        "di",
        "tri",
        "tetra",
        "penta",
        "hexa",
        "hepta",
        "octa",
        "nona",
        "deca",
        "undeca",
        "dodeca",
        "trideca",
        "tetradeca",
        "pentadeca",
        "hexadeca",
        "heptadeca",
        "octadeca",
        "nonadeca",
        "eicosa",
    ];
    for (index, name) in NAMES.iter().enumerate() {
        if oligomeric.len() == name.len() + 5 && oligomeric.starts_with(name) {
            return index as i32 + 1;
        }
    }
    let digits = oligomeric
        .trim_start_matches(|character: char| character.is_ascii_whitespace())
        .bytes()
        .take_while(u8::is_ascii_digit);
    let mut number = 0_i32;
    for digit in digits {
        number = number
            .wrapping_mul(10)
            .wrapping_add(i32::from(digit - b'0'));
    }
    number
}

// BEGIN GEMMI CPP FUNCTION gemmi::string_append_sep
// Gemmi✔️✔️: if (!str.empty())
// Gemmi✔️✔️:   str += sep;
// Gemmi✔️✔️: str += item;
// END GEMMI CPP FUNCTION
fn append_with_comma(output: &mut String, item: &str) {
    if !output.is_empty() {
        output.push(',');
    }
    output.push_str(item);
}

// BEGIN GEMMI CPP FUNCTION gemmi::Transform::approx
// Gemmi✔️✔️: return mat.approx(o.mat, epsilon) && vec.approx(o.vec, epsilon);
// END GEMMI CPP FUNCTION
// BEGIN GEMMI CPP FUNCTION gemmi::Mat33::approx
// Gemmi✔️✔️: for (int i = 0; i < 3; ++i)
// Gemmi✔️✔️:   for (int j = 0; j < 3; ++j)
// Gemmi✔️✔️:     if (std::fabs(a[i][j] - other.a[i][j]) > epsilon)
// Gemmi✔️✔️:       return false;
// Gemmi✔️✔️: return true;
// END GEMMI CPP FUNCTION
// BEGIN GEMMI CPP FUNCTION gemmi::Vec3_::approx
// Gemmi✔️✔️: return std::fabs(x - o.x) <= epsilon &&
// Gemmi✔️✔️:        std::fabs(y - o.y) <= epsilon &&
// Gemmi✔️✔️:        std::fabs(z - o.z) <= epsilon;
// END GEMMI CPP FUNCTION
fn transform_approx(left: &BioTransform, right: &BioTransform, epsilon: f64) -> bool {
    for row in 0..3 {
        for column in 0..3 {
            if (left.mat[row][column] - right.mat[row][column]).abs() > epsilon {
                return false;
            }
        }
    }
    (left.vec[0] - right.vec[0]).abs() <= epsilon
        && (left.vec[1] - right.vec[1]).abs() <= epsilon
        && (left.vec[2] - right.vec[2]).abs() <= epsilon
}

// BEGIN GEMMI CPP FUNCTION gemmi::Transform::is_identity
// Gemmi✔️✔️: return mat.is_identity() && vec.x == 0. && vec.y == 0. && vec.z == 0.;
// END GEMMI CPP FUNCTION
// BEGIN GEMMI CPP FUNCTION gemmi::Mat33::is_identity
// Gemmi✔️✔️: return a[0][0] == 1 && a[0][1] == 0 && a[0][2] == 0 &&
// Gemmi✔️✔️:        a[1][0] == 0 && a[1][1] == 1 && a[1][2] == 0 &&
// Gemmi✔️✔️:        a[2][0] == 0 && a[2][1] == 0 && a[2][2] == 1;
// END GEMMI CPP FUNCTION
fn transform_is_exact_identity(transform: &BioTransform) -> bool {
    transform.mat == [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        && transform.vec == [0.0, 0.0, 0.0]
}

fn active_loop_index(block: &CifBlock, tag: &str) -> Result<usize, BioWriteError> {
    block
        .entries
        .iter()
        .find_map(|entry| match entry {
            CifEntry::Loop(index) if block.loops[*index].tags.iter().any(|item| item == tag) => {
                Some(*index)
            }
            CifEntry::Pair(_) | CifEntry::Loop(_) | CifEntry::Erased => None,
        })
        .ok_or(BioWriteError::Invariant(
            "initialized CIF loop is not active in the block",
        ))
}

fn chain_name_text(
    chain: &crate::bio::ChainRow,
    field: &'static str,
) -> Result<String, BioWriteError> {
    match chain.source.auth_chain_id.or(chain.source.label_asym_id) {
        Some(ref id) => Ok(pdb_chain_text(id, field)?.to_string()),
        None => Ok(String::new()),
    }
}

fn generator_subchains(
    structure: &BioStructure,
    generator: &crate::bio::BioAssemblyGenerator,
) -> Result<String, BioWriteError> {
    let mut output = String::new();
    for name in &generator.subchains {
        append_with_comma(&mut output, name);
    }
    if !output.is_empty() {
        return Ok(output);
    }
    let Some(model) = structure.models.first() else {
        return Ok(output);
    };
    let chain_start = model.chain_span.start as usize;
    let chain_end = model.chain_span.end() as usize;
    let chains = structure
        .chains
        .get(chain_start..chain_end)
        .ok_or(BioWriteError::Invariant(
            "model chain span is out of bounds",
        ))?;
    for chain in chains {
        let chain_name = chain_name_text(chain, "assembly chain identifier")?;
        if !generator.chains.iter().any(|name| *name == chain_name) {
            continue;
        }
        let residue_start = chain.residue_span.start as usize;
        let residue_end = chain.residue_span.end() as usize;
        let residues =
            structure
                .residues
                .get(residue_start..residue_end)
                .ok_or(BioWriteError::Invariant(
                    "chain residue span is out of bounds",
                ))?;
        let mut previous = None;
        for residue in residues {
            if previous == Some(residue.source.subchain_id) {
                continue;
            }
            previous = Some(residue.source.subchain_id);
            let subchain = match residue.source.subchain_id {
                Some(ref id) => pdb_chain_text(id, "assembly subchain identifier")?,
                None => "",
            };
            append_with_comma(&mut output, subchain);
        }
    }
    Ok(output)
}

// BEGIN GEMMI CPP FUNCTION gemmi::write_assemblies
// Gemmi✔️✔️: void write_assemblies(const Structure& st, cif::Block& block) {
// Gemmi✔️✔️:   block.items.reserve(block.items.size() + 4); // avoid re-allocation
// Gemmi✔️✔️:   cif::Loop& a_loop = block.init_mmcif_loop("_pdbx_struct_assembly.",
// Gemmi✔️✔️:       {"id", "details", "method_details",
// Gemmi✔️✔️:        "oligomeric_details", "oligomeric_count"});
// Gemmi✔️✔️:   cif::Loop& prop_loop = block.init_mmcif_loop("_pdbx_struct_assembly_prop.",
// Gemmi✔️✔️:       {"biol_id", "type", "value"});
// Gemmi✔️✔️:   cif::Loop& gen_loop = block.init_mmcif_loop("_pdbx_struct_assembly_gen.",
// Gemmi✔️✔️:       {"assembly_id", "oper_expression", "asym_id_list"});
// Gemmi✔️✔️:   cif::Loop& oper_loop = block.init_mmcif_loop("_pdbx_struct_oper_list.",
// Gemmi✔️✔️:       {"id", "type",
// Gemmi✔️✔️:        "matrix[1][1]", "matrix[1][2]", "matrix[1][3]", "vector[1]",
// Gemmi✔️✔️:        "matrix[2][1]", "matrix[2][2]", "matrix[2][3]", "vector[2]",
// Gemmi✔️✔️:        "matrix[3][1]", "matrix[3][2]", "matrix[3][3]", "vector[3]"});
// Gemmi✔️✔️:   std::vector<const Assembly::Operator*> distinct_oper;
// Gemmi✔️✔️:   for (const Assembly& as : st.assemblies) {
// Gemmi✔️✔️:     std::string how_defined = "?";
// Gemmi✔️✔️:     if (as.author_determined && as.software_determined)
// Gemmi✔️✔️:       how_defined = "author_and_software_defined_assembly";
// Gemmi✔️✔️:     else if (as.author_determined)
// Gemmi✔️✔️:       how_defined = "author_defined_assembly";
// Gemmi✔️✔️:     else if (as.software_determined)
// Gemmi✔️✔️:       how_defined = "software_defined_assembly";
// Gemmi✔️✔️:     else if (as.special_kind == Assembly::SpecialKind::CompleteIcosahedral)
// Gemmi✔️✔️:       how_defined = "'complete icosahedral assembly'";
// Gemmi✔️✔️:     else if (as.special_kind == Assembly::SpecialKind::RepresentativeHelical)
// Gemmi✔️✔️:       how_defined = "'representative helical assembly'";
// Gemmi✔️✔️:     else if (as.special_kind == Assembly::SpecialKind::CompletePoint)
// Gemmi✔️✔️:       how_defined = "'complete point assembly'";
// Gemmi✔️✔️:     std::string oligomer = to_lower(as.oligomeric_details);
// Gemmi✔️✔️:     int nmer = as.oligomeric_count != 0 ? as.oligomeric_count
// Gemmi✔️✔️:                                         : xmeric_to_number(oligomer);
// Gemmi✔️✔️:     // _pdbx_struct_assembly
// Gemmi✔️✔️:     a_loop.add_row({as.name,
// Gemmi✔️✔️:                     how_defined,
// Gemmi✔️✔️:                     string_or_qmark(as.software_name),
// Gemmi✔️✔️:                     string_or_qmark(oligomer),
// Gemmi✔️✔️:                     nmer == 0 ? "?" : std::to_string(nmer)});
// Gemmi✔️✔️:
// Gemmi✔️✔️:     // _pdbx_struct_assembly_prop
// Gemmi✔️✔️:     if (!std::isnan(as.absa))
// Gemmi✔️✔️:       prop_loop.add_row({as.name, "'ABSA (A^2)'", to_str(as.absa)});
// Gemmi✔️✔️:     if (!std::isnan(as.ssa))
// Gemmi✔️✔️:       prop_loop.add_row({as.name, "'SSA (A^2)'", to_str(as.ssa)});
// Gemmi✔️✔️:     if (!std::isnan(as.more))
// Gemmi✔️✔️:       prop_loop.add_row({as.name, "MORE", to_str(as.more)});
// Gemmi✔️✔️:
// Gemmi✔️✔️:     // _pdbx_struct_assembly_gen and _pdbx_struct_oper_list
// Gemmi✔️✔️:     for (const Assembly::Gen& gen : as.generators) {
// Gemmi✔️✔️:       std::string subchain_str;
// Gemmi✔️✔️:       for (const std::string& name : gen.subchains)
// Gemmi✔️✔️:         string_append_sep(subchain_str, ',', name);
// Gemmi✔️✔️:       if (subchain_str.empty()) // chain names to subchain names
// Gemmi✔️✔️:         for (const Chain& chain : st.models[0].chains)
// Gemmi✔️✔️:           if (in_vector(chain.name, gen.chains))
// Gemmi✔️✔️:             for (const auto& sub : chain.subchains())
// Gemmi✔️✔️:               string_append_sep(subchain_str, ',', sub.front().subchain);
// Gemmi✔️✔️:       std::string oper_str;
// Gemmi✔️✔️:       for (const Assembly::Operator& oper : gen.operators) {
// Gemmi✔️✔️:         size_t k = 0;
// Gemmi✔️✔️:         for (; k != distinct_oper.size(); ++k)
// Gemmi✔️✔️:           if (distinct_oper[k]->transform.approx(oper.transform, 1e-9))
// Gemmi✔️✔️:             break;
// Gemmi✔️✔️:         string_append_sep(oper_str, ',', std::to_string(k+1));
// Gemmi✔️✔️:         if (k != distinct_oper.size())
// Gemmi✔️✔️:           continue;
// Gemmi✔️✔️:         distinct_oper.emplace_back(&oper);
// Gemmi✔️✔️:         oper_loop.values.emplace_back(std::to_string(k+1));
// Gemmi✔️✔️:         if (!oper.type.empty()) {
// Gemmi✔️✔️:           oper_loop.values.emplace_back(cif::quote(oper.type));
// Gemmi✔️✔️:         } else if (oper.transform.is_identity()) {
// Gemmi✔️✔️:           oper_loop.values.emplace_back("'identity operation'");
// Gemmi✔️✔️:         } else if (as.author_determined || as.software_determined) {
// Gemmi✔️✔️:           oper_loop.values.emplace_back("'crystal symmetry operation'");
// Gemmi✔️✔️:         } else {
// Gemmi✔️✔️:           oper_loop.values.emplace_back(".");
// Gemmi✔️✔️:         }
// Gemmi✔️✔️:         for (int i = 0; i < 3; ++i) {
// Gemmi✔️✔️:           for (int j = 0; j < 3; ++j)
// Gemmi✔️✔️:             oper_loop.values.emplace_back(to_str(oper.transform.mat[i][j]));
// Gemmi✔️✔️:           oper_loop.values.emplace_back(to_str(oper.transform.vec.at(i)));
// Gemmi✔️✔️:         }
// Gemmi✔️✔️:       }
// Gemmi✔️✔️:       gen_loop.add_row({as.name,
// Gemmi✔️✔️:                         oper_str.empty() ? "." : oper_str,
// Gemmi✔️✔️:                         subchain_str.empty() ? "?" : subchain_str});
// Gemmi✔️✔️:     }
// Gemmi✔️✔️:   }
// Gemmi✔️✔️: }
// END GEMMI CPP FUNCTION
pub(super) fn write_assemblies(
    structure: &BioStructure,
    block: &mut CifBlock,
) -> Result<(), BioWriteError> {
    block.entries.reserve(4);
    block.loops.reserve(4);
    block.init_mmcif_loop(
        "_pdbx_struct_assembly.",
        &[
            "id",
            "details",
            "method_details",
            "oligomeric_details",
            "oligomeric_count",
        ],
    );
    block.init_mmcif_loop("_pdbx_struct_assembly_prop.", &["biol_id", "type", "value"]);
    block.init_mmcif_loop(
        "_pdbx_struct_assembly_gen.",
        &["assembly_id", "oper_expression", "asym_id_list"],
    );
    block.init_mmcif_loop(
        "_pdbx_struct_oper_list.",
        &[
            "id",
            "type",
            "matrix[1][1]",
            "matrix[1][2]",
            "matrix[1][3]",
            "vector[1]",
            "matrix[2][1]",
            "matrix[2][2]",
            "matrix[2][3]",
            "vector[2]",
            "matrix[3][1]",
            "matrix[3][2]",
            "matrix[3][3]",
            "vector[3]",
        ],
    );
    let assembly_loop = active_loop_index(block, "_pdbx_struct_assembly.id")?;
    let property_loop = active_loop_index(block, "_pdbx_struct_assembly_prop.biol_id")?;
    let generator_loop = active_loop_index(block, "_pdbx_struct_assembly_gen.assembly_id")?;
    let operator_loop = active_loop_index(block, "_pdbx_struct_oper_list.id")?;

    let mut distinct_operators = Vec::<&crate::bio::BioAssemblyOperator>::new();
    for assembly in &structure.assemblies {
        let how_defined = if assembly.author_determined && assembly.software_determined {
            "author_and_software_defined_assembly"
        } else if assembly.author_determined {
            "author_defined_assembly"
        } else if assembly.software_determined {
            "software_defined_assembly"
        } else {
            match assembly.special_kind {
                BioAssemblySpecialKind::CompleteIcosahedral => "'complete icosahedral assembly'",
                BioAssemblySpecialKind::RepresentativeHelical => {
                    "'representative helical assembly'"
                }
                BioAssemblySpecialKind::CompletePoint => "'complete point assembly'",
                BioAssemblySpecialKind::NA => "?",
            }
        };
        let oligomer = assembly.oligomeric_details.to_ascii_lowercase();
        let nmer = if assembly.oligomeric_count != 0 {
            assembly.oligomeric_count
        } else {
            xmeric_to_number(&oligomer)
        };
        let assembly_values = [
            assembly.name.clone(),
            how_defined.to_string(),
            if assembly.software_name.is_empty() {
                "?".to_string()
            } else {
                quote_cif_value(&assembly.software_name)
            },
            if oligomer.is_empty() {
                "?".to_string()
            } else {
                quote_cif_value(&oligomer)
            },
            if nmer == 0 {
                "?".to_string()
            } else {
                nmer.to_string()
            },
        ];
        for value in assembly_values {
            push_cif_value(&mut block.loops[assembly_loop], value);
        }

        for (type_, value) in [
            ("'ABSA (A^2)'", assembly.absa),
            ("'SSA (A^2)'", assembly.ssa),
            ("MORE", assembly.more),
        ] {
            if let Some(value) = value.filter(|value| !value.is_nan()) {
                push_cif_value(&mut block.loops[property_loop], assembly.name.clone());
                push_cif_value(&mut block.loops[property_loop], type_.to_string());
                push_cif_value(&mut block.loops[property_loop], gemmi_to_string_f64(value));
            }
        }

        for generator in &assembly.generators {
            let subchain_string = generator_subchains(structure, generator)?;
            let mut operator_string = String::new();
            for operator in &generator.operators {
                let distinct_index = distinct_operators
                    .iter()
                    .position(|distinct| {
                        transform_approx(&distinct.transform, &operator.transform, 1e-9)
                    })
                    .unwrap_or(distinct_operators.len());
                append_with_comma(&mut operator_string, &(distinct_index + 1).to_string());
                if distinct_index != distinct_operators.len() {
                    continue;
                }
                distinct_operators.push(operator);
                push_cif_value(
                    &mut block.loops[operator_loop],
                    (distinct_index + 1).to_string(),
                );
                let operator_type = if !operator.type_.is_empty() {
                    quote_cif_value(&operator.type_)
                } else if transform_is_exact_identity(&operator.transform) {
                    "'identity operation'".to_string()
                } else if assembly.author_determined || assembly.software_determined {
                    "'crystal symmetry operation'".to_string()
                } else {
                    ".".to_string()
                };
                push_cif_value(&mut block.loops[operator_loop], operator_type);
                for row in 0..3 {
                    for column in 0..3 {
                        push_cif_value(
                            &mut block.loops[operator_loop],
                            gemmi_to_string_f64(operator.transform.mat[row][column]),
                        );
                    }
                    push_cif_value(
                        &mut block.loops[operator_loop],
                        gemmi_to_string_f64(operator.transform.vec[row]),
                    );
                }
            }
            push_cif_value(&mut block.loops[generator_loop], assembly.name.clone());
            push_cif_value(
                &mut block.loops[generator_loop],
                if operator_string.is_empty() {
                    ".".to_string()
                } else {
                    operator_string
                },
            );
            push_cif_value(
                &mut block.loops[generator_loop],
                if subchain_string.is_empty() {
                    "?".to_string()
                } else {
                    subchain_string
                },
            );
        }
    }
    Ok(())
}

fn identity_transform() -> BioTransform {
    BioTransform {
        mat: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        vec: [0.0, 0.0, 0.0],
    }
}

fn push_transform_values(loop_: &mut super::cif::CifLoop, transform: &BioTransform) {
    for row in 0..3 {
        for column in 0..3 {
            push_cif_value(loop_, gemmi_to_string_f64(transform.mat[row][column]));
        }
        push_cif_value(loop_, gemmi_to_string_f64(transform.vec[row]));
    }
}

// BEGIN GEMMI CPP FUNCTION gemmi::write_ncs_oper
// Gemmi✔️✔️: void write_ncs_oper(const Structure& st, cif::Block& block) {
// Gemmi✔️✔️:   // _struct_ncs_oper (MTRIX)
// Gemmi✔️✔️:   if (st.ncs.empty())
// Gemmi✔️✔️:     return;
// Gemmi✔️✔️:   cif::Loop& ncs_oper = block.init_mmcif_loop("_struct_ncs_oper.",
// Gemmi✔️✔️:       {"id", "code",
// Gemmi✔️✔️:        "matrix[1][1]", "matrix[1][2]", "matrix[1][3]", "vector[1]",
// Gemmi✔️✔️:        "matrix[2][1]", "matrix[2][2]", "matrix[2][3]", "vector[2]",
// Gemmi✔️✔️:        "matrix[3][1]", "matrix[3][2]", "matrix[3][3]", "vector[3]"});
// Gemmi✔️✔️:   auto add_op = [&ncs_oper](const NcsOp& op) {
// Gemmi✔️✔️:     ncs_oper.values.emplace_back(op.id);
// Gemmi✔️✔️:     ncs_oper.values.emplace_back(op.given ? "given" : "generate");
// Gemmi✔️✔️:     for (int i = 0; i < 3; ++i) {
// Gemmi✔️✔️:       for (int j = 0; j < 3; ++j)
// Gemmi✔️✔️:         ncs_oper.values.emplace_back(to_str(op.tr.mat[i][j]));
// Gemmi✔️✔️:       ncs_oper.values.emplace_back(to_str(op.tr.vec.at(i)));
// Gemmi✔️✔️:     }
// Gemmi✔️✔️:   };
// Gemmi✔️✔️:   auto identity = st.info.find("_struct_ncs_oper.id");
// Gemmi✔️✔️:   if (identity != st.info.end() &&
// Gemmi✔️✔️:       !in_vector_f([&](const NcsOp& op) { return op.id == identity->second; }, st.ncs))
// Gemmi✔️✔️:     add_op(NcsOp{identity->second, true, {}});
// Gemmi✔️✔️:   for (const NcsOp& op : st.ncs)
// Gemmi✔️✔️:     add_op(op);
// Gemmi✔️✔️: }
// END GEMMI CPP FUNCTION
pub(super) fn write_ncs_oper(structure: &BioStructure, block: &mut CifBlock) {
    if structure.ncs_operators.is_empty() {
        return;
    }
    let loop_ = block.init_mmcif_loop(
        "_struct_ncs_oper.",
        &[
            "id",
            "code",
            "matrix[1][1]",
            "matrix[1][2]",
            "matrix[1][3]",
            "vector[1]",
            "matrix[2][1]",
            "matrix[2][2]",
            "matrix[2][3]",
            "vector[2]",
            "matrix[3][1]",
            "matrix[3][2]",
            "matrix[3][3]",
            "vector[3]",
        ],
    );
    if let Some(identity_id) = structure.ncs_oper_identity_id.as_deref()
        && !structure
            .ncs_operators
            .iter()
            .any(|operator| operator.id == identity_id)
    {
        push_cif_value(loop_, identity_id.to_string());
        push_cif_value(loop_, "given".to_string());
        push_transform_values(loop_, &identity_transform());
    }
    for operator in &structure.ncs_operators {
        push_cif_value(loop_, operator.id.clone());
        push_cif_value(
            loop_,
            if operator.given { "given" } else { "generate" }.to_string(),
        );
        push_transform_values(loop_, &operator.transform);
    }
}

#[derive(Debug, Clone, Copy)]
struct WriterCra {
    residue_index: usize,
    atom_index: Option<usize>,
}

// BEGIN GEMMI CPP FUNCTION gemmi::Model::find_cra
// Gemmi✔️✔️: for (Chain& chain : chains)
// Gemmi✔️✔️:   if (chain.name == address.chain_name) {
// Gemmi✔️✔️:     for (Residue& res : chain.residues)
// Gemmi✔️✔️:       if (address.res_id.matches_noseg(res) &&
// Gemmi✔️✔️:           (ignore_segment || address.res_id.segment == res.segment)) {
// Gemmi✔️✔️:         Atom *at = nullptr;
// Gemmi✔️✔️:         if (!address.atom_name.empty())
// Gemmi✔️✔️:           at = res.find_atom(address.atom_name, address.altloc);
// Gemmi✔️✔️:         return {&chain, &res, at};
// Gemmi✔️✔️:       }
// Gemmi✔️✔️:   }
// Gemmi✔️✔️: return {nullptr, nullptr, nullptr};
// Gemmi✔️✔️: bool altloc_matches(char request) const {
// Gemmi✔️✔️:   return request == '*' || altloc == '\0' || altloc == request;
// Gemmi✔️✔️: }
// END GEMMI CPP FUNCTION
fn find_cra(
    structure: &BioStructure,
    model_index: usize,
    address: &BioAtomAddress,
) -> Result<Option<WriterCra>, BioWriteError> {
    let model = structure
        .models
        .get(model_index)
        .ok_or(BioWriteError::Invariant("model index is out of bounds"))?;
    let chains = structure
        .chains
        .get(model.chain_span.start as usize..model.chain_span.end() as usize)
        .ok_or(BioWriteError::Invariant(
            "model chain span is out of bounds",
        ))?;
    for chain in chains {
        if chain_name_text(chain, "atom address chain identifier")? != address.chain_name {
            continue;
        }
        let residue_start = chain.residue_span.start as usize;
        let residues = structure
            .residues
            .get(residue_start..chain.residue_span.end() as usize)
            .ok_or(BioWriteError::Invariant(
                "chain residue span is out of bounds",
            ))?;
        for (residue_offset, residue) in residues.iter().enumerate() {
            if residue.source.seq_id != address.seq_id
                || residue_name_text(residue)? != address.residue_name
            {
                continue;
            }
            let residue_index = residue_start + residue_offset;
            let atom_index = if address.atom_name.is_empty() {
                None
            } else {
                let atom_start = residue.atom_span.start as usize;
                let atoms = structure
                    .atoms
                    .get(atom_start..residue.atom_span.end() as usize)
                    .ok_or(BioWriteError::Invariant(
                        "residue atom span is out of bounds",
                    ))?;
                let mut found = None;
                for (atom_offset, atom) in atoms.iter().enumerate() {
                    if atom_name_text(atom.name)? == address.atom_name
                        && (atom.altloc.is_none() || atom.altloc == address.altloc)
                    {
                        found = Some(atom_start + atom_offset);
                        break;
                    }
                }
                found
            };
            return Ok(Some(WriterCra {
                residue_index,
                atom_index,
            }));
        }
    }
    Ok(None)
}

#[derive(Debug, Clone, Copy)]
struct NearestImage {
    distance_sq: f64,
    pbc_shift: [i32; 3],
    symmetry_index: usize,
}

// BEGIN GEMMI CPP FUNCTIONS gemmi::Mat33::multiply / gemmi::Transform::apply
// Gemmi✔️✔️: Vec3 multiply(const Vec3& p) const {
// Gemmi✔️✔️:   return {a[0][0] * p.x + a[0][1] * p.y + a[0][2] * p.z,
// Gemmi✔️✔️:           a[1][0] * p.x + a[1][1] * p.y + a[1][2] * p.z,
// Gemmi✔️✔️:           a[2][0] * p.x + a[2][1] * p.y + a[2][2] * p.z};
// Gemmi✔️✔️: }
// Gemmi✔️✔️: Vec3 apply(const Vec3& x) const { return mat.multiply(x) + vec; }
// END GEMMI CPP FUNCTIONS
fn transform_apply(transform: &BioTransform, value: [f64; 3]) -> [f64; 3] {
    let multiplied = [
        transform.mat[0][0] * value[0]
            + transform.mat[0][1] * value[1]
            + transform.mat[0][2] * value[2],
        transform.mat[1][0] * value[0]
            + transform.mat[1][1] * value[1]
            + transform.mat[1][2] * value[2],
        transform.mat[2][0] * value[0]
            + transform.mat[2][1] * value[1]
            + transform.mat[2][2] * value[2],
    ];
    [
        multiplied[0] + transform.vec[0],
        multiplied[1] + transform.vec[1],
        multiplied[2] + transform.vec[2],
    ]
}

// BEGIN GEMMI CPP FUNCTIONS gemmi::Vec3_::length_sq / gemmi::Vec3_::dist_sq
// Gemmi✔️✔️: Real length_sq() const { return x * x + y * y + z * z; }
// Gemmi✔️✔️: Real dist_sq(const Vec3_& o) const { return (*this - o).length_sq(); }
// END GEMMI CPP FUNCTIONS
fn vector_distance_sq(left: [f64; 3], right: [f64; 3]) -> f64 {
    let difference = [left[0] - right[0], left[1] - right[1], left[2] - right[2]];
    difference[0] * difference[0] + difference[1] * difference[1] + difference[2] * difference[2]
}

// BEGIN GEMMI CPP FUNCTION gemmi::UnitCell::is_crystal
// Gemmi✔️✔️: bool is_crystal() const { return a != 1.0 && frac.mat[0][0] != 1.0; }
// END GEMMI CPP FUNCTION
fn crystal_is_crystal(structure: &BioStructure) -> bool {
    structure
        .crystal
        .as_ref()
        .is_some_and(|crystal| crystal.cell.a != 1.0 && crystal.frac.mat[0][0] != 1.0)
}

// BEGIN GEMMI CPP FUNCTION gemmi::UnitCell::search_pbc_images
// Gemmi✔️✔️: int neg_shift[3] = {0, 0, 0};
// Gemmi✔️✔️: if (is_crystal()) {
// Gemmi✔️✔️:   for (int j = 0; j < 3; ++j)
// Gemmi✔️✔️:     neg_shift[j] = iround(diff.at(j));
// Gemmi✔️✔️:   diff.x -= neg_shift[0];
// Gemmi✔️✔️:   diff.y -= neg_shift[1];
// Gemmi✔️✔️:   diff.z -= neg_shift[2];
// Gemmi✔️✔️: }
// Gemmi✔️✔️: Position orth_diff = orthogonalize_difference(diff);
// Gemmi✔️✔️: double dsq = orth_diff.length_sq();
// Gemmi✔️✔️: if (dsq < image.dist_sq) {
// Gemmi✔️✔️:   image.dist_sq = dsq;
// Gemmi✔️✔️:   for (int j = 0; j < 3; ++j)
// Gemmi✔️✔️:     image.pbc_shift[j] = -neg_shift[j];
// Gemmi✔️✔️:   return true;
// Gemmi✔️✔️: }
// Gemmi✔️✔️: return false;
// Gemmi✔️✔️: inline int iround(double d) { return static_cast<int>(std::round(d)); }
// END GEMMI CPP FUNCTION
fn search_pbc_images(
    structure: &BioStructure,
    mut difference: [f64; 3],
    image: &mut NearestImage,
) -> bool {
    let mut negative_shift = [0_i32; 3];
    if crystal_is_crystal(structure) {
        for axis in 0..3 {
            negative_shift[axis] = difference[axis].round() as i32;
            difference[axis] -= f64::from(negative_shift[axis]);
        }
    }
    let orthogonal_difference = structure.crystal.as_ref().map_or(difference, |crystal| {
        let mut transform = crystal.orth;
        transform.vec = [0.0; 3];
        transform_apply(&transform, difference)
    });
    let distance_sq = orthogonal_difference
        .iter()
        .map(|value| value * value)
        .sum();
    if distance_sq < image.distance_sq {
        image.distance_sq = distance_sq;
        for axis in 0..3 {
            image.pbc_shift[axis] = -negative_shift[axis];
        }
        return true;
    }
    false
}

// BEGIN GEMMI CPP FUNCTION gemmi::UnitCell::find_nearest_image
// Gemmi✔️✔️: NearestImage image;
// Gemmi✔️✔️: if (asu == Asu::Different)
// Gemmi✔️✔️:   image.dist_sq = INFINITY;
// Gemmi✔️✔️: else
// Gemmi✔️✔️:   image.dist_sq = ref.dist_sq(pos);
// Gemmi✔️✔️: if (asu == Asu::Same)
// Gemmi✔️✔️:   return image;
// Gemmi✔️✔️: Fractional fpos = fractionalize(pos);
// Gemmi✔️✔️: Fractional fref = fractionalize(ref);
// Gemmi✔️✔️: search_pbc_images(fpos - fref, image);
// Gemmi✔️✔️: if (asu == Asu::Different &&
// Gemmi✔️✔️:     image.pbc_shift[0] == 0 && image.pbc_shift[1] == 0 && image.pbc_shift[2] == 0)
// Gemmi✔️✔️:   image.dist_sq = INFINITY;
// Gemmi✔️✔️: for (int n = 0; n != static_cast<int>(images.size()); ++n)
// Gemmi✔️✔️:   if (search_pbc_images(images[n].apply(fpos) - fref, image))
// Gemmi✔️✔️:     image.sym_idx = n + 1;
// Gemmi✔️✔️: return image;
// END GEMMI CPP FUNCTION
fn find_nearest_image(
    structure: &BioStructure,
    reference: [f64; 3],
    position: [f64; 3],
    asu: BioAsu,
) -> NearestImage {
    let mut image = NearestImage {
        distance_sq: if asu == BioAsu::Different {
            f64::INFINITY
        } else {
            vector_distance_sq(reference, position)
        },
        pbc_shift: [0; 3],
        symmetry_index: 0,
    };
    if asu == BioAsu::Same {
        return image;
    }
    let fractional_position = structure
        .crystal
        .as_ref()
        .map_or(position, |crystal| transform_apply(&crystal.frac, position));
    let fractional_reference = structure.crystal.as_ref().map_or(reference, |crystal| {
        transform_apply(&crystal.frac, reference)
    });
    search_pbc_images(
        structure,
        [
            fractional_position[0] - fractional_reference[0],
            fractional_position[1] - fractional_reference[1],
            fractional_position[2] - fractional_reference[2],
        ],
        &mut image,
    );
    if asu == BioAsu::Different && image.pbc_shift == [0; 3] {
        image.distance_sq = f64::INFINITY;
    }
    if let Some(crystal) = structure.crystal.as_ref() {
        for (index, transform) in crystal.cell_images.iter().enumerate() {
            let transformed = transform_apply(transform, fractional_position);
            if search_pbc_images(
                structure,
                [
                    transformed[0] - fractional_reference[0],
                    transformed[1] - fractional_reference[1],
                    transformed[2] - fractional_reference[2],
                ],
                &mut image,
            ) {
                image.symmetry_index = index + 1;
            }
        }
    }
    image
}

// BEGIN GEMMI CPP FUNCTION gemmi::NearestImage::symmetry_code
// Gemmi✔️✔️: std::string s = std::to_string(sym_idx + 1);
// Gemmi✔️✔️: if (underscore)
// Gemmi✔️✔️:   s += '_';
// Gemmi✔️✔️: if (unsigned(5 + pbc_shift[0]) <= 9 &&
// Gemmi✔️✔️:     unsigned(5 + pbc_shift[1]) <= 9 &&
// Gemmi✔️✔️:     unsigned(5 + pbc_shift[2]) <= 9) {  // normal, quick path
// Gemmi✔️✔️:   for (int shift : pbc_shift)
// Gemmi✔️✔️:     s += char('5' + shift);
// Gemmi✔️✔️: } else {                                // problematic, non-standard path
// Gemmi✔️✔️:   for (int i = 0; i < 3; ++i) {
// Gemmi✔️✔️:     if (i != 0 && underscore)
// Gemmi✔️✔️:       s += '_';
// Gemmi✔️✔️:     s += std::to_string(5 + pbc_shift[i]);
// Gemmi✔️✔️:   }
// Gemmi✔️✔️: }
// Gemmi✔️✔️: return s;
// END GEMMI CPP FUNCTION
fn symmetry_code(image: NearestImage) -> String {
    let mut output = format!("{}_", image.symmetry_index + 1);
    if image
        .pbc_shift
        .iter()
        .all(|shift| (0..=9).contains(&(5 + shift)))
    {
        for shift in image.pbc_shift {
            output.push(char::from((5 + shift) as u8 + b'0'));
        }
    } else {
        for (index, shift) in image.pbc_shift.iter().enumerate() {
            if index != 0 {
                output.push('_');
            }
            output.push_str(&(5 + shift).to_string());
        }
    }
    output
}

// BEGIN GEMMI CPP FUNCTION gemmi::to_str_prec
// Gemmi✔️✔️: int len = d > -1e8 && d < 1e8 ? sprintf_z(buf, "%.*f", Prec, d)
// Gemmi✔️✔️:                                 : sprintf_z(buf, "%g", d);
// Gemmi✔️✔️: return std::string(buf, len > 0 ? len : 0);
// END GEMMI CPP FUNCTION
fn distance_to_string(distance: f64) -> String {
    if distance > -1e8 && distance < 1e8 {
        format_fixed(distance, 4)
    } else {
        format_general(distance, 6)
    }
}

fn connection_type_text(type_: BioConnectionType) -> &'static str {
    match type_ {
        BioConnectionType::Covale => "covale",
        BioConnectionType::Disulf => "disulf",
        BioConnectionType::Hydrog => "hydrog",
        BioConnectionType::MetalC => "metalc",
        BioConnectionType::Unknown => ".",
    }
}

fn optional_quoted(value: &str) -> String {
    if value.is_empty() {
        "?".to_string()
    } else {
        quote_cif_value(value)
    }
}

fn residue_writer_values(
    structure: &BioStructure,
    cra: WriterCra,
    address: &BioAtomAddress,
) -> Result<[String; 8], BioWriteError> {
    let residue = &structure.residues[cra.residue_index];
    let atom = cra.atom_index.map(|index| &structure.atoms[index]);
    let atom_name = match atom {
        Some(atom) => quote_cif_value(&atom_name_text(atom.name)?),
        None => "?".to_string(),
    };
    Ok([
        residue
            .source
            .subchain_id
            .map_or_else(|| ".".to_string(), |id| quote_cif_value(id.as_str())),
        residue_name_text(residue)?.to_string(),
        residue
            .source
            .label_seq_id
            .map_or_else(|| ".".to_string(), |number| number.to_string()),
        atom_name,
        atom.map_or_else(
            || "?".to_string(),
            |atom| {
                atom.altloc
                    .map_or_else(|| "?".to_string(), |alt| char::from(alt.0).to_string())
            },
        ),
        quote_cif_value(&address.chain_name),
        residue
            .source
            .seq_id
            .map_or_else(|| "?".to_string(), |seq_id| seq_id.seq_num.to_string()),
        address
            .seq_id
            .and_then(|seq_id| seq_id.ins_code)
            .map_or_else(|| "?".to_string(), |code| char::from(code).to_string()),
    ])
}

// BEGIN GEMMI CPP FUNCTION gemmi::write_struct_conn
// Gemmi✔️✔️: void write_struct_conn(const Structure& st, cif::Block& block) {
// Gemmi✔️✔️:   // example:
// Gemmi✔️✔️:   // disulf1 disulf A CYS 3  SG ? 3 ? 1_555 A CYS 18 SG ? 18 ?  1_555 ? 2.045
// Gemmi✔️✔️:   std::array<bool,(int)Connection::Type::Unknown+1> type_ids{};
// Gemmi✔️✔️:   bool use_ccp4_link_id = false;
// Gemmi✔️✔️:   for (const Connection& con : st.connections)
// Gemmi✔️✔️:     if (!con.link_id.empty())
// Gemmi✔️✔️:       use_ccp4_link_id = true;
// Gemmi✔️✔️:   cif::Loop& conn_loop = block.init_mmcif_loop("_struct_conn.",
// Gemmi✔️✔️:       {"id", "conn_type_id",
// Gemmi✔️✔️:        "ptnr1_label_asym_id", "ptnr1_label_comp_id", "ptnr1_label_seq_id",
// Gemmi✔️✔️:        "ptnr1_label_atom_id", "pdbx_ptnr1_label_alt_id", "ptnr1_auth_asym_id",
// Gemmi✔️✔️:        "ptnr1_auth_seq_id", "pdbx_ptnr1_PDB_ins_code", "ptnr1_symmetry",
// Gemmi✔️✔️:        "ptnr2_label_asym_id", "ptnr2_label_comp_id", "ptnr2_label_seq_id",
// Gemmi✔️✔️:        "ptnr2_label_atom_id", "pdbx_ptnr2_label_alt_id", "ptnr2_auth_asym_id",
// Gemmi✔️✔️:        "ptnr2_auth_seq_id", "pdbx_ptnr2_PDB_ins_code", "ptnr2_symmetry",
// Gemmi✔️✔️:        "details", "pdbx_dist_value"});
// Gemmi✔️✔️:   if (use_ccp4_link_id)
// Gemmi✔️✔️:     conn_loop.tags.push_back("_struct_conn.ccp4_link_id");
// Gemmi✔️✔️:   for (const Connection& con : st.connections) {
// Gemmi✔️✔️:     const_CRA cra1 = st.models[0].find_cra(con.partner1, true);
// Gemmi✔️✔️:     const_CRA cra2 = st.models[0].find_cra(con.partner2, true);
// Gemmi✔️✔️:     if (!cra1.residue || !cra2.residue)
// Gemmi✔️✔️:       continue;
// Gemmi✔️✔️:     const Atom* at1 = cra1.atom;
// Gemmi✔️✔️:     const Atom* at2 = cra2.atom;
// Gemmi✔️✔️:     std::string im_pdb_symbol = "?", im_dist_str = "?";
// Gemmi✔️✔️:     if (at1 && at2) {
// Gemmi✔️✔️:       NearestImage im = st.cell.find_nearest_image(at1->pos, at2->pos, con.asu);
// Gemmi✔️✔️:       im_pdb_symbol = im.symmetry_code(true);
// Gemmi✔️✔️:       im_dist_str = to_str_prec<4>(im.dist());
// Gemmi✔️✔️:     }
// Gemmi✔️✔️:     auto& v = conn_loop.values;
// Gemmi✔️✔️:     v.emplace_back(string_or_qmark(con.name));            // id
// Gemmi✔️✔️:     v.emplace_back(connection_type_to_string(con.type));  // conn_type_id
// Gemmi✔️✔️:     v.emplace_back(subchain_or_dot(*cra1.residue));       // ptnr1_label_asym_id
// Gemmi✔️✔️:     v.emplace_back(cra1.residue->name);                   // ptnr1_label_comp_id
// Gemmi✔️✔️:     v.emplace_back(cra1.residue->label_seq.str('.'));     // ptnr1_label_seq_id
// Gemmi✔️✔️:     v.emplace_back(at1 ? cif::quote(at1->name) : "?");    // ptnr1_label_atom_id
// Gemmi✔️✔️:     v.emplace_back(1, at1 ? at1->altloc_or('?') : '?');   // pdbx_ptnr1_label_alt_id
// Gemmi✔️✔️:     v.emplace_back(qchain(con.partner1.chain_name));      // ptnr1_auth_asym_id
// Gemmi✔️✔️:     v.emplace_back(cra1.residue->seqid.num.str());        // ptnr1_auth_seq_id
// Gemmi✔️✔️:     v.emplace_back(pdbx_icode(con.partner1.res_id));      // ptnr1_PDB_ins_code
// Gemmi✔️✔️:     v.emplace_back("1_555");                              // ptnr1_symmetry
// Gemmi✔️✔️:     v.emplace_back(subchain_or_dot(*cra2.residue));       // ptnr2_label_asym_id
// Gemmi✔️✔️:     v.emplace_back(cra2.residue->name);                   // ptnr2_label_comp_id
// Gemmi✔️✔️:     v.emplace_back(cra2.residue->label_seq.str('.'));     // ptnr2_label_seq_id
// Gemmi✔️✔️:     v.emplace_back(at2 ? cif::quote(at2->name) : "?");    // ptnr2_label_atom_id
// Gemmi✔️✔️:     v.emplace_back(1, at2 ? at2->altloc_or('?') : '?');   // pdbx_ptnr2_label_alt_id
// Gemmi✔️✔️:     v.emplace_back(qchain(con.partner2.chain_name));      // ptnr2_auth_asym_id
// Gemmi✔️✔️:     v.emplace_back(cra2.residue->seqid.num.str());        // ptnr2_auth_seq_id
// Gemmi✔️✔️:     v.emplace_back(pdbx_icode(con.partner2.res_id));      // ptnr2_PDB_ins_code
// Gemmi✔️✔️:     v.emplace_back(im_pdb_symbol);                        // ptnr2_symmetry
// Gemmi✔️✔️:     v.emplace_back("?");                                  // details
// Gemmi✔️✔️:     v.emplace_back(im_dist_str);                          // pdbx_dist_value
// Gemmi✔️✔️:     if (use_ccp4_link_id)
// Gemmi✔️✔️:       v.emplace_back(string_or_qmark(con.link_id));       // ccp4_link_id
// Gemmi✔️✔️:     type_ids[int(con.type)] = true;
// Gemmi✔️✔️:   }
// Gemmi✔️✔️:
// Gemmi✔️✔️:   cif::Loop& type_loop = block.init_mmcif_loop("_struct_conn_type.", {"id"});
// Gemmi✔️✔️:   for (int i = 0; i < (int)type_ids.size() - 1; ++i)
// Gemmi✔️✔️:     if (type_ids[i])
// Gemmi✔️✔️:       type_loop.add_row({connection_type_to_string((Connection::Type)i)});
// Gemmi✔️✔️: }
// END GEMMI CPP FUNCTION
pub(super) fn write_struct_conn(
    structure: &BioStructure,
    block: &mut CifBlock,
) -> Result<(), BioWriteError> {
    let use_ccp4_link_id = structure
        .connections
        .iter()
        .any(|connection| !connection.link_id.is_empty());
    let conn_loop = block.init_mmcif_loop(
        "_struct_conn.",
        &[
            "id",
            "conn_type_id",
            "ptnr1_label_asym_id",
            "ptnr1_label_comp_id",
            "ptnr1_label_seq_id",
            "ptnr1_label_atom_id",
            "pdbx_ptnr1_label_alt_id",
            "ptnr1_auth_asym_id",
            "ptnr1_auth_seq_id",
            "pdbx_ptnr1_PDB_ins_code",
            "ptnr1_symmetry",
            "ptnr2_label_asym_id",
            "ptnr2_label_comp_id",
            "ptnr2_label_seq_id",
            "ptnr2_label_atom_id",
            "pdbx_ptnr2_label_alt_id",
            "ptnr2_auth_asym_id",
            "ptnr2_auth_seq_id",
            "pdbx_ptnr2_PDB_ins_code",
            "ptnr2_symmetry",
            "details",
            "pdbx_dist_value",
        ],
    );
    if use_ccp4_link_id {
        conn_loop.tags.push("_struct_conn.ccp4_link_id".to_string());
    }
    let mut used_types = [false; 5];
    for connection in &structure.connections {
        let Some(left) = find_cra(structure, 0, &connection.partner1)? else {
            continue;
        };
        let Some(right) = find_cra(structure, 0, &connection.partner2)? else {
            continue;
        };
        let (symmetry, distance) = match (left.atom_index, right.atom_index) {
            (Some(left_atom), Some(right_atom)) => {
                let left_position = structure
                    .coordinates
                    .positions
                    .get(left_atom)
                    .copied()
                    .ok_or(BioWriteError::Invariant(
                        "connection atom coordinate is out of bounds",
                    ))?;
                let right_position = structure
                    .coordinates
                    .positions
                    .get(right_atom)
                    .copied()
                    .ok_or(BioWriteError::Invariant(
                        "connection atom coordinate is out of bounds",
                    ))?;
                let image =
                    find_nearest_image(structure, left_position, right_position, connection.asu);
                (
                    symmetry_code(image),
                    distance_to_string(image.distance_sq.sqrt()),
                )
            }
            _ => ("?".to_string(), "?".to_string()),
        };
        let left_values = residue_writer_values(structure, left, &connection.partner1)?;
        let right_values = residue_writer_values(structure, right, &connection.partner2)?;
        let mut row = Vec::with_capacity(if use_ccp4_link_id { 23 } else { 22 });
        row.push(optional_quoted(&connection.name));
        row.push(connection_type_text(connection.type_).to_string());
        row.extend(left_values[..8].iter().cloned());
        row.push("1_555".to_string());
        row.extend(right_values[..8].iter().cloned());
        row.push(symmetry);
        row.push("?".to_string());
        row.push(distance);
        if use_ccp4_link_id {
            row.push(optional_quoted(&connection.link_id));
        }
        for value in row {
            push_cif_value(conn_loop, value);
        }
        used_types[match connection.type_ {
            BioConnectionType::Covale => 0,
            BioConnectionType::Disulf => 1,
            BioConnectionType::Hydrog => 2,
            BioConnectionType::MetalC => 3,
            BioConnectionType::Unknown => 4,
        }] = true;
    }

    let type_loop = block.init_mmcif_loop("_struct_conn_type.", &["id"]);
    for (index, type_) in [
        BioConnectionType::Covale,
        BioConnectionType::Disulf,
        BioConnectionType::Hydrog,
        BioConnectionType::MetalC,
    ]
    .into_iter()
    .enumerate()
    {
        if used_types[index] {
            push_cif_value(type_loop, connection_type_text(type_).to_string());
        }
    }
    Ok(())
}

// BEGIN GEMMI CPP FUNCTION gemmi::write_cispeps
// Gemmi✔️✔️: void write_cispeps(const Structure& st, cif::Block& block) {
// Gemmi✔️✔️:   cif::Loop* prot_cis_loop = nullptr;
// Gemmi✔️✔️:   int pdbx_id = 0;
// Gemmi✔️✔️:   for (const CisPep& cispep : st.cispeps) {
// Gemmi✔️✔️:     const Model* model = &st.models[0];
// Gemmi✔️✔️:     if (st.models.size() > 1) {
// Gemmi✔️✔️:       model = st.find_model(cispep.model_num);
// Gemmi✔️✔️:       if (!model)
// Gemmi✔️✔️:         continue;
// Gemmi✔️✔️:     }
// Gemmi✔️✔️:     const_CRA cra1 = model->find_cra(cispep.partner_c, true);
// Gemmi✔️✔️:     const_CRA cra2 = model->find_cra(cispep.partner_n, true);
// Gemmi✔️✔️:     if (!cra1.residue || !cra2.residue)
// Gemmi✔️✔️:       continue;
// Gemmi✔️✔️:     if (!prot_cis_loop)
// Gemmi✔️✔️:       prot_cis_loop = &block.init_mmcif_loop("_struct_mon_prot_cis.",
// Gemmi✔️✔️:           {"pdbx_id", "pdbx_PDB_model_num",
// Gemmi✔️✔️:            "label_asym_id", "label_seq_id", "label_comp_id",
// Gemmi✔️✔️:            "auth_asym_id", "auth_seq_id", "pdbx_PDB_ins_code",
// Gemmi✔️✔️:            "pdbx_label_asym_id_2", "pdbx_label_seq_id_2", "pdbx_label_comp_id_2",
// Gemmi✔️✔️:            "pdbx_auth_asym_id_2", "pdbx_auth_seq_id_2", "pdbx_PDB_ins_code_2",
// Gemmi✔️✔️:            "label_alt_id", "pdbx_omega_angle"});
// Gemmi✔️✔️:     auto& v = prot_cis_loop->values;
// Gemmi✔️✔️:     v.emplace_back(std::to_string(++pdbx_id));            // pdbx_id
// Gemmi✔️✔️:     v.emplace_back(std::to_string(model->num));           // pdbx_PDB_model_num
// Gemmi✔️✔️:     v.emplace_back(subchain_or_dot(*cra1.residue));       // label_asym_id
// Gemmi✔️✔️:     v.emplace_back(cra1.residue->label_seq.str('.'));     // label_seq_id
// Gemmi✔️✔️:     v.emplace_back(cra1.residue->name);                   // label_comp_id
// Gemmi✔️✔️:     v.emplace_back(qchain(cispep.partner_c.chain_name));  // auth_asym_id
// Gemmi✔️✔️:     v.emplace_back(cispep.partner_c.res_id.seqid.num.str()); // auth_seq_id
// Gemmi✔️✔️:     v.emplace_back(pdbx_icode(cispep.partner_c.res_id));  // pdbx_PDB_ins_code
// Gemmi✔️✔️:     v.emplace_back(subchain_or_dot(*cra2.residue));       // pdbx_label_asym_id_2
// Gemmi✔️✔️:     v.emplace_back(cra2.residue->label_seq.str('.'));     // pdbx_label_seq_id_2
// Gemmi✔️✔️:     v.emplace_back(cra2.residue->name);                   // pdbx_label_comp_id_2
// Gemmi✔️✔️:     v.emplace_back(qchain(cispep.partner_n.chain_name));  // pdbx_auth_asym_id_2
// Gemmi✔️✔️:     v.emplace_back(cispep.partner_n.res_id.seqid.num.str()); // auth_seq_id_2
// Gemmi✔️✔️:     v.emplace_back(pdbx_icode(cispep.partner_n.res_id));  // pdbx_PDB_ins_code_2
// Gemmi✔️✔️:     v.emplace_back(1, cispep.only_altloc ? cispep.only_altloc : '.');
// Gemmi✔️✔️:     v.emplace_back(number_or_qmark(cispep.reported_angle));
// Gemmi✔️✔️:   }
// Gemmi✔️✔️: }
// END GEMMI CPP FUNCTION
pub(super) fn write_cispeps(
    structure: &BioStructure,
    block: &mut CifBlock,
) -> Result<(), BioWriteError> {
    let mut loop_index = None;
    let mut pdbx_id = 0_i32;
    for cispep in &structure.cispeps {
        let model_index = if structure.models.len() > 1 {
            let Some(index) = structure
                .models
                .iter()
                .position(|model| model.source_model_number == Some(cispep.model_num))
            else {
                continue;
            };
            index
        } else if structure.models.is_empty() {
            return Err(BioWriteError::Invariant(
                "cis-peptide writing requires at least one model",
            ));
        } else {
            0
        };
        let Some(left) = find_cra(structure, model_index, &cispep.partner_c)? else {
            continue;
        };
        let Some(right) = find_cra(structure, model_index, &cispep.partner_n)? else {
            continue;
        };
        let index = match loop_index {
            Some(index) => index,
            None => {
                block.init_mmcif_loop(
                    "_struct_mon_prot_cis.",
                    &[
                        "pdbx_id",
                        "pdbx_PDB_model_num",
                        "label_asym_id",
                        "label_seq_id",
                        "label_comp_id",
                        "auth_asym_id",
                        "auth_seq_id",
                        "pdbx_PDB_ins_code",
                        "pdbx_label_asym_id_2",
                        "pdbx_label_seq_id_2",
                        "pdbx_label_comp_id_2",
                        "pdbx_auth_asym_id_2",
                        "pdbx_auth_seq_id_2",
                        "pdbx_PDB_ins_code_2",
                        "label_alt_id",
                        "pdbx_omega_angle",
                    ],
                );
                let index = active_loop_index(block, "_struct_mon_prot_cis.pdbx_id")?;
                loop_index = Some(index);
                index
            }
        };
        pdbx_id += 1;
        let model_number = structure.models[model_index]
            .source_model_number
            .unwrap_or(0);
        let left_values = residue_writer_values(structure, left, &cispep.partner_c)?;
        let right_values = residue_writer_values(structure, right, &cispep.partner_n)?;
        let row = [
            pdbx_id.to_string(),
            model_number.to_string(),
            left_values[0].clone(),
            left_values[2].clone(),
            left_values[1].clone(),
            left_values[5].clone(),
            cispep
                .partner_c
                .seq_id
                .map_or_else(|| "?".to_string(), |seq_id| seq_id.seq_num.to_string()),
            left_values[7].clone(),
            right_values[0].clone(),
            right_values[2].clone(),
            right_values[1].clone(),
            right_values[5].clone(),
            cispep
                .partner_n
                .seq_id
                .map_or_else(|| "?".to_string(), |seq_id| seq_id.seq_num.to_string()),
            right_values[7].clone(),
            cispep.only_altloc.map_or_else(
                || ".".to_string(),
                |altloc: AltLocLabel| char::from(altloc.0).to_string(),
            ),
            number_or_qmark(cispep.reported_angle),
        ];
        for value in row {
            push_cif_value(&mut block.loops[index], value);
        }
    }
    Ok(())
}

fn string_or_dot(value: &str) -> String {
    // Gemmi✔️✔️: return s.empty() ? "." : cif::quote(s);
    if value.is_empty() {
        ".".to_string()
    } else {
        quote_cif_value(value)
    }
}

fn string_or_qmark(value: &str) -> String {
    // Gemmi✔️✔️: return s.empty() ? "?" : cif::quote(s);
    if value.is_empty() {
        "?".to_string()
    } else {
        quote_cif_value(value)
    }
}

fn entity_kind_text(kind: EntityKind) -> &'static str {
    // Gemmi✔️✔️: case EntityType::Polymer: return "polymer";
    // Gemmi✔️✔️: case EntityType::Branched: return "branched";
    // Gemmi✔️✔️: case EntityType::NonPolymer: return "non-polymer";
    // Gemmi✔️✔️: case EntityType::Water: return "water";
    // Gemmi✔️✔️: default /*EntityType::Unknown*/: return "?";
    match kind {
        EntityKind::Polymer => "polymer",
        EntityKind::Branched => "branched",
        EntityKind::NonPolymer => "non-polymer",
        EntityKind::Water => "water",
        EntityKind::Unknown => "?",
    }
}

fn polymer_kind_text(kind: PolymerKind) -> Option<&'static str> {
    // Gemmi✔️✔️: case PolymerType::PeptideL: return "polypeptide(L)";
    // Gemmi✔️✔️: case PolymerType::Dna: return "polydeoxyribonucleotide";
    // Gemmi✔️✔️: case PolymerType::Rna: return "polyribonucleotide";
    // Gemmi✔️✔️: case PolymerType::DnaRnaHybrid:
    // Gemmi✔️✔️:   return "'polydeoxyribonucleotide/polyribonucleotide hybrid'";
    // Gemmi✔️✔️: case PolymerType::SaccharideD: return "polysaccharide(D)";
    // Gemmi✔️✔️: case PolymerType::Other: return "other";
    // Gemmi✔️✔️: default /*PolymerType::Unknown*/: return "?";
    match kind {
        PolymerKind::Peptide => Some("polypeptide(L)"),
        PolymerKind::DNA => Some("polydeoxyribonucleotide"),
        PolymerKind::RNA => Some("polyribonucleotide"),
        PolymerKind::NucleicAcidHybrid => {
            Some("'polydeoxyribonucleotide/polyribonucleotide hybrid'")
        }
        PolymerKind::Saccharide => Some("polysaccharide(D)"),
        PolymerKind::PeptideLike => Some("other"),
        PolymerKind::NonPolymer | PolymerKind::Water | PolymerKind::Unknown => None,
    }
}

fn first_monomer(value: &str) -> &str {
    // Gemmi✔️✔️: return mon_list.substr(0, mon_list.find(','));
    value.split_once(',').map_or(value, |(first, _)| first)
}

fn pdbx_one_letter_code(sequence: &[String], polymer_kind: PolymerKind) -> String {
    // Gemmi✔️✔️: std::string r;
    // Gemmi✔️✔️: for (const std::string& item : seq) {
    // Gemmi✔️✔️:   std::string code = Entity::first_mon(item);
    // Gemmi✔️✔️:   const ResidueInfo ri = find_tabulated_residue(code);
    // Gemmi✔️✔️:   if (ri.is_standard() && ri.kind == kind)
    // Gemmi✔️✔️:     r += ri.one_letter_code;
    // Gemmi✔️✔️:   else
    // Gemmi✔️✔️:     cat_to(r, '(', code, ')');
    // Gemmi✔️✔️: }
    // Gemmi✔️✔️: return r;
    let mut result = String::new();
    for item in sequence {
        let code = first_monomer(item);
        let info = crate::bio::resinfo::find_tabulated_residue(code);
        let same_kind = match polymer_kind {
            PolymerKind::Peptide | PolymerKind::PeptideLike | PolymerKind::Saccharide => {
                info.is_amino_acid()
            }
            PolymerKind::DNA => info.is_dna(),
            PolymerKind::RNA | PolymerKind::NucleicAcidHybrid => info.is_rna(),
            PolymerKind::NonPolymer | PolymerKind::Water | PolymerKind::Unknown => false,
        };
        if info.is_standard() && same_kind {
            result.push(info.one_letter_code);
        } else {
            result.push('(');
            result.push_str(code);
            result.push(')');
        }
    }
    result
}

fn add_mmcif_rows(
    block: &mut CifBlock,
    category: &str,
    suffixes: &[&str],
    rows: Vec<Vec<String>>,
) -> Result<(), BioWriteError> {
    let loop_ = block.init_mmcif_loop(category, suffixes);
    for row in rows {
        loop_.add_row(row).map_err(BioWriteError::Invariant)?;
    }
    Ok(())
}

fn first_model_subchain_to_chain(
    structure: &BioStructure,
) -> Result<Vec<(PdbChainId, String)>, BioWriteError> {
    let Some(model) = structure.models.first() else {
        return Ok(Vec::new());
    };
    let mut result = Vec::new();
    for chain_index in model.chain_span.start..model.chain_span.end() {
        let chain = structure
            .chains
            .get(chain_index as usize)
            .ok_or(BioWriteError::Invariant(
                "model chain span is out of bounds",
            ))?;
        let chain_name = chain.source.auth_chain_id.as_ref().map_or_else(
            || Ok(String::new()),
            |id| pdb_chain_text(id, "author chain identifier").map(str::to_string),
        )?;
        for residue_index in chain.residue_span.start..chain.residue_span.end() {
            let residue =
                structure
                    .residues
                    .get(residue_index as usize)
                    .ok_or(BioWriteError::Invariant(
                        "chain residue span is out of bounds",
                    ))?;
            if let Some(subchain) = residue.source.subchain_id
                && !result.iter().any(|(candidate, _)| *candidate == subchain)
            {
                result.push((subchain, chain_name.clone()));
            }
        }
    }
    Ok(result)
}

fn subchain_residues<'a>(structure: &'a BioStructure, subchain: PdbChainId) -> Vec<&'a ResidueRow> {
    structure
        .residues
        .iter()
        .filter(|residue| residue.source.subchain_id == Some(subchain))
        .collect()
}

fn label_from_auth(residues: &[&ResidueRow], auth: crate::bio::PdbSeqId) -> Option<i32> {
    residues
        .iter()
        .find(|residue| residue.source.seq_id == Some(auth))
        .and_then(|residue| residue.source.label_seq_id)
}

fn auth_from_label(residues: &[&ResidueRow], label: Option<i32>) -> Option<crate::bio::PdbSeqId> {
    let label = label?;
    residues
        .iter()
        .find(|residue| residue.source.label_seq_id == Some(label))
        .and_then(|residue| residue.source.seq_id)
}

fn seq_number_or_dot(value: Option<crate::bio::PdbSeqId>) -> String {
    value.map_or_else(|| ".".to_string(), |seq_id| seq_id.seq_num.to_string())
}

fn pdbx_icode_value(value: Option<crate::bio::PdbSeqId>) -> String {
    value
        .and_then(|seq_id| seq_id.ins_code)
        .map_or_else(|| "?".to_string(), |code| char::from(code).to_string())
}

fn get_number_obs(refinement: &crate::bio::BioRefinementInfo) -> Option<i32> {
    // Gemmi✔️✔️: int nobs = ref.reflection_count;
    // Gemmi✔️✔️: if (nobs == -1 && ref.rfree_set_count >= 0 && ref.work_set_count >= 0)
    // Gemmi✔️✔️:   nobs = ref.work_set_count + ref.rfree_set_count;
    // Gemmi✔️✔️: return nobs;
    refinement
        .reflection_count
        .or_else(|| Some(refinement.work_set_count? + refinement.rfree_set_count?))
}

fn get_number_work(refinement: &crate::bio::BioRefinementInfo) -> Option<i32> {
    // Gemmi✔️✔️: int nwork = ref.work_set_count;
    // Gemmi✔️✔️: if (nwork == -1 && ref.rfree_set_count >= 0 && ref.reflection_count >= 0)
    // Gemmi✔️✔️:   nwork = ref.reflection_count - ref.rfree_set_count;
    // Gemmi✔️✔️: return nwork;
    refinement
        .work_set_count
        .or_else(|| Some(refinement.reflection_count? - refinement.rfree_set_count?))
}

fn get_bin_number_obs(bin: &crate::bio::BioRefinementBin) -> Option<i32> {
    bin.reflection_count
        .or_else(|| Some(bin.work_set_count? + bin.rfree_set_count?))
}

fn get_bin_number_work(bin: &crate::bio::BioRefinementBin) -> Option<i32> {
    bin.work_set_count
        .or_else(|| Some(bin.reflection_count? - bin.rfree_set_count?))
}

// BEGIN GEMMI CPP FUNCTION gemmi::update_mmcif_block (primary metadata categories)
// Gemmi✔️✔️: if (st.models.empty())
// Gemmi✔️✔️:   return;
// Gemmi✔️✔️: if (groups.block_name)
// Gemmi✔️✔️:   block.name = is_valid_block_name(st.name) ? st.name : "model";
// Gemmi✔️✔️: auto e_id = st.info.find("_entry.id");
// Gemmi✔️✔️: std::string id = cif::quote(e_id != st.info.end() ? e_id->second : block.name);
// Gemmi✔️✔️: if (groups.entry)
// Gemmi✔️✔️:   block.set_pair("_entry.id", id);
// Gemmi✔️✔️: else if (const std::string* val = block.find_value("_entry.id"))
// Gemmi✔️✔️:   id = *val;
pub(super) fn write_primary_mmcif_categories(
    structure: &BioStructure,
    block: &mut CifBlock,
    groups: MmcifOutputGroups,
) -> Result<String, BioWriteError> {
    if structure.models.is_empty() {
        return Ok(String::new());
    }
    if groups.block_name {
        block.name = if is_valid_block_name(&structure.name) {
            structure.name.clone()
        } else {
            "model".to_string()
        };
    }
    let mut id = quote_cif_value(
        structure
            .metadata
            .entry_id
            .as_deref()
            .unwrap_or(&block.name),
    );
    if groups.entry {
        block.set_pair("_entry.id", id.clone());
    } else if let Some(value) = block
        .items
        .iter()
        .find(|item| item.tag.eq_ignore_ascii_case("_entry.id"))
    {
        id = value.value.value.clone();
    }

    // Gemmi✔️✔️: if (groups.database_status) {
    // Gemmi✔️✔️:   auto initial_date = st.info.find("_pdbx_database_status.recvd_initial_deposition_date");
    // Gemmi✔️✔️:   if (initial_date != st.info.end() && !initial_date->second.empty()) {
    // Gemmi✔️✔️:     span.set_pair("_pdbx_database_status.entry_id", id);
    // Gemmi✔️✔️:     span.set_pair(initial_date->first, initial_date->second);
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️: }
    if groups.database_status
        && let Some(date) = structure
            .metadata
            .received_initial_deposition_date
            .as_deref()
            .filter(|date| !date.is_empty())
    {
        block.set_pair_in_category(
            "_pdbx_database_status.",
            "_pdbx_database_status.entry_id",
            id.clone(),
        );
        block.set_pair_in_category(
            "_pdbx_database_status.",
            "_pdbx_database_status.recvd_initial_deposition_date",
            date.to_string(),
        );
    }

    // Gemmi✔️✔️: if (groups.author && !st.meta.authors.empty()) {
    // Gemmi✔️✔️:   cif::Loop& loop = block.init_mmcif_loop("_audit_author.", {"pdbx_ordinal", "name"});
    // Gemmi✔️✔️:   int n = 0;
    // Gemmi✔️✔️:   for (const std::string& author : st.meta.authors)
    // Gemmi✔️✔️:     loop.add_row({std::to_string(++n), cif::quote(author)});
    // Gemmi✔️✔️: }
    if groups.author && !structure.metadata.authors.is_empty() {
        let rows = structure
            .metadata
            .authors
            .iter()
            .enumerate()
            .map(|(index, author)| vec![(index + 1).to_string(), quote_cif_value(author)])
            .collect();
        add_mmcif_rows(block, "_audit_author.", &["pdbx_ordinal", "name"], rows)?;
    }

    // Gemmi✔️✔️: if (groups.cell) {
    // Gemmi✔️✔️:   cell_span.set_pair("_cell.entry_id", id);
    // Gemmi✔️✔️:   write_cell_parameters(st.cell, cell_span);
    // Gemmi✔️✔️:   auto z_pdb = st.info.find("_cell.Z_PDB");
    // Gemmi✔️✔️:   if (z_pdb != st.info.end())
    // Gemmi✔️✔️:     cell_span.set_pair(z_pdb->first, z_pdb->second);
    // Gemmi✔️✔️: }
    if groups.cell {
        let default_crystal;
        let crystal = if let Some(crystal) = &structure.crystal {
            crystal
        } else {
            default_crystal = super::default_crystal_info();
            &default_crystal
        };
        block.set_pair_in_category("_cell.", "_cell.entry_id", id.clone());
        write_cell_parameters(crystal.cell, block);
        if let Some(z_pdb) = &crystal.z_pdb {
            block.set_pair_in_category("_cell.", "_cell.Z_PDB", z_pdb.clone());
        }
    }

    // Gemmi✔️✔️: if (groups.symmetry) {
    // Gemmi✔️✔️:   span.set_pair("_symmetry.entry_id", id);
    // Gemmi✔️✔️:   span.set_pair("_symmetry.space_group_name_H-M",
    // Gemmi✔️✔️:                  cif::quote(st.spacegroup_hm));
    // Gemmi✔️✔️:   if (const SpaceGroup* sg = st.find_spacegroup())
    // Gemmi✔️✔️:     span.set_pair("_symmetry.Int_Tables_number", std::to_string(sg->number));
    if groups.symmetry {
        let default_crystal;
        let crystal = if let Some(crystal) = &structure.crystal {
            crystal
        } else {
            default_crystal = super::default_crystal_info();
            &default_crystal
        };
        block.set_pair_in_category("_symmetry.", "_symmetry.entry_id", id.clone());
        block.set_pair_in_category(
            "_symmetry.",
            "_symmetry.space_group_name_H-M",
            quote_cif_value(crystal.spacegroup_hm.as_deref().unwrap_or("")),
        );
        if let Some(spacegroup) = super::find_crystal_spacegroup(crystal)
            && let Some(index) = crate::io::gemmi_spacegroup_table::GEMMI_SPACEGROUPS
                .iter()
                .position(|candidate| std::ptr::eq(candidate, spacegroup))
        {
            block.set_pair_in_category(
                "_symmetry.",
                "_symmetry.Int_Tables_number",
                crate::io::gemmi_spacegroup_table::GEMMI_SPACEGROUP_NUMBERS[index].to_string(),
            );
        }
    }

    // Gemmi✔️✔️: if (groups.entity) {
    // Gemmi✔️✔️:   cif::Loop& entity_loop = block.init_mmcif_loop("_entity.", {"id", "type"});
    // Gemmi✔️✔️:   for (const Entity& ent : st.entities)
    // Gemmi✔️✔️:     entity_loop.add_row({qchain(ent.name),
    // Gemmi✔️✔️:                          entity_type_to_string(ent.entity_type)});
    // Gemmi✔️✔️: }
    if groups.entity {
        let rows = structure
            .entities
            .iter()
            .map(|entity| {
                vec![
                    quote_cif_value(&entity.source.source_entity_id),
                    entity_kind_text(entity.kind).to_string(),
                ]
            })
            .collect();
        add_mmcif_rows(block, "_entity.", &["id", "type"], rows)?;
    }

    let subchain_to_chain = if groups.entity_poly || groups.struct_ref {
        first_model_subchain_to_chain(structure)?
    } else {
        Vec::new()
    };

    // Gemmi✔️✔️: if (groups.entity_poly) {
    // Gemmi✔️✔️:   cif::Loop& ent_poly_loop = block.init_mmcif_loop("_entity_poly.",
    // Gemmi✔️✔️:       {"entity_id", "type", "pdbx_strand_id", "pdbx_seq_one_letter_code"});
    // Gemmi✔️✔️:   for (const Entity& ent : st.entities)
    // Gemmi✔️✔️:     if (ent.entity_type == EntityType::Polymer) {
    // Gemmi✔️✔️:       if (ent.polymer_type == PolymerType::Unknown)
    // Gemmi✔️✔️:         continue;
    // Gemmi✔️✔️:       std::string seq1 = pdbx_one_letter_code(ent.full_sequence, kind);
    // Gemmi✔️✔️:       std::string strand_ids;
    // Gemmi✔️✔️:       for (const std::string& sub : ent.subchains) {
    // Gemmi✔️✔️:         auto strand_id = subs_to_strands.find(sub);
    // Gemmi✔️✔️:         if (strand_id != subs_to_strands.end()) {
    // Gemmi✔️✔️:           if (!strand_ids.empty()) strand_ids += ',';
    // Gemmi✔️✔️:           strand_ids += strand_id->second;
    // Gemmi✔️✔️:         }
    // Gemmi✔️✔️:       }
    // Gemmi✔️✔️:       ent_poly_loop.add_row({qchain(ent.name), polymer_type_to_string(ent.polymer_type),
    // Gemmi✔️✔️:                              string_or_qmark(strand_ids), string_or_qmark(seq1)});
    // Gemmi✔️✔️:     }
    // Gemmi✔️✔️: }
    if groups.entity_poly {
        let mut rows = Vec::new();
        for entity in &structure.entities {
            if entity.kind != EntityKind::Polymer {
                continue;
            }
            let Some(polymer_type) = polymer_kind_text(entity.polymer_kind) else {
                continue;
            };
            let strand_ids = entity
                .subchains
                .iter()
                .filter_map(|subchain| {
                    subchain_to_chain
                        .iter()
                        .find(|(candidate, _)| candidate == subchain)
                        .map(|(_, chain)| chain.as_str())
                })
                .collect::<Vec<_>>()
                .join(",");
            let sequence = pdbx_one_letter_code(&entity.sequence, entity.polymer_kind);
            rows.push(vec![
                quote_cif_value(&entity.source.source_entity_id),
                polymer_type.to_string(),
                string_or_qmark(&strand_ids),
                string_or_qmark(&sequence),
            ]);
        }
        add_mmcif_rows(
            block,
            "_entity_poly.",
            &[
                "entity_id",
                "type",
                "pdbx_strand_id",
                "pdbx_seq_one_letter_code",
            ],
            rows,
        )?;
    }

    write_struct_ref_categories(structure, block, groups, &subchain_to_chain, &id)?;
    write_chem_comp_category(structure, block, groups)?;
    write_experiment_categories(structure, block, groups, &id)?;
    write_reflection_categories(structure, block, groups, &id)?;
    write_refinement_categories(structure, block, groups, &id)?;
    write_title_keyword_categories(structure, block, groups, &id);
    Ok(id)
}

fn write_struct_ref_categories(
    structure: &BioStructure,
    block: &mut CifBlock,
    groups: MmcifOutputGroups,
    subchain_to_chain: &[(PdbChainId, String)],
    entry_id: &str,
) -> Result<(), BioWriteError> {
    if !groups.struct_ref {
        return Ok(());
    }
    // Gemmi✔️✔️: int counter = 0;
    // Gemmi✔️✔️: int counter2 = 0;
    // Gemmi✔️✔️: for (const Entity& ent : st.entities)
    // Gemmi✔️✔️:   for (const Entity::DbRef& dbref : ent.dbrefs) {
    // Gemmi✔️✔️:     ref_loop.add_row({std::to_string(++counter),
    // Gemmi✔️✔️:                       qchain(ent.name),
    // Gemmi✔️✔️:                       string_or_dot(dbref.db_name),
    // Gemmi✔️✔️:                       string_or_dot(dbref.id_code),
    // Gemmi✔️✔️:                       string_or_qmark(dbref.accession_code),
    // Gemmi✔️✔️:                       string_or_qmark(dbref.isoform)});
    let mut ref_rows = Vec::new();
    let mut seq_rows = Vec::new();
    let mut counter = 0_usize;
    let mut align_counter = 0_usize;
    for entity in &structure.entities {
        for dbref in &entity.dbrefs {
            counter += 1;
            ref_rows.push(vec![
                counter.to_string(),
                quote_cif_value(&entity.source.source_entity_id),
                string_or_dot(&dbref.db_name),
                string_or_dot(&dbref.id_code),
                string_or_qmark(&dbref.accession_code),
                string_or_qmark(&dbref.isoform),
            ]);
            // Gemmi✔️✔️: for (const std::string& subchain : ent.subchains) {
            // Gemmi✔️✔️:   auto strand_id = subs_to_strands.find(subchain);
            // Gemmi✔️✔️:   if (strand_id == subs_to_strands.end())
            // Gemmi✔️✔️:     continue;
            for subchain in &entity.subchains {
                let Some((_, strand_id)) = subchain_to_chain
                    .iter()
                    .find(|(candidate, _)| candidate == subchain)
                else {
                    continue;
                };
                let residues = subchain_residues(structure, *subchain);
                // Gemmi✔️✔️: Residue::OptionalNum label_begin = dbref.label_seq_begin;
                // Gemmi✔️✔️: Residue::OptionalNum label_end = dbref.label_seq_end;
                // Gemmi✔️✔️: if (!label_begin || !label_end) {
                // Gemmi✔️✔️:   ConstResidueSpan span = st.models[0].get_subchain(subchain);
                // Gemmi✔️✔️:   try {
                // Gemmi✔️✔️:     label_begin = span.auth_seq_id_to_label(dbref.seq_begin);
                // Gemmi✔️✔️:     label_end = span.auth_seq_id_to_label(dbref.seq_end);
                // Gemmi✔️✔️:   } catch (const std::out_of_range&) {}
                // Gemmi✔️✔️: }
                let label_begin = dbref.label_seq_begin.or_else(|| {
                    dbref
                        .seq_begin
                        .and_then(|auth| label_from_auth(&residues, auth))
                });
                let label_end = dbref.label_seq_end.or_else(|| {
                    dbref
                        .seq_end
                        .and_then(|auth| label_from_auth(&residues, auth))
                });
                // Gemmi✔️✔️: SeqId begin = dbref.seq_begin;
                // Gemmi✔️✔️: SeqId end = dbref.seq_end;
                // Gemmi✔️✔️: if (!begin.num || !end.num) {
                // Gemmi✔️✔️:   if (const Chain* chain = st.models[0].find_chain(strand_id->second))
                // Gemmi✔️✔️:     if (ConstResidueGroup polymer = chain->get_polymer()) {
                // Gemmi✔️✔️:       begin = polymer.label_seq_id_to_auth(dbref.label_seq_begin);
                // Gemmi✔️✔️:       end = polymer.label_seq_id_to_auth(dbref.label_seq_end);
                // Gemmi✔️✔️:     }
                // Gemmi✔️✔️: }
                let begin = dbref
                    .seq_begin
                    .or_else(|| auth_from_label(&residues, dbref.label_seq_begin));
                let end = dbref
                    .seq_end
                    .or_else(|| auth_from_label(&residues, dbref.label_seq_end));
                align_counter += 1;
                // Gemmi✔️✔️: seq_loop.add_row({std::to_string(++counter2),
                // Gemmi✔️✔️:                   std::to_string(counter),
                // Gemmi✔️✔️:                   strand_id->second,
                // Gemmi✔️✔️:                   id,
                // Gemmi✔️✔️:                   label_begin.str(), label_end.str(),
                // Gemmi✔️✔️:                   string_or_qmark(dbref.accession_code),
                // Gemmi✔️✔️:                   dbref.db_begin.num.str(), dbref.db_end.num.str(),
                // Gemmi✔️✔️:                   begin.num.str(), pdbx_icode(begin),
                // Gemmi✔️✔️:                   end.num.str(), pdbx_icode(end)});
                seq_rows.push(vec![
                    align_counter.to_string(),
                    counter.to_string(),
                    strand_id.clone(),
                    entry_id.to_string(),
                    label_begin.map_or_else(|| ".".to_string(), |value| value.to_string()),
                    label_end.map_or_else(|| ".".to_string(), |value| value.to_string()),
                    string_or_qmark(&dbref.accession_code),
                    seq_number_or_dot(dbref.db_begin),
                    seq_number_or_dot(dbref.db_end),
                    seq_number_or_dot(begin),
                    pdbx_icode_value(begin),
                    seq_number_or_dot(end),
                    pdbx_icode_value(end),
                ]);
            }
        }
    }
    add_mmcif_rows(
        block,
        "_struct_ref.",
        &[
            "id",
            "entity_id",
            "db_name",
            "db_code",
            "pdbx_db_accession",
            "pdbx_db_isoform",
        ],
        ref_rows,
    )?;
    add_mmcif_rows(
        block,
        "_struct_ref_seq.",
        &[
            "align_id",
            "ref_id",
            "pdbx_strand_id",
            "pdbx_PDB_id_code",
            "seq_align_beg",
            "seq_align_end",
            "pdbx_db_accession",
            "db_align_beg",
            "db_align_end",
            "pdbx_auth_seq_align_beg",
            "pdbx_seq_align_beg_ins_code",
            "pdbx_auth_seq_align_end",
            "pdbx_seq_align_end_ins_code",
        ],
        seq_rows,
    )
}

fn write_chem_comp_category(
    structure: &BioStructure,
    block: &mut CifBlock,
    groups: MmcifOutputGroups,
) -> Result<(), BioWriteError> {
    if !groups.chem_comp {
        return Ok(());
    }
    // Gemmi✔️✔️: std::set<std::string> resnames;
    // Gemmi✔️✔️: for (const Model& model : st.models)
    // Gemmi✔️✔️:   for (const Chain& chain : model.chains)
    // Gemmi✔️✔️:     for (const Residue& res : chain.residues)
    // Gemmi✔️✔️:       resnames.insert(res.name);
    // Gemmi✔️✔️: for (const Entity& ent : st.entities)
    // Gemmi✔️✔️:   for (const std::string& item : ent.full_sequence)
    // Gemmi✔️✔️:     resnames.insert(Entity::first_mon(item));
    let mut names = BTreeSet::new();
    for residue in &structure.residues {
        names.insert(residue_name_text(residue)?.to_string());
    }
    for entity in &structure.entities {
        for item in &entity.sequence {
            names.insert(first_monomer(item).to_string());
        }
    }
    let has_shortened = !structure.shortened_ccd_codes.is_empty();
    let mut rows = Vec::with_capacity(names.len());
    // Gemmi✔️✔️: for (const std::string& name : resnames) {
    // Gemmi✔️✔️:   chem_comp_loop.values.push_back(cif::quote(name));
    // Gemmi✔️✔️:   chem_comp_loop.values.push_back(".");
    // Gemmi✔️✔️:   if (!st.shortened_ccd_codes.empty()) {
    // Gemmi✔️✔️:     chem_comp_loop.values.push_back(cif::quote(name));
    // Gemmi✔️✔️:     for (const auto& old_new : st.shortened_ccd_codes)
    // Gemmi✔️✔️:       if (old_new.second == name)
    // Gemmi✔️✔️:         chem_comp_loop.values.back() = old_new.first;
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️: }
    for name in names {
        let mut row = vec![quote_cif_value(&name), ".".to_string()];
        if has_shortened {
            let three_letter = structure
                .shortened_ccd_codes
                .iter()
                .filter(|(_, shortened)| shortened == &name)
                .map(|(original, _)| original.as_str())
                .last()
                .map_or_else(|| quote_cif_value(&name), str::to_string);
            row.push(three_letter);
        }
        rows.push(row);
    }
    let suffixes = if has_shortened {
        vec!["id", "type", "three_letter_code"]
    } else {
        vec!["id", "type"]
    };
    add_mmcif_rows(block, "_chem_comp.", &suffixes, rows)
}

fn write_experiment_categories(
    structure: &BioStructure,
    block: &mut CifBlock,
    groups: MmcifOutputGroups,
    entry_id: &str,
) -> Result<(), BioWriteError> {
    if groups.exptl {
        // Gemmi✔️✔️: if (!st.meta.experiments.empty()) {
        // Gemmi✔️✔️:   cif::Loop& loop = block.init_mmcif_loop("_exptl.",
        // Gemmi✔️✔️:                                           {"entry_id", "method", "crystals_number"});
        // Gemmi✔️✔️:   for (const ExperimentInfo& exper : st.meta.experiments)
        // Gemmi✔️✔️:     loop.add_row({id, cif::quote(exper.method),
        // Gemmi✔️✔️:                   int_or_qmark(exper.number_of_crystals)});
        // Gemmi✔️✔️: } else {
        // Gemmi✔️✔️:   auto exptl_method = st.info.find("_exptl.method");
        // Gemmi✔️✔️:   if (exptl_method != st.info.end()) {
        // Gemmi✔️✔️:     for (const std::string& m : gemmi::split_str(exptl_method->second, "; "))
        // Gemmi✔️✔️:       loop.add_row({id, cif::quote(m)});
        // Gemmi✔️✔️:   }
        // Gemmi✔️✔️: }
        if !structure.metadata.experiments.is_empty() {
            let rows = structure
                .metadata
                .experiments
                .iter()
                .map(|experiment| {
                    vec![
                        entry_id.to_string(),
                        quote_cif_value(&experiment.method),
                        int_or_qmark(experiment.number_of_crystals),
                    ]
                })
                .collect();
            add_mmcif_rows(
                block,
                "_exptl.",
                &["entry_id", "method", "crystals_number"],
                rows,
            )?;
        } else if let Some(methods) = &structure.metadata.experimental_method {
            let rows = methods
                .split("; ")
                .map(|method| vec![entry_id.to_string(), quote_cif_value(method)])
                .collect();
            add_mmcif_rows(block, "_exptl.", &["entry_id", "method"], rows)?;
        }
        // Gemmi✔️✔️: if (!st.meta.crystals.empty()) {
        // Gemmi✔️✔️:   for (const CrystalInfo& cryst : st.meta.crystals)
        // Gemmi✔️✔️:     loop.add_row({cryst.id, string_or_qmark(cryst.description)});
        // Gemmi✔️✔️: }
        if !structure.metadata.experiment_crystals.is_empty() {
            let rows = structure
                .metadata
                .experiment_crystals
                .iter()
                .map(|crystal| vec![crystal.id.clone(), string_or_qmark(&crystal.description)])
                .collect();
            add_mmcif_rows(block, "_exptl_crystal.", &["id", "description"], rows)?;
        }
        // Gemmi✔️✔️: if (std::any_of(st.meta.crystals.begin(), st.meta.crystals.end(),
        // Gemmi✔️✔️:       [](const CrystalInfo& c) { return !c.ph_range.empty() || !std::isnan(c.ph); })) {
        // Gemmi✔️✔️:   for (const CrystalInfo& crystal : st.meta.crystals)
        // Gemmi✔️✔️:     grow_loop.add_row({cif::quote(crystal.id), number_or_qmark(crystal.ph),
        // Gemmi✔️✔️:                        string_or_qmark(crystal.ph_range)});
        // Gemmi✔️✔️: }
        if structure
            .metadata
            .experiment_crystals
            .iter()
            .any(|crystal| !crystal.ph_range.is_empty() || crystal.ph.is_some())
        {
            let rows = structure
                .metadata
                .experiment_crystals
                .iter()
                .map(|crystal| {
                    vec![
                        quote_cif_value(&crystal.id),
                        number_or_qmark(crystal.ph),
                        string_or_qmark(&crystal.ph_range),
                    ]
                })
                .collect();
            add_mmcif_rows(
                block,
                "_exptl_crystal_grow.",
                &["crystal_id", "pH", "pdbx_pH_range"],
                rows,
            )?;
        }
    }
    if groups.diffrn
        && structure
            .metadata
            .experiment_crystals
            .iter()
            .any(|crystal| !crystal.diffractions.is_empty())
    {
        write_diffraction_categories(structure, block)?;
    }
    Ok(())
}

fn write_diffraction_categories(
    structure: &BioStructure,
    block: &mut CifBlock,
) -> Result<(), BioWriteError> {
    // Gemmi✔️✔️: for (const CrystalInfo& cryst : st.meta.crystals)
    // Gemmi✔️✔️:   for (const DiffractionInfo& diffr : cryst.diffractions)
    // Gemmi✔️✔️:     loop.add_row({diffr.id, cryst.id, number_or_qmark(diffr.temperature)});
    let diffractions = structure
        .metadata
        .experiment_crystals
        .iter()
        .flat_map(|crystal| {
            crystal
                .diffractions
                .iter()
                .map(move |diffraction| (crystal, diffraction))
        })
        .collect::<Vec<_>>();
    add_mmcif_rows(
        block,
        "_diffrn.",
        &["id", "crystal_id", "ambient_temp"],
        diffractions
            .iter()
            .map(|(crystal, diffraction)| {
                vec![
                    diffraction.id.clone(),
                    crystal.id.clone(),
                    number_or_qmark(diffraction.temperature),
                ]
            })
            .collect(),
    )?;
    // Gemmi✔️✔️: det_loop.add_row({diffr.id,
    // Gemmi✔️✔️:                   string_or_qmark(diffr.collection_date),
    // Gemmi✔️✔️:                   string_or_qmark(diffr.detector),
    // Gemmi✔️✔️:                   string_or_qmark(diffr.detector_make),
    // Gemmi✔️✔️:                   string_or_qmark(diffr.optics)});
    add_mmcif_rows(
        block,
        "_diffrn_detector.",
        &[
            "diffrn_id",
            "pdbx_collection_date",
            "detector",
            "type",
            "details",
        ],
        diffractions
            .iter()
            .map(|(_, diffraction)| {
                vec![
                    diffraction.id.clone(),
                    string_or_qmark(&diffraction.collection_date),
                    string_or_qmark(&diffraction.detector),
                    string_or_qmark(&diffraction.detector_make),
                    string_or_qmark(&diffraction.optics),
                ]
            })
            .collect(),
    )?;
    // Gemmi✔️✔️: rad_loop.add_row({diffr.id,
    // Gemmi✔️✔️:                   string_or_qmark(diffr.scattering_type),
    // Gemmi✔️✔️:                   std::string(1, diffr.mono_or_laue ? diffr.mono_or_laue : '?'),
    // Gemmi✔️✔️:                   string_or_qmark(diffr.monochromator)});
    add_mmcif_rows(
        block,
        "_diffrn_radiation.",
        &[
            "diffrn_id",
            "pdbx_scattering_type",
            "pdbx_monochromatic_or_laue_m_l",
            "monochromator",
        ],
        diffractions
            .iter()
            .map(|(_, diffraction)| {
                vec![
                    diffraction.id.clone(),
                    string_or_qmark(&diffraction.scattering_type),
                    diffraction.mono_or_laue.unwrap_or('?').to_string(),
                    string_or_qmark(&diffraction.monochromator),
                ]
            })
            .collect(),
    )?;
    // Gemmi✔️✔️: source_loop.add_row({diffr.id,
    // Gemmi✔️✔️:                      string_or_qmark(diffr.source),
    // Gemmi✔️✔️:                      string_or_qmark(diffr.source_type),
    // Gemmi✔️✔️:                      string_or_qmark(diffr.synchrotron),
    // Gemmi✔️✔️:                      string_or_qmark(diffr.beamline),
    // Gemmi✔️✔️:                      string_or_qmark(diffr.wavelengths)});
    add_mmcif_rows(
        block,
        "_diffrn_source.",
        &[
            "diffrn_id",
            "source",
            "type",
            "pdbx_synchrotron_site",
            "pdbx_synchrotron_beamline",
            "pdbx_wavelength_list",
        ],
        diffractions
            .iter()
            .map(|(_, diffraction)| {
                vec![
                    diffraction.id.clone(),
                    string_or_qmark(&diffraction.source),
                    string_or_qmark(&diffraction.source_type),
                    string_or_qmark(&diffraction.synchrotron),
                    string_or_qmark(&diffraction.beamline),
                    string_or_qmark(&diffraction.wavelengths),
                ]
            })
            .collect(),
    )
}

fn write_reflection_categories(
    structure: &BioStructure,
    block: &mut CifBlock,
    groups: MmcifOutputGroups,
    entry_id: &str,
) -> Result<(), BioWriteError> {
    if !groups.reflns || structure.metadata.experiments.is_empty() {
        return Ok(());
    }
    // Gemmi✔️✔️: int n = 0;
    // Gemmi✔️✔️: for (const ExperimentInfo& exper : st.meta.experiments)
    // Gemmi✔️✔️:   loop.add_row({id,
    // Gemmi✔️✔️:                 std::to_string(++n),
    // Gemmi✔️✔️:                 string_or_dot(join_str(exper.diffraction_ids, ",")),
    // Gemmi✔️✔️:                 int_or_qmark(exper.unique_reflections),
    // Gemmi✔️✔️:                 number_or_qmark(exper.reflections.resolution_high),
    // Gemmi✔️✔️:                 number_or_qmark(exper.reflections.resolution_low),
    // Gemmi✔️✔️:                 number_or_qmark(exper.reflections.completeness),
    // Gemmi✔️✔️:                 number_or_qmark(exper.reflections.redundancy),
    // Gemmi✔️✔️:                 number_or_qmark(exper.reflections.r_merge),
    // Gemmi✔️✔️:                 number_or_qmark(exper.reflections.r_sym),
    // Gemmi✔️✔️:                 number_or_qmark(exper.reflections.mean_I_over_sigma)});
    let rows = structure
        .metadata
        .experiments
        .iter()
        .enumerate()
        .map(|(index, experiment)| {
            vec![
                entry_id.to_string(),
                (index + 1).to_string(),
                string_or_dot(&experiment.diffraction_ids.join(",")),
                int_or_qmark(experiment.unique_reflections),
                number_or_qmark(experiment.reflections.resolution_high),
                number_or_qmark(experiment.reflections.resolution_low),
                number_or_qmark(experiment.reflections.completeness),
                number_or_qmark(experiment.reflections.redundancy),
                number_or_qmark(experiment.reflections.r_merge),
                number_or_qmark(experiment.reflections.r_sym),
                number_or_qmark(experiment.reflections.mean_i_over_sigma),
            ]
        })
        .collect();
    add_mmcif_rows(
        block,
        "_reflns.",
        &[
            "entry_id",
            "pdbx_ordinal",
            "pdbx_diffrn_id",
            "number_obs",
            "d_resolution_high",
            "d_resolution_low",
            "percent_possible_obs",
            "pdbx_redundancy",
            "pdbx_Rmerge_I_obs",
            "pdbx_Rsym_value",
            "pdbx_netI_over_sigmaI",
        ],
        rows,
    )?;
    let mut shell_rows = Vec::new();
    let mut ordinal = 0_usize;
    // Gemmi✔️✔️: for (const ExperimentInfo& exper : st.meta.experiments) {
    // Gemmi✔️✔️:   std::string diffrn_id = string_or_dot(join_str(exper.diffraction_ids, ","));
    // Gemmi✔️✔️:   for (const ReflectionsInfo& shell : exper.shells) {
    // Gemmi✔️✔️:     shell_loop->add_row({std::to_string(++n), diffrn_id,
    // Gemmi✔️✔️:                          number_or_qmark(shell.resolution_high),
    // Gemmi✔️✔️:                          number_or_qmark(shell.resolution_low),
    // Gemmi✔️✔️:                          number_or_qmark(shell.completeness),
    // Gemmi✔️✔️:                          number_or_qmark(shell.redundancy),
    // Gemmi✔️✔️:                          number_or_qmark(shell.r_merge),
    // Gemmi✔️✔️:                          number_or_qmark(shell.r_sym),
    // Gemmi✔️✔️:                          number_or_qmark(shell.mean_I_over_sigma)});
    for experiment in &structure.metadata.experiments {
        let diffraction_id = string_or_dot(&experiment.diffraction_ids.join(","));
        for shell in &experiment.shells {
            ordinal += 1;
            shell_rows.push(vec![
                ordinal.to_string(),
                diffraction_id.clone(),
                number_or_qmark(shell.resolution_high),
                number_or_qmark(shell.resolution_low),
                number_or_qmark(shell.completeness),
                number_or_qmark(shell.redundancy),
                number_or_qmark(shell.r_merge),
                number_or_qmark(shell.r_sym),
                number_or_qmark(shell.mean_i_over_sigma),
            ]);
        }
    }
    if !shell_rows.is_empty() {
        add_mmcif_rows(
            block,
            "_reflns_shell.",
            &[
                "pdbx_ordinal",
                "pdbx_diffrn_id",
                "d_res_high",
                "d_res_low",
                "percent_possible_all",
                "pdbx_redundancy",
                "Rmerge_I_obs",
                "pdbx_Rsym_value",
                "meanI_over_sigI_obs",
            ],
            shell_rows,
        )?;
    }
    Ok(())
}

fn write_refinement_categories(
    structure: &BioStructure,
    block: &mut CifBlock,
    groups: MmcifOutputGroups,
    entry_id: &str,
) -> Result<(), BioWriteError> {
    if !groups.refine || structure.metadata.refinement.is_empty() {
        return Ok(());
    }
    // Gemmi✔️✔️: bool has_shell_fsc = false;
    // Gemmi✔️✔️: bool has_shell_ffcc = false;
    // Gemmi✔️✔️: bool has_shell_iicc = false;
    // Gemmi✔️✔️: for (const RefinementInfo& ref : st.meta.refinement)
    // Gemmi✔️✔️:   for (const BasicRefinementInfo& bin : ref.bins) {
    // Gemmi✔️✔️:     if (!std::isnan(bin.fsc_work) || !std::isnan(bin.fsc_free)) has_shell_fsc = true;
    // Gemmi✔️✔️:     if (!std::isnan(bin.cc_fo_fc_work) || !std::isnan(bin.cc_fo_fc_free)) has_shell_ffcc = true;
    // Gemmi✔️✔️:     if (!std::isnan(bin.cc_intensity_work) || !std::isnan(bin.cc_intensity_free)) has_shell_iicc = true;
    // Gemmi✔️✔️:   }
    let bins = structure
        .metadata
        .refinement
        .iter()
        .flat_map(|refinement| &refinement.bins)
        .collect::<Vec<_>>();
    let has_shell_fsc = bins
        .iter()
        .any(|bin| bin.fsc_work.is_some() || bin.fsc_free.is_some());
    let has_shell_ffcc = bins
        .iter()
        .any(|bin| bin.cc_fo_fc_work.is_some() || bin.cc_fo_fc_free.is_some());
    let has_shell_iicc = bins
        .iter()
        .any(|bin| bin.cc_intensity_work.is_some() || bin.cc_intensity_free.is_some());
    let refinements = &structure.metadata.refinement;
    let has_rfree_count = refinements
        .iter()
        .any(|value| value.rfree_set_count.is_some());
    let has_r_all = refinements.iter().any(|value| value.r_all.is_some());
    let has_r_work = refinements.iter().any(|value| value.r_work.is_some());
    let has_r_free = refinements.iter().any(|value| value.r_free.is_some());
    let has_cross_validation = refinements
        .iter()
        .any(|value| !value.cross_validation_method.is_empty());
    let has_rfree_selection = refinements
        .iter()
        .any(|value| !value.rfree_selection_method.is_empty());
    let has_mean_b = refinements.iter().any(|value| value.mean_b.is_some());
    let has_aniso = refinements.iter().any(|value| !value.aniso_b.u11.is_nan());
    let has_dpi_blow_r = refinements.iter().any(|value| value.dpi_blow_r.is_some());
    let has_dpi_blow_rfree = refinements
        .iter()
        .any(|value| value.dpi_blow_rfree.is_some());
    let has_dpi_cruickshank_r = refinements
        .iter()
        .any(|value| value.dpi_cruickshank_r.is_some());
    let has_dpi_cruickshank_rfree = refinements
        .iter()
        .any(|value| value.dpi_cruickshank_rfree.is_some());
    let has_cc_fo_fc_work = refinements
        .iter()
        .any(|value| value.cc_fo_fc_work.is_some());
    let has_cc_fo_fc_free = refinements
        .iter()
        .any(|value| value.cc_fo_fc_free.is_some());
    let has_fsc_work = refinements.iter().any(|value| value.fsc_work.is_some());
    let has_fsc_free = refinements.iter().any(|value| value.fsc_free.is_some());
    let has_cc_intensity_work = refinements
        .iter()
        .any(|value| value.cc_intensity_work.is_some());
    let has_cc_intensity_free = refinements
        .iter()
        .any(|value| value.cc_intensity_free.is_some());

    let mut refine_tags = vec![
        "entry_id",
        "pdbx_refine_id",
        "ls_d_res_high",
        "ls_d_res_low",
        "ls_percent_reflns_obs",
        "ls_number_reflns_obs",
        "ls_number_reflns_R_work",
    ];
    let optional_tags = [
        (has_rfree_count, "ls_number_reflns_R_free"),
        (has_r_all, "ls_R_factor_obs"),
        (has_r_work, "ls_R_factor_R_work"),
        (has_r_free, "ls_R_factor_R_free"),
        (has_cross_validation, "pdbx_ls_cross_valid_method"),
        (has_rfree_selection, "pdbx_R_Free_selection_details"),
        (has_mean_b, "B_iso_mean"),
    ];
    for (present, tag) in optional_tags {
        if present {
            refine_tags.push(tag);
        }
    }
    if has_aniso {
        refine_tags.extend([
            "aniso_B[1][1]",
            "aniso_B[2][2]",
            "aniso_B[3][3]",
            "aniso_B[1][2]",
            "aniso_B[1][3]",
            "aniso_B[2][3]",
        ]);
    }
    for (present, tag) in [
        (has_dpi_blow_r, "pdbx_overall_SU_R_Blow_DPI"),
        (has_dpi_blow_rfree, "pdbx_overall_SU_R_free_Blow_DPI"),
        (has_dpi_cruickshank_r, "overall_SU_R_Cruickshank_DPI"),
        (
            has_dpi_cruickshank_rfree,
            "pdbx_overall_SU_R_free_Cruickshank_DPI",
        ),
        (has_cc_fo_fc_work, "correlation_coeff_Fo_to_Fc"),
        (has_cc_fo_fc_free, "correlation_coeff_Fo_to_Fc_free"),
        (has_fsc_work, "pdbx_average_fsc_work"),
        (has_fsc_free, "pdbx_average_fsc_free"),
        (has_cc_intensity_work, "correlation_coeff_I_to_Fcsqd_work"),
        (has_cc_intensity_free, "correlation_coeff_I_to_Fcsqd_free"),
    ] {
        if present {
            refine_tags.push(tag);
        }
    }
    let has_solved_by = structure
        .metadata
        .solved_by
        .as_deref()
        .is_some_and(|value| !value.is_empty());
    let has_starting_model = structure
        .metadata
        .starting_model
        .as_deref()
        .is_some_and(|value| !value.is_empty());
    if has_solved_by {
        refine_tags.push("pdbx_method_to_determine_struct");
    }
    if has_starting_model {
        refine_tags.push("pdbx_starting_model");
    }

    let mut refine_rows = Vec::with_capacity(refinements.len());
    let mut analyze_rows = Vec::new();
    let mut restraint_rows = Vec::new();
    let mut shell_rows = Vec::new();
    for refinement in refinements {
        // Gemmi✔️✔️: loop.add_values({id, cif::quote(ref.id),
        // Gemmi✔️✔️:                  number_or_dot(ref.resolution_high), number_or_dot(ref.resolution_low),
        // Gemmi✔️✔️:                  number_or_dot(ref.completeness), int_or_dot(get_number_obs(ref)),
        // Gemmi✔️✔️:                  int_or_qmark(get_number_work(ref))});
        let mut row = vec![
            entry_id.to_string(),
            quote_cif_value(&refinement.id),
            number_or_dot(refinement.resolution_high),
            number_or_dot(refinement.resolution_low),
            number_or_dot(refinement.completeness),
            int_or_dot(get_number_obs(refinement)),
            int_or_qmark(get_number_work(refinement)),
        ];
        if has_rfree_count {
            row.push(int_or_dot(refinement.rfree_set_count));
        }
        if has_r_all {
            row.push(number_or_qmark(refinement.r_all));
        }
        if has_r_work {
            row.push(number_or_qmark(refinement.r_work));
        }
        if has_r_free {
            row.push(number_or_qmark(refinement.r_free));
        }
        if has_cross_validation {
            row.push(string_or_qmark(&refinement.cross_validation_method));
        }
        if has_rfree_selection {
            row.push(string_or_qmark(&refinement.rfree_selection_method));
        }
        if has_mean_b {
            row.push(number_or_qmark(refinement.mean_b));
        }
        if has_aniso {
            row.extend([
                number_or_qmark(Some(refinement.aniso_b.u11)),
                number_or_qmark(Some(refinement.aniso_b.u22)),
                number_or_qmark(Some(refinement.aniso_b.u33)),
                number_or_qmark(Some(refinement.aniso_b.u12)),
                number_or_qmark(Some(refinement.aniso_b.u13)),
                number_or_qmark(Some(refinement.aniso_b.u23)),
            ]);
        }
        if has_dpi_blow_r {
            row.push(number_or_qmark(refinement.dpi_blow_r));
        }
        if has_dpi_blow_rfree {
            row.push(number_or_qmark(refinement.dpi_blow_rfree));
        }
        if has_dpi_cruickshank_r {
            row.push(number_or_qmark(refinement.dpi_cruickshank_r));
        }
        if has_dpi_cruickshank_rfree {
            row.push(number_or_qmark(refinement.dpi_cruickshank_rfree));
        }
        if has_cc_fo_fc_work {
            row.push(number_or_qmark(refinement.cc_fo_fc_work));
        }
        if has_cc_fo_fc_free {
            row.push(number_or_qmark(refinement.cc_fo_fc_free));
        }
        if has_fsc_work {
            row.push(number_or_qmark(refinement.fsc_work));
        }
        if has_fsc_free {
            row.push(number_or_qmark(refinement.fsc_free));
        }
        if has_cc_intensity_work {
            row.push(number_or_qmark(refinement.cc_intensity_work));
        }
        if has_cc_intensity_free {
            row.push(number_or_qmark(refinement.cc_intensity_free));
        }
        if has_solved_by {
            row.push(string_or_qmark(
                structure.metadata.solved_by.as_deref().unwrap_or(""),
            ));
        }
        if has_starting_model {
            row.push(string_or_qmark(
                structure.metadata.starting_model.as_deref().unwrap_or(""),
            ));
        }
        refine_rows.push(row);

        // Gemmi✔️✔️: if (!std::isnan(ref.luzzati_error))
        // Gemmi✔️✔️:   analyze_loop.add_row({id, cif::quote(ref.id), number_or_qmark(ref.luzzati_error)});
        if refinement.luzzati_error.is_some() {
            analyze_rows.push(vec![
                entry_id.to_string(),
                quote_cif_value(&refinement.id),
                number_or_qmark(refinement.luzzati_error),
            ]);
        }
        // Gemmi✔️✔️: for (const RefinementInfo::Restr& restr : ref.restr_stats)
        // Gemmi✔️✔️:   restr_loop.add_row({cif::quote(ref.id), cif::quote(restr.name),
        // Gemmi✔️✔️:                       int_or_qmark(restr.count), number_or_qmark(restr.weight),
        // Gemmi✔️✔️:                       string_or_qmark(restr.function), number_or_qmark(restr.dev_ideal)});
        for restraint in &refinement.restr_stats {
            restraint_rows.push(vec![
                quote_cif_value(&refinement.id),
                quote_cif_value(&restraint.name),
                int_or_qmark(restraint.count),
                number_or_qmark(restraint.weight),
                string_or_qmark(&restraint.function),
                number_or_qmark(restraint.dev_ideal),
            ]);
        }
        for bin in &refinement.bins {
            let mut shell = vec![
                quote_cif_value(&refinement.id),
                number_or_dot(bin.resolution_high),
                number_or_qmark(bin.resolution_low),
                number_or_qmark(bin.completeness),
                int_or_qmark(get_bin_number_obs(bin)),
                int_or_qmark(get_bin_number_work(bin)),
                int_or_qmark(bin.rfree_set_count),
                number_or_qmark(bin.r_all),
                number_or_qmark(bin.r_work),
                number_or_qmark(bin.r_free),
            ];
            if has_shell_fsc {
                shell.extend([number_or_qmark(bin.fsc_work), number_or_qmark(bin.fsc_free)]);
            }
            if has_shell_ffcc {
                shell.extend([
                    number_or_qmark(bin.cc_fo_fc_work),
                    number_or_qmark(bin.cc_fo_fc_free),
                ]);
            }
            if has_shell_iicc {
                shell.extend([
                    number_or_qmark(bin.cc_intensity_work),
                    number_or_qmark(bin.cc_intensity_free),
                ]);
            }
            shell_rows.push(shell);
        }
    }
    add_mmcif_rows(block, "_refine.", &refine_tags, refine_rows)?;
    add_mmcif_rows(
        block,
        "_refine_analyze.",
        &["entry_id", "pdbx_refine_id", "Luzzati_coordinate_error_obs"],
        analyze_rows,
    )?;
    add_mmcif_rows(
        block,
        "_refine_ls_restr.",
        &[
            "pdbx_refine_id",
            "type",
            "number",
            "weight",
            "pdbx_restraint_function",
            "dev_ideal",
        ],
        restraint_rows,
    )?;
    let mut shell_tags = vec![
        "pdbx_refine_id",
        "d_res_high",
        "d_res_low",
        "percent_reflns_obs",
        "number_reflns_obs",
        "number_reflns_R_work",
        "number_reflns_R_free",
        "R_factor_obs",
        "R_factor_R_work",
        "R_factor_R_free",
    ];
    if has_shell_fsc {
        shell_tags.extend(["pdbx_fsc_work", "pdbx_fsc_free"]);
    }
    if has_shell_ffcc {
        shell_tags.extend([
            "correlation_coeff_Fo_to_Fc",
            "correlation_coeff_Fo_to_Fc_free",
        ]);
    }
    if has_shell_iicc {
        shell_tags.extend([
            "correlation_coeff_I_to_Fcsqd_work",
            "correlation_coeff_I_to_Fcsqd_free",
        ]);
    }
    add_mmcif_rows(block, "_refine_ls_shell.", &shell_tags, shell_rows)
}

fn write_title_keyword_categories(
    structure: &BioStructure,
    block: &mut CifBlock,
    groups: MmcifOutputGroups,
    entry_id: &str,
) {
    if !groups.title_keywords {
        return;
    }
    // Gemmi✔️✔️: auto title = st.info.find("_struct.title");
    // Gemmi✔️✔️: if (title != st.info.end()) {
    // Gemmi✔️✔️:   span.set_pair("_struct.entry_id", id);
    // Gemmi✔️✔️:   span.set_pair(title->first, cif::quote(title->second));
    // Gemmi✔️✔️: }
    if let Some(title) = &structure.metadata.title {
        block.set_pair_in_category("_struct.", "_struct.entry_id", entry_id.to_string());
        block.set_pair_in_category("_struct.", "_struct.title", quote_cif_value(title));
    }
    // Gemmi✔️✔️: if (pdbx_keywords != st.info.end() || keywords != st.info.end())
    // Gemmi✔️✔️:   span.set_pair("_struct_keywords.entry_id", id);
    // Gemmi✔️✔️: if (pdbx_keywords != st.info.end())
    // Gemmi✔️✔️:   span.set_pair(pdbx_keywords->first, cif::quote(pdbx_keywords->second));
    // Gemmi✔️✔️: if (keywords != st.info.end())
    // Gemmi✔️✔️:   span.set_pair(keywords->first, cif::quote(keywords->second));
    if structure.metadata.pdbx_keywords.is_some() || structure.metadata.keywords.is_some() {
        block.set_pair_in_category(
            "_struct_keywords.",
            "_struct_keywords.entry_id",
            entry_id.to_string(),
        );
    }
    if let Some(keywords) = &structure.metadata.pdbx_keywords {
        block.set_pair_in_category(
            "_struct_keywords.",
            "_struct_keywords.pdbx_keywords",
            quote_cif_value(keywords),
        );
    }
    if let Some(keywords) = &structure.metadata.keywords {
        block.set_pair_in_category(
            "_struct_keywords.",
            "_struct_keywords.text",
            quote_cif_value(keywords),
        );
    }
}

fn helix_class_number(class: crate::bio::BioHelixClass) -> i32 {
    use crate::bio::BioHelixClass;
    match class {
        BioHelixClass::UnknownHelix => 0,
        BioHelixClass::RAlpha => 1,
        BioHelixClass::ROmega => 2,
        BioHelixClass::RPi => 3,
        BioHelixClass::RGamma => 4,
        BioHelixClass::R310 => 5,
        BioHelixClass::LAlpha => 6,
        BioHelixClass::LOmega => 7,
        BioHelixClass::LGamma => 8,
        BioHelixClass::Helix27 => 9,
        BioHelixClass::HelixPolyProlineNone => 10,
    }
}

fn software_classification_text(
    classification: crate::bio::BioSoftwareClassification,
) -> &'static str {
    use crate::bio::BioSoftwareClassification;
    // BEGIN GEMMI CPP FUNCTION gemmi::software_classification_to_string
    // Gemmi✔️✔️: case SoftwareItem::DataCollection: return "data collection";
    // Gemmi✔️✔️: case SoftwareItem::DataExtraction: return "data extraction";
    // Gemmi✔️✔️: case SoftwareItem::DataProcessing: return "data processing";
    // Gemmi✔️✔️: case SoftwareItem::DataReduction:  return "data reduction";
    // Gemmi✔️✔️: case SoftwareItem::DataScaling:    return "data scaling";
    // Gemmi✔️✔️: case SoftwareItem::ModelBuilding:  return "model building";
    // Gemmi✔️✔️: case SoftwareItem::Phasing:        return "phasing";
    // Gemmi✔️✔️: case SoftwareItem::Refinement:     return "refinement";
    // Gemmi✔️✔️: case SoftwareItem::Unspecified:    return "";
    // END GEMMI CPP FUNCTION
    match classification {
        BioSoftwareClassification::DataCollection => "data collection",
        BioSoftwareClassification::DataExtraction => "data extraction",
        BioSoftwareClassification::DataProcessing => "data processing",
        BioSoftwareClassification::DataReduction => "data reduction",
        BioSoftwareClassification::DataScaling => "data scaling",
        BioSoftwareClassification::ModelBuilding => "model building",
        BioSoftwareClassification::Phasing => "phasing",
        BioSoftwareClassification::Refinement => "refinement",
        BioSoftwareClassification::Unspecified => "",
    }
}

fn secondary_residue_values(
    structure: &BioStructure,
    cra: WriterCra,
    address: &BioAtomAddress,
) -> Result<[String; 6], BioWriteError> {
    let residue = &structure.residues[cra.residue_index];
    Ok([
        quote_cif_value(&address.chain_name),
        residue
            .source
            .subchain_id
            .map_or_else(|| ".".to_string(), |id| quote_cif_value(id.as_str())),
        residue_name_text(residue)?.to_string(),
        residue
            .source
            .label_seq_id
            .map_or_else(|| "?".to_string(), |number| number.to_string()),
        residue
            .source
            .seq_id
            .map_or_else(|| "?".to_string(), |seq_id| seq_id.seq_num.to_string()),
        pdbx_icode_value(address.seq_id),
    ])
}

// BEGIN GEMMI CPP FUNCTION gemmi::update_mmcif_block (secondary categories)
pub(super) fn write_secondary_mmcif_categories(
    structure: &BioStructure,
    block: &mut CifBlock,
    groups: MmcifOutputGroups,
    entry_id: &str,
) -> Result<(), BioWriteError> {
    if structure.models.is_empty() {
        return Ok(());
    }

    // Gemmi✔️✔️: if (groups.ncs)
    // Gemmi✔️✔️:   write_ncs_oper(st, block);
    if groups.ncs {
        write_ncs_oper(structure, block);
    }

    // Gemmi✔️✔️: if (groups.struct_asym) {
    // Gemmi✔️✔️:   cif::Loop& asym_loop = block.init_mmcif_loop("_struct_asym.",
    // Gemmi✔️✔️:                                                {"id", "entity_id"});
    // Gemmi✔️✔️:   for (const Chain& chain : st.models[0].chains)
    // Gemmi✔️✔️:     for (ConstResidueSpan& sub : chain.subchains()) {
    // Gemmi✔️✔️:       const std::string& sub_id = sub.subchain_id();
    // Gemmi✔️✔️:       if (!sub_id.empty()) {
    // Gemmi✔️✔️:         const Entity* ent = find_entity_of_subchain(sub_id, st.entities);
    // Gemmi✔️✔️:         asym_loop.add_row({sub_id, (ent ? qchain(ent->name) : "?")});
    // Gemmi✔️✔️:       }
    // Gemmi✔️✔️:     }
    // Gemmi✔️✔️: }
    if groups.struct_asym {
        let model = &structure.models[0];
        let mut rows = Vec::new();
        for chain_index in model.chain_span.start..model.chain_span.end() {
            let chain =
                structure
                    .chains
                    .get(chain_index as usize)
                    .ok_or(BioWriteError::Invariant(
                        "model chain span is out of bounds",
                    ))?;
            let mut previous = None;
            for residue_index in chain.residue_span.start..chain.residue_span.end() {
                let residue = structure.residues.get(residue_index as usize).ok_or(
                    BioWriteError::Invariant("chain residue span is out of bounds"),
                )?;
                let subchain = residue.source.subchain_id;
                if subchain == previous {
                    continue;
                }
                previous = subchain;
                let Some(subchain) = subchain else {
                    continue;
                };
                let entity_id = structure
                    .entities
                    .iter()
                    .find(|entity| entity.subchains.contains(&subchain))
                    .map_or_else(
                        || "?".to_string(),
                        |entity| quote_cif_value(&entity.source.source_entity_id),
                    );
                rows.push(vec![subchain.as_str().to_string(), entity_id]);
            }
        }
        add_mmcif_rows(block, "_struct_asym.", &["id", "entity_id"], rows)?;
    }

    // Gemmi✔️✔️: bool nontrivial_origx = st.has_origx && !st.origx.is_identity();
    // Gemmi✔️✔️: if (groups.origx && nontrivial_origx) { // _database_PDB_matrix (ORIGX)
    // Gemmi✔️✔️:   cif::ItemSpan span(block.items, "_database_PDB_matrix.");
    // Gemmi✔️✔️:   span.set_pair("_database_PDB_matrix.entry_id", id);
    // Gemmi✔️✔️:   for (int i = 0; i < 3; ++i) {
    // Gemmi✔️✔️:     for (int j = 0; j < 3; ++j) {
    // Gemmi✔️✔️:       span.set_pair(tag_mat, to_str(st.origx.mat[i][j]));
    // Gemmi✔️✔️:     }
    // Gemmi✔️✔️:     span.set_pair(tag_vec, to_str(st.origx.vec.at(i)));
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️: }
    let nontrivial_origx = structure.has_origx && !transform_is_exact_identity(&structure.origx);
    if groups.origx && nontrivial_origx {
        block.set_pair_in_category(
            "_database_PDB_matrix.",
            "_database_PDB_matrix.entry_id",
            entry_id.to_string(),
        );
        for row in 0..3 {
            for column in 0..3 {
                block.set_pair_in_category(
                    "_database_PDB_matrix.",
                    &format!("_database_PDB_matrix.origx[{}][{}]", row + 1, column + 1),
                    gemmi_to_string_f64(structure.origx.mat[row][column]),
                );
            }
            block.set_pair_in_category(
                "_database_PDB_matrix.",
                &format!("_database_PDB_matrix.origx_vector[{}]", row + 1),
                gemmi_to_string_f64(structure.origx.vec[row]),
            );
        }
    }

    // Gemmi✔️✔️: if (groups.struct_conf && !st.helices.empty()) {
    // Gemmi✔️✔️:   int count = 0;
    // Gemmi✔️✔️:   for (const Helix& helix : st.helices) {
    // Gemmi✔️✔️:     const_CRA cra1 = st.models[0].find_cra(helix.start);
    // Gemmi✔️✔️:     const_CRA cra2 = st.models[0].find_cra(helix.end);
    // Gemmi✔️✔️:     if (!cra1.residue || !cra2.residue)
    // Gemmi✔️✔️:       continue;
    // Gemmi✔️✔️:     struct_conf_loop.add_row({
    // Gemmi✔️✔️:       "HELX_P",
    // Gemmi✔️✔️:       "H" + std::to_string(++count),
    // Gemmi✔️✔️:       std::to_string((int)helix.pdb_helix_class),
    // Gemmi✔️✔️:       int_or_qmark(helix.length)
    // Gemmi✔️✔️:     });
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️:   if (count != 0)
    // Gemmi✔️✔️:     block.set_pair("_struct_conf_type.id", "HELX_P");
    // Gemmi✔️✔️: }
    if groups.struct_conf && !structure.helices.is_empty() {
        let mut rows = Vec::new();
        for helix in &structure.helices {
            let Some(begin) = find_cra(structure, 0, &helix.start)? else {
                continue;
            };
            let Some(end) = find_cra(structure, 0, &helix.end)? else {
                continue;
            };
            let begin_values = secondary_residue_values(structure, begin, &helix.start)?;
            let end_values = secondary_residue_values(structure, end, &helix.end)?;
            let count = rows.len() + 1;
            rows.push(vec![
                "HELX_P".to_string(),
                format!("H{count}"),
                begin_values[0].clone(),
                begin_values[1].clone(),
                begin_values[2].clone(),
                begin_values[3].clone(),
                begin_values[4].clone(),
                begin_values[5].clone(),
                end_values[0].clone(),
                end_values[1].clone(),
                end_values[2].clone(),
                end_values[3].clone(),
                end_values[4].clone(),
                end_values[5].clone(),
                helix_class_number(helix.helix_class).to_string(),
                if helix.length == -1 {
                    "?".to_string()
                } else {
                    helix.length.to_string()
                },
            ]);
        }
        if !rows.is_empty() {
            block.set_pair("_struct_conf_type.id", "HELX_P".to_string());
        }
        add_mmcif_rows(
            block,
            "_struct_conf.",
            &[
                "conf_type_id",
                "id",
                "beg_auth_asym_id",
                "beg_label_asym_id",
                "beg_label_comp_id",
                "beg_label_seq_id",
                "beg_auth_seq_id",
                "pdbx_beg_PDB_ins_code",
                "end_auth_asym_id",
                "end_label_asym_id",
                "end_label_comp_id",
                "end_label_seq_id",
                "end_auth_seq_id",
                "pdbx_end_PDB_ins_code",
                "pdbx_PDB_helix_class",
                "pdbx_PDB_helix_length",
            ],
            rows,
        )?;
    }

    write_secondary_sheet_categories(structure, block, groups)?;
    write_secondary_tail_categories(structure, block, groups, entry_id, nontrivial_origx)
}
// END GEMMI CPP FUNCTION

fn write_secondary_sheet_categories(
    structure: &BioStructure,
    block: &mut CifBlock,
    groups: MmcifOutputGroups,
) -> Result<(), BioWriteError> {
    // Gemmi✔️✔️: if (groups.struct_sheet && !st.sheets.empty()) {
    if !groups.struct_sheet || structure.sheets.is_empty() {
        return Ok(());
    }

    // Gemmi✔️✔️: cif::Loop& sheet_loop = block.init_mmcif_loop("_struct_sheet.",
    // Gemmi✔️✔️:                                                  {"id", "number_strands"});
    // Gemmi✔️✔️: for (const Sheet& sheet : st.sheets)
    // Gemmi✔️✔️:   sheet_loop.add_row({string_or_dot(sheet.name),
    // Gemmi✔️✔️:                       std::to_string(sheet.strands.size())});
    add_mmcif_rows(
        block,
        "_struct_sheet.",
        &["id", "number_strands"],
        structure
            .sheets
            .iter()
            .map(|sheet| vec![string_or_dot(&sheet.name), sheet.strands.len().to_string()])
            .collect(),
    )?;

    // Gemmi✔️✔️: cif::Loop& order_loop = block.init_mmcif_loop("_struct_sheet_order.",
    // Gemmi✔️✔️:                 {"sheet_id", "range_id_1", "range_id_2", "sense"});
    // Gemmi✔️✔️: for (const Sheet& sheet : st.sheets)
    // Gemmi✔️✔️:   for (size_t i = 1; i < sheet.strands.size(); ++i) {
    // Gemmi✔️✔️:     const Sheet::Strand& strand = sheet.strands[i];
    // Gemmi✔️✔️:     if (strand.sense != 0)
    // Gemmi✔️✔️:       order_loop.add_row({string_or_dot(sheet.name),
    // Gemmi✔️✔️:                           std::to_string(i), std::to_string(i+1),
    // Gemmi✔️✔️:                           strand.sense > 0 ? "parallel" : "anti-parallel"});
    // Gemmi✔️✔️:   }
    let mut order_rows = Vec::new();
    for sheet in &structure.sheets {
        for (index, strand) in sheet.strands.iter().enumerate().skip(1) {
            if strand.sense != 0 {
                order_rows.push(vec![
                    string_or_dot(&sheet.name),
                    index.to_string(),
                    (index + 1).to_string(),
                    if strand.sense > 0 {
                        "parallel"
                    } else {
                        "anti-parallel"
                    }
                    .to_string(),
                ]);
            }
        }
    }
    add_mmcif_rows(
        block,
        "_struct_sheet_order.",
        &["sheet_id", "range_id_1", "range_id_2", "sense"],
        order_rows,
    )?;

    // Gemmi✔️✔️: for (const Sheet& sheet : st.sheets)
    // Gemmi✔️✔️:   for (size_t i = 0; i < sheet.strands.size(); ++i) {
    // Gemmi✔️✔️:     const Sheet::Strand& strand = sheet.strands[i];
    // Gemmi✔️✔️:     const_CRA cra1 = st.models[0].find_cra(strand.start);
    // Gemmi✔️✔️:     const_CRA cra2 = st.models[0].find_cra(strand.end);
    // Gemmi✔️✔️:     if (!cra1.residue || !cra2.residue)
    // Gemmi✔️✔️:       continue;
    // Gemmi✔️✔️:     range_loop.add_row({
    // Gemmi✔️✔️:       string_or_dot(sheet.name),
    // Gemmi✔️✔️:       std::to_string(i+1),
    // Gemmi✔️✔️:     });
    // Gemmi✔️✔️:   }
    let mut range_rows = Vec::new();
    for sheet in &structure.sheets {
        for (index, strand) in sheet.strands.iter().enumerate() {
            let Some(begin) = find_cra(structure, 0, &strand.start)? else {
                continue;
            };
            let Some(end) = find_cra(structure, 0, &strand.end)? else {
                continue;
            };
            let begin_values = secondary_residue_values(structure, begin, &strand.start)?;
            let end_values = secondary_residue_values(structure, end, &strand.end)?;
            range_rows.push(vec![
                string_or_dot(&sheet.name),
                (index + 1).to_string(),
                begin_values[0].clone(),
                begin_values[1].clone(),
                begin_values[2].clone(),
                begin_values[3].clone(),
                begin_values[4].clone(),
                begin_values[5].clone(),
                end_values[0].clone(),
                end_values[1].clone(),
                end_values[2].clone(),
                end_values[3].clone(),
                end_values[4].clone(),
                end_values[5].clone(),
            ]);
        }
    }
    add_mmcif_rows(
        block,
        "_struct_sheet_range.",
        &[
            "sheet_id",
            "id",
            "beg_auth_asym_id",
            "beg_label_asym_id",
            "beg_label_comp_id",
            "beg_label_seq_id",
            "beg_auth_seq_id",
            "pdbx_beg_PDB_ins_code",
            "end_auth_asym_id",
            "end_label_asym_id",
            "end_label_comp_id",
            "end_label_seq_id",
            "end_auth_seq_id",
            "pdbx_end_PDB_ins_code",
        ],
        range_rows,
    )?;

    // Gemmi✔️✔️: for (const Sheet& sheet : st.sheets)
    // Gemmi✔️✔️:   for (size_t i = 1; i < sheet.strands.size(); ++i) {
    // Gemmi✔️✔️:     const Sheet::Strand& strand = sheet.strands[i];
    // Gemmi✔️✔️:     if (strand.hbond_atom2.atom_name.empty())
    // Gemmi✔️✔️:       continue;
    // Gemmi✔️✔️:     // hbond_atomN is not a full atom "address": altloc is missing
    // Gemmi✔️✔️:     const_CRA cra1 = st.models[0].find_cra(strand.hbond_atom1);
    // Gemmi✔️✔️:     const_CRA cra2 = st.models[0].find_cra(strand.hbond_atom2);
    // Gemmi✔️✔️:     if (!cra1.residue || !cra2.residue)
    // Gemmi✔️✔️:       continue;
    let mut hbond_rows = Vec::new();
    for sheet in &structure.sheets {
        for (index, strand) in sheet.strands.iter().enumerate().skip(1) {
            if strand.hbond_atom2.atom_name.is_empty() {
                continue;
            }
            let Some(left) = find_cra(structure, 0, &strand.hbond_atom1)? else {
                continue;
            };
            let Some(right) = find_cra(structure, 0, &strand.hbond_atom2)? else {
                continue;
            };
            let left_values = secondary_residue_values(structure, left, &strand.hbond_atom1)?;
            let right_values = secondary_residue_values(structure, right, &strand.hbond_atom2)?;
            hbond_rows.push(vec![
                string_or_dot(&sheet.name),
                index.to_string(),
                (index + 1).to_string(),
                left_values[0].clone(),
                left_values[1].clone(),
                left_values[2].clone(),
                left_values[3].clone(),
                left_values[4].clone(),
                left_values[5].clone(),
                quote_cif_value(&strand.hbond_atom1.atom_name),
                right_values[0].clone(),
                right_values[1].clone(),
                right_values[2].clone(),
                right_values[3].clone(),
                right_values[4].clone(),
                right_values[5].clone(),
                quote_cif_value(&strand.hbond_atom2.atom_name),
            ]);
        }
    }
    add_mmcif_rows(
        block,
        "_pdbx_struct_sheet_hbond.",
        &[
            "sheet_id",
            "range_id_1",
            "range_id_2",
            "range_1_auth_asym_id",
            "range_1_label_asym_id",
            "range_1_label_comp_id",
            "range_1_label_seq_id",
            "range_1_auth_seq_id",
            "range_1_PDB_ins_code",
            "range_1_label_atom_id",
            "range_2_auth_asym_id",
            "range_2_label_asym_id",
            "range_2_label_comp_id",
            "range_2_label_seq_id",
            "range_2_auth_seq_id",
            "range_2_PDB_ins_code",
            "range_2_label_atom_id",
        ],
        hbond_rows,
    )?;
    Ok(())
}

fn write_secondary_tail_categories(
    structure: &BioStructure,
    block: &mut CifBlock,
    groups: MmcifOutputGroups,
    entry_id: &str,
    nontrivial_origx: bool,
) -> Result<(), BioWriteError> {
    // Gemmi✔️✔️: // _pdbx_struct_assembly* and _struct_biol are REMARK 300/350 in PDB
    // Gemmi✔️✔️: if (groups.struct_biol && !st.meta.remark_300_detail.empty()) {
    // Gemmi✔️✔️:   cif::ItemSpan span(block.items, "_struct_biol.");
    // Gemmi✔️✔️:   span.set_pair("_struct_biol.id", "1");
    // Gemmi✔️✔️:   span.set_pair("_struct_biol.details", cif::quote(st.meta.remark_300_detail));
    // Gemmi✔️✔️: }
    if groups.struct_biol
        && let Some(details) = structure
            .metadata
            .remark_300_detail
            .as_deref()
            .filter(|details| !details.is_empty())
    {
        block.set_pair_in_category("_struct_biol.", "_struct_biol.id", "1".to_string());
        block.set_pair_in_category(
            "_struct_biol.",
            "_struct_biol.details",
            quote_cif_value(details),
        );
    }

    // Gemmi✔️✔️: if (groups.assembly && !st.assemblies.empty())
    // Gemmi✔️✔️:   write_assemblies(st, block);
    if groups.assembly && !structure.assemblies.is_empty() {
        write_assemblies(structure, block)?;
    }
    // Gemmi✔️✔️: if (groups.conn)
    // Gemmi✔️✔️:   write_struct_conn(st, block);
    if groups.conn {
        write_struct_conn(structure, block)?;
    }
    // Gemmi✔️✔️: if (groups.cis)  // _struct_mon_prot_cis
    // Gemmi✔️✔️:   write_cispeps(st, block);
    if groups.cis {
        write_cispeps(structure, block)?;
    }

    // Gemmi✔️✔️: // _pdbx_struct_mod_residue (MODRES)
    // Gemmi✔️✔️: if (groups.modres && !st.mod_residues.empty()) {
    // Gemmi✔️✔️:   bool use_ccp4_mod_id = false;
    // Gemmi✔️✔️:   for (const ModRes& modres : st.mod_residues)
    // Gemmi✔️✔️:     if (!modres.mod_id.empty())
    // Gemmi✔️✔️:       use_ccp4_mod_id = true;
    if groups.modres && !structure.mod_residues.is_empty() {
        let use_ccp4_mod_id = structure
            .mod_residues
            .iter()
            .any(|modres| !modres.mod_id.is_empty());
        let mut tags = vec![
            "id",
            "auth_asym_id",
            "auth_seq_id",
            "PDB_ins_code",
            "auth_comp_id",
            "label_comp_id",
            "parent_comp_id",
            "details",
        ];
        if use_ccp4_mod_id {
            tags.push("ccp4_mod_id");
        }
        let mut rows = Vec::with_capacity(structure.mod_residues.len());
        for (index, modres) in structure.mod_residues.iter().enumerate() {
            // Gemmi✔️✔️: loop.add_values({std::to_string(++counter),
            // Gemmi✔️✔️:                  qchain(modres.chain_name),
            // Gemmi✔️✔️:                  modres.res_id.seqid.num.str(),
            // Gemmi✔️✔️:                  pdbx_icode(modres.res_id),
            // Gemmi✔️✔️:                  string_or_dot(modres.res_id.name),
            // Gemmi✔️✔️:                  string_or_qmark(modres.res_id.name),
            // Gemmi✔️✔️:                  string_or_qmark(modres.parent_comp_id),
            // Gemmi✔️✔️:                  string_or_qmark(modres.details)});
            let mut row = vec![
                (index + 1).to_string(),
                quote_cif_value(&modres.chain_name),
                modres.res_id.seq_num.to_string(),
                pdbx_icode_value(Some(modres.res_id)),
                string_or_dot(&modres.residue_name),
                string_or_qmark(&modres.residue_name),
                string_or_qmark(&modres.parent_comp_id),
                string_or_qmark(&modres.details),
            ];
            // Gemmi✔️✔️: if (use_ccp4_mod_id)
            // Gemmi✔️✔️:   loop.values.push_back(string_or_qmark(modres.mod_id));
            if use_ccp4_mod_id {
                row.push(string_or_qmark(&modres.mod_id));
            }
            rows.push(row);
        }
        add_mmcif_rows(block, "_pdbx_struct_mod_residue.", &tags, rows)?;
    }

    // Gemmi✔️✔️: // _atom_sites (SCALE)
    // Gemmi✔️✔️: if (groups.scale && (nontrivial_origx || st.cell.explicit_matrices)) {
    // Gemmi✔️✔️:   cif::ItemSpan span(block.items, "_atom_sites.");
    // Gemmi✔️✔️:   span.set_pair("_atom_sites.entry_id", id);
    if groups.scale
        && (nontrivial_origx
            || structure
                .crystal
                .as_ref()
                .is_some_and(|crystal| crystal.explicit_matrices))
    {
        let default_crystal;
        let crystal = if let Some(crystal) = &structure.crystal {
            crystal
        } else {
            default_crystal = super::default_crystal_info();
            &default_crystal
        };
        block.set_pair_in_category("_atom_sites.", "_atom_sites.entry_id", entry_id.to_string());
        // Gemmi✔️✔️: for (int i = 0; i < 3; ++i) {
        // Gemmi✔️✔️:   span.set_pair(matrix_idx + "[1]", to_str(frac.mat[i][0]));
        // Gemmi✔️✔️:   span.set_pair(matrix_idx + "[2]", to_str(frac.mat[i][1]));
        // Gemmi✔️✔️:   span.set_pair(matrix_idx + "[3]", to_str(frac.mat[i][2]));
        // Gemmi✔️✔️:   span.set_pair(cat(prefix, "vector", idx), to_str(frac.vec.at(i)));
        // Gemmi✔️✔️: }
        for row in 0..3 {
            for column in 0..3 {
                block.set_pair_in_category(
                    "_atom_sites.",
                    &format!(
                        "_atom_sites.fract_transf_matrix[{}][{}]",
                        row + 1,
                        column + 1
                    ),
                    gemmi_to_string_f64(crystal.frac.mat[row][column]),
                );
            }
            block.set_pair_in_category(
                "_atom_sites.",
                &format!("_atom_sites.fract_transf_vector[{}]", row + 1),
                gemmi_to_string_f64(crystal.frac.vec[row]),
            );
        }
    }

    // Gemmi✔️✔️: // _atom_type
    // Gemmi✔️✔️: if (groups.atom_type) {
    // Gemmi✔️✔️:   std::array<bool, (int)El::END> types{};
    // Gemmi✔️✔️:   for (const Model& model : st.models)
    // Gemmi✔️✔️:     for (const Chain& chain : model.chains)
    // Gemmi✔️✔️:       for (const Residue& res : chain.residues)
    // Gemmi✔️✔️:         for (const Atom& atom : res.atoms)
    // Gemmi✔️✔️:           types[atom.element.ordinal()] = true;
    // Gemmi✔️✔️:   for (int i = 0; i < (int)El::END; ++i)
    // Gemmi✔️✔️:     if (types[i])
    // Gemmi✔️✔️:       atom_type_loop.add_row({Element((El)i).uname()});
    // Gemmi✔️✔️: }
    if groups.atom_type {
        let mut types = [false; 119];
        for atom in &structure.atoms {
            types[usize::from(atom.element.atomic_number())] = true;
        }
        let rows = types
            .iter()
            .enumerate()
            .filter(|(_, present)| **present)
            .map(|(atomic_number, _)| {
                let symbol = if atomic_number == 0 {
                    "X".to_string()
                } else {
                    crate::Element::from_atomic_number(atomic_number as u8)
                        .expect("BioStructure atoms always contain valid elements")
                        .symbol()
                        .to_ascii_uppercase()
                };
                vec![symbol]
            })
            .collect();
        add_mmcif_rows(block, "_atom_type.", &["symbol"], rows)?;
    }

    // Gemmi✔️✔️: if (groups.entity_poly_seq) {
    // Gemmi✔️✔️:   for (const Entity& ent : st.entities)
    // Gemmi✔️✔️:     if (ent.entity_type == EntityType::Polymer) {
    // Gemmi✔️✔️:       // SEQRES from PDB doesn't record microheterogeneity.
    // Gemmi✔️✔️:       std::string hetero_no = ent.reflects_microhetero ? "n" : "?";
    let mut poly_rows = Vec::new();
    if groups.entity_poly_seq {
        for entity in &structure.entities {
            if entity.kind != EntityKind::Polymer {
                continue;
            }
            let hetero_no = if entity.reflects_microhetero {
                "n"
            } else {
                "?"
            };
            for (index, monomers) in entity.sequence.iter().enumerate() {
                let hetero = if monomers.contains(',') {
                    "y"
                } else {
                    hetero_no
                };
                for monomer in monomers.split(',') {
                    // Gemmi✔️✔️: poly_loop.add_row({qchain(ent.name), num,
                    // Gemmi✔️✔️:                    mon_ids.substr(start),
                    // Gemmi✔️✔️:                    start == 0 ? hetero_no : "y"});
                    poly_rows.push(vec![
                        quote_cif_value(&entity.source.source_entity_id),
                        (index + 1).to_string(),
                        monomer.to_string(),
                        hetero.to_string(),
                    ]);
                }
            }
        }
        add_mmcif_rows(
            block,
            "_entity_poly_seq.",
            &["entity_id", "num", "mon_id", "hetero"],
            poly_rows,
        )?;
    }

    // Gemmi✔️✔️: if (groups.atoms)
    // Gemmi✔️✔️:   add_cif_atoms(st, block, groups.group_pdb, groups.auth_all);
    if groups.atoms {
        add_cif_atoms(structure, block, groups.group_pdb, groups.auth_all)?;
    }

    write_tls_categories(structure, block, groups)?;
    write_software_category(structure, block, groups)
}

fn write_tls_categories(
    structure: &BioStructure,
    block: &mut CifBlock,
    groups: MmcifOutputGroups,
) -> Result<(), BioWriteError> {
    // Gemmi✔️✔️: if (groups.tls && st.meta.get_tls_groups() != nullptr) {
    if !groups.tls
        || !structure
            .metadata
            .refinement
            .iter()
            .any(|refinement| !refinement.tls_groups.is_empty())
    {
        return Ok(());
    }

    // Gemmi✔️✔️: // pdbx_refine_id doesn't make sense here, but it's required
    // Gemmi✔️✔️: // by the mmCIF spec. In joint refinement, TLS constraints can't be
    // Gemmi✔️✔️: // specific to a dataset, because they constrain the shared model.
    // Gemmi✔️✔️: for (const RefinementInfo& ref : st.meta.refinement)
    // Gemmi✔️✔️:   for (const TlsGroup& tls : ref.tls_groups) {
    let mut tls_rows = Vec::new();
    for refinement in &structure.metadata.refinement {
        for tls in &refinement.tls_groups {
            // Gemmi✔️✔️: loop.add_row({string_or_dot(tls.id), cif::quote(ref.id),
            // Gemmi✔️✔️:               q(tls.origin.x), q(tls.origin.y), q(tls.origin.z),
            // Gemmi✔️✔️:               q(T.u11), q(T.u22), q(T.u33), q(T.u12), q(T.u13), q(T.u23),
            // Gemmi✔️✔️:               q(L.u11), q(L.u22), q(L.u33), q(L.u12), q(L.u13), q(L.u23),
            // Gemmi✔️✔️:               q(S[0][0]), q(S[0][1]), q(S[0][2]),
            // Gemmi✔️✔️:               q(S[1][0]), q(S[1][1]), q(S[1][2]),
            // Gemmi✔️✔️:               q(S[2][0]), q(S[2][1]), q(S[2][2])});
            tls_rows.push(vec![
                string_or_dot(&tls.id),
                quote_cif_value(&refinement.id),
                number_or_qmark(Some(tls.origin[0])),
                number_or_qmark(Some(tls.origin[1])),
                number_or_qmark(Some(tls.origin[2])),
                number_or_qmark(Some(tls.t[0][0])),
                number_or_qmark(Some(tls.t[1][1])),
                number_or_qmark(Some(tls.t[2][2])),
                number_or_qmark(Some(tls.t[0][1])),
                number_or_qmark(Some(tls.t[0][2])),
                number_or_qmark(Some(tls.t[1][2])),
                number_or_qmark(Some(tls.l[0][0])),
                number_or_qmark(Some(tls.l[1][1])),
                number_or_qmark(Some(tls.l[2][2])),
                number_or_qmark(Some(tls.l[0][1])),
                number_or_qmark(Some(tls.l[0][2])),
                number_or_qmark(Some(tls.l[1][2])),
                number_or_qmark(Some(tls.s[0][0])),
                number_or_qmark(Some(tls.s[0][1])),
                number_or_qmark(Some(tls.s[0][2])),
                number_or_qmark(Some(tls.s[1][0])),
                number_or_qmark(Some(tls.s[1][1])),
                number_or_qmark(Some(tls.s[1][2])),
                number_or_qmark(Some(tls.s[2][0])),
                number_or_qmark(Some(tls.s[2][1])),
                number_or_qmark(Some(tls.s[2][2])),
            ]);
        }
    }
    add_mmcif_rows(
        block,
        "_pdbx_refine_tls.",
        &[
            "id",
            "pdbx_refine_id",
            "origin_x",
            "origin_y",
            "origin_z",
            "T[1][1]",
            "T[2][2]",
            "T[3][3]",
            "T[1][2]",
            "T[1][3]",
            "T[2][3]",
            "L[1][1]",
            "L[2][2]",
            "L[3][3]",
            "L[1][2]",
            "L[1][3]",
            "L[2][3]",
            "S[1][1]",
            "S[1][2]",
            "S[1][3]",
            "S[2][1]",
            "S[2][2]",
            "S[2][3]",
            "S[3][1]",
            "S[3][2]",
            "S[3][3]",
        ],
        tls_rows,
    )?;

    // Gemmi✔️✔️: int counter = 1;
    // Gemmi✔️✔️: for (const RefinementInfo& ref : st.meta.refinement)
    // Gemmi✔️✔️:   for (const TlsGroup& tls : ref.tls_groups)
    // Gemmi✔️✔️:     for (const TlsGroup::Selection& sel : tls.selections)
    let mut selection_rows = Vec::new();
    for refinement in &structure.metadata.refinement {
        for tls in &refinement.tls_groups {
            for selection in &tls.selections {
                let begin = selection.res_begin;
                let end = selection.res_end;
                // Gemmi✔️✔️: group_loop.add_row({std::to_string(counter++),
                // Gemmi✔️✔️:                     string_or_dot(tls.id),
                // Gemmi✔️✔️:                     cif::quote(ref.id),
                // Gemmi✔️✔️:                     string_or_qmark(sel.chain),
                // Gemmi✔️✔️:                     sel.res_begin.num.str(),
                // Gemmi✔️✔️:                     pdbx_icode(sel.res_begin),
                // Gemmi✔️✔️:                     string_or_qmark(sel.chain),
                // Gemmi✔️✔️:                     sel.res_end.num.str(),
                // Gemmi✔️✔️:                     pdbx_icode(sel.res_end),
                // Gemmi✔️✔️:                     string_or_qmark(sel.details)});
                selection_rows.push(vec![
                    (selection_rows.len() + 1).to_string(),
                    string_or_dot(&tls.id),
                    quote_cif_value(&refinement.id),
                    string_or_qmark(&selection.chain),
                    begin.map_or_else(|| "?".to_string(), |value| value.seq_num.to_string()),
                    pdbx_icode_value(begin),
                    string_or_qmark(&selection.chain),
                    end.map_or_else(|| "?".to_string(), |value| value.seq_num.to_string()),
                    pdbx_icode_value(end),
                    string_or_qmark(&selection.details),
                ]);
            }
        }
    }
    add_mmcif_rows(
        block,
        "_pdbx_refine_tls_group.",
        &[
            "id",
            "refine_tls_id",
            "pdbx_refine_id",
            "beg_auth_asym_id",
            "beg_auth_seq_id",
            "beg_PDB_ins_code",
            "end_auth_asym_id",
            "end_auth_seq_id",
            "end_PDB_ins_code",
            "selection_details",
        ],
        selection_rows,
    )
}

fn write_software_category(
    structure: &BioStructure,
    block: &mut CifBlock,
    groups: MmcifOutputGroups,
) -> Result<(), BioWriteError> {
    // Gemmi✔️✔️: if (groups.software && !st.meta.software.empty()) {
    if !groups.software || structure.metadata.software.is_empty() {
        return Ok(());
    }
    // Gemmi✔️✔️: bool write_all_fields = false;
    // Gemmi✔️✔️: for (const SoftwareItem& item : st.meta.software)
    // Gemmi✔️✔️:   if (!item.date.empty() || !item.description.empty() ||
    // Gemmi✔️✔️:       !item.contact_author.empty() || !item.contact_author_email.empty())
    // Gemmi✔️✔️:     write_all_fields = true;
    let write_all_fields = structure.metadata.software.iter().any(|item| {
        !item.date.is_empty()
            || !item.description.is_empty()
            || !item.contact_author.is_empty()
            || !item.contact_author_email.is_empty()
    });
    let mut tags = vec!["pdbx_ordinal", "classification", "name", "version"];
    if write_all_fields {
        // Gemmi✔️✔️: loop.tags.insert(loop.tags.end(),
        // Gemmi✔️✔️:                  {"_software.date", "_software.description",
        // Gemmi✔️✔️:                   "_software.contact_author", "_software.contact_author_email"});
        tags.extend([
            "date",
            "description",
            "contact_author",
            "contact_author_email",
        ]);
    }
    let mut rows = Vec::with_capacity(structure.metadata.software.len());
    for (index, item) in structure.metadata.software.iter().enumerate() {
        // Gemmi✔️✔️: loop.add_values({
        // Gemmi✔️✔️:     std::to_string(++ordinal),
        // Gemmi✔️✔️:     cif::quote(software_classification_to_string(item.classification)),
        // Gemmi✔️✔️:     cif::quote(item.name),
        // Gemmi✔️✔️:     string_or_dot(item.version)});
        let mut row = vec![
            (index + 1).to_string(),
            quote_cif_value(software_classification_text(item.classification)),
            quote_cif_value(&item.name),
            string_or_dot(&item.version),
        ];
        if write_all_fields {
            // Gemmi✔️✔️: loop.add_values({
            // Gemmi✔️✔️:     string_or_qmark(item.date),
            // Gemmi✔️✔️:     string_or_qmark(item.description),
            // Gemmi✔️✔️:     string_or_qmark(item.contact_author),
            // Gemmi✔️✔️:     string_or_qmark(item.contact_author_email)});
            row.extend([
                string_or_qmark(&item.date),
                string_or_qmark(&item.description),
                string_or_qmark(&item.contact_author),
                string_or_qmark(&item.contact_author_email),
            ]);
        }
        rows.push(row);
    }
    add_mmcif_rows(block, "_software.", &tags, rows)
}

// BEGIN GEMMI CPP FUNCTION gemmi::update_mmcif_block
// Gemmi✔️✔️: if (st.models.empty())
// Gemmi✔️✔️:   return;
// END GEMMI CPP FUNCTION
pub(super) fn update_mmcif_block(
    structure: &BioStructure,
    block: &mut CifBlock,
    groups: MmcifOutputGroups,
) -> Result<(), BioWriteError> {
    if structure.models.is_empty() {
        return Ok(());
    }
    let entry_id = write_primary_mmcif_categories(structure, block, groups)?;
    write_secondary_mmcif_categories(structure, block, groups, &entry_id)
}

// BEGIN GEMMI CPP FUNCTION gemmi::make_mmcif_document
// Gemmi✔️✔️: cif::Document make_mmcif_document(const Structure& st, MmcifOutputGroups groups) {
// Gemmi✔️✔️:   cif::Document doc;
// Gemmi✔️✔️:   doc.blocks.resize(1);
// Gemmi✔️✔️:   update_mmcif_block(st, doc.blocks[0], groups);
// Gemmi✔️✔️:   return doc;
// Gemmi✔️✔️: }
// END GEMMI CPP FUNCTION
pub(super) fn make_mmcif_document(
    structure: &BioStructure,
    groups: MmcifOutputGroups,
) -> Result<CifDocument, BioWriteError> {
    let mut document = CifDocument {
        blocks: vec![CifBlock::default()],
    };
    update_mmcif_block(structure, &mut document.blocks[0], groups)?;
    Ok(document)
}

// BEGIN GEMMI CPP FUNCTION gemmi::make_mmcif_block
// Gemmi✔️✔️: cif::Block make_mmcif_block(const Structure& st, MmcifOutputGroups groups) {
// Gemmi✔️✔️:   cif::Block block;
// Gemmi✔️✔️:   update_mmcif_block(st, block, groups);
// Gemmi✔️✔️:   return block;
// Gemmi✔️✔️: }
// END GEMMI CPP FUNCTION
pub(super) fn make_mmcif_block(
    structure: &BioStructure,
    groups: MmcifOutputGroups,
) -> Result<CifBlock, BioWriteError> {
    let mut block = CifBlock::default();
    update_mmcif_block(structure, &mut block, groups)?;
    Ok(block)
}

// BEGIN GEMMI CPP FUNCTION gemmi::make_mmcif_headers
// Gemmi✔️✔️: cif::Block make_mmcif_headers(const Structure& st) {
// Gemmi✔️✔️:   MmcifOutputGroups groups(true);
// Gemmi✔️✔️:   groups.atoms = false;
// Gemmi✔️✔️:   return make_mmcif_block(st, groups);
// Gemmi✔️✔️: }
// END GEMMI CPP FUNCTION
pub(super) fn make_mmcif_headers(structure: &BioStructure) -> Result<CifBlock, BioWriteError> {
    let mut groups = MmcifOutputGroups::all(true);
    groups.atoms = false;
    make_mmcif_block(structure, groups)
}

// BEGIN GEMMI CPP FUNCTION gemmi::add_minimal_mmcif_data
// Gemmi✔️✔️: void add_minimal_mmcif_data(const Structure& st, cif::Block& block) {
// Gemmi✔️✔️:   cif::ItemSpan cell_span(block.items, "_cell.");
// Gemmi✔️✔️:   write_cell_parameters(st.cell, cell_span);
// Gemmi✔️✔️:   block.set_pair("_symmetry.space_group_name_H-M", cif::quote(st.spacegroup_hm));
// Gemmi✔️✔️:   write_ncs_oper(st, block);
// Gemmi✔️✔️:   add_cif_atoms(st, block, /*use_group_pdb=*/false, /*auth_all=*/false);
// Gemmi✔️✔️: }
// END GEMMI CPP FUNCTION
pub(super) fn add_minimal_mmcif_data(
    structure: &BioStructure,
    block: &mut CifBlock,
) -> Result<(), BioWriteError> {
    let default_crystal;
    let crystal = if let Some(crystal) = &structure.crystal {
        crystal
    } else {
        default_crystal = super::default_crystal_info();
        &default_crystal
    };
    write_cell_parameters(crystal.cell, block);
    block.set_pair(
        "_symmetry.space_group_name_H-M",
        quote_cif_value(crystal.spacegroup_hm.as_deref().unwrap_or("")),
    );
    write_ncs_oper(structure, block);
    add_cif_atoms(structure, block, false, false)
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::bio::{
        BioAssembly, BioAssemblyGenerator, BioAssemblyOperator, BioCisPep, BioConnection,
        BioConnectionType, BioDiffractionInfo, BioEntityDbRef, BioExperimentCrystalInfo,
        BioExperimentInfo, BioHelix, BioHelixClass, BioModRes, BioNcsOperator, BioRefinementInfo,
        BioReflectionsInfo, BioSheet, BioSheetStrand, BioSoftwareClassification, BioSoftwareItem,
        BioTlsGroup, BioTlsSelection, PdbChainId, PdbSeqId,
    };
    use crate::io::bio::cif::{CifEntry, CifLoop, CifToken};

    fn loop_with_tag<'a>(block: &'a CifBlock, tag: &str) -> &'a CifLoop {
        block
            .loops
            .iter()
            .find(|loop_| loop_.tags.iter().any(|candidate| candidate == tag))
            .unwrap_or_else(|| panic!("missing CIF loop containing {tag}"))
    }

    fn loop_rows(loop_: &CifLoop) -> Vec<Vec<&str>> {
        assert_eq!(loop_.values.len() % loop_.tags.len(), 0);
        loop_
            .values
            .chunks(loop_.tags.len())
            .map(|row| row.iter().map(|token| token.value.as_str()).collect())
            .collect()
    }

    fn parse_structure(text: &str) -> BioStructure {
        BioStructure::from_mmcif_str(text, "writer-test.cif").unwrap()
    }

    fn active_has_tag(block: &CifBlock, tag: &str) -> bool {
        block.entries.iter().any(|entry| match entry {
            CifEntry::Pair(index) => block.items[*index].tag == tag,
            CifEntry::Loop(index) => block.loops[*index].tags.iter().any(|value| value == tag),
            CifEntry::Erased => false,
        })
    }

    fn full_output_structure() -> BioStructure {
        let mut structure = parse_structure(
            r#"
data_full
_entry.id FULL
_cell.length_a 10
_cell.length_b 11
_cell.length_c 12
_cell.angle_alpha 90
_cell.angle_beta 90
_cell.angle_gamma 90
_symmetry.space_group_name_H-M 'P 1'
loop_
_entity.id
_entity.type
1 polymer
loop_
_entity_poly.entity_id
_entity_poly.type
_entity_poly.pdbx_strand_id
1 polypeptide(L) X
loop_
_entity_poly_seq.entity_id
_entity_poly_seq.num
_entity_poly_seq.mon_id
_entity_poly_seq.hetero
1 1 ALA n
1 2 'GLY,SER' y
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
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
ATOM 1 C C . ALA A 1 1 0 0 0 10 ALA X
ATOM 2 N N . GLY A 1 2 1 0 0 20 GLY X
"#,
        );
        structure.name = "full".to_string();
        structure.metadata.entry_id = Some("FULL".to_string());
        structure.metadata.received_initial_deposition_date = Some("2026-08-23".to_string());
        structure.metadata.authors = vec!["A. AUTHOR".to_string()];
        structure.metadata.title = Some("Complete writer fixture".to_string());
        structure.metadata.pdbx_keywords = Some("TEST".to_string());
        structure.metadata.keywords = Some("writer,metadata".to_string());
        structure.metadata.remark_300_detail = Some("biological detail".to_string());
        structure.entities[0].dbrefs.push(BioEntityDbRef {
            db_name: "UNP".to_string(),
            accession_code: "P00001".to_string(),
            id_code: "TEST_HUMAN".to_string(),
            seq_begin: Some(PdbSeqId {
                seq_num: 10,
                ins_code: None,
            }),
            seq_end: Some(PdbSeqId {
                seq_num: 20,
                ins_code: None,
            }),
            db_begin: Some(PdbSeqId {
                seq_num: 1,
                ins_code: None,
            }),
            db_end: Some(PdbSeqId {
                seq_num: 2,
                ins_code: None,
            }),
            label_seq_begin: Some(1),
            label_seq_end: Some(2),
            ..BioEntityDbRef::default()
        });
        structure.metadata.experiments = vec![BioExperimentInfo {
            method: "X-RAY DIFFRACTION".to_string(),
            number_of_crystals: Some(1),
            unique_reflections: Some(100),
            diffraction_ids: vec!["D1".to_string()],
            reflections: BioReflectionsInfo {
                resolution_high: Some(1.5),
                completeness: Some(99.0),
                ..BioReflectionsInfo::default()
            },
            ..BioExperimentInfo::default()
        }];
        structure.metadata.experiment_crystals = vec![BioExperimentCrystalInfo {
            id: "C1".to_string(),
            description: "test crystal".to_string(),
            ph: Some(7.0),
            diffractions: vec![BioDiffractionInfo {
                id: "D1".to_string(),
                temperature: Some(100.0),
                source: "synchrotron".to_string(),
                detector: "pixel".to_string(),
                ..BioDiffractionInfo::default()
            }],
            ..BioExperimentCrystalInfo::default()
        }];
        structure.metadata.refinement = vec![BioRefinementInfo {
            id: "X-RAY".to_string(),
            resolution_high: Some(1.5),
            reflection_count: Some(100),
            work_set_count: Some(95),
            tls_groups: vec![BioTlsGroup {
                id: "1".to_string(),
                origin: [1.0, 2.0, 3.0],
                selections: vec![BioTlsSelection {
                    chain: "X".to_string(),
                    res_begin: Some(PdbSeqId {
                        seq_num: 10,
                        ins_code: None,
                    }),
                    res_end: Some(PdbSeqId {
                        seq_num: 20,
                        ins_code: Some(b'A'),
                    }),
                    details: "chain selection".to_string(),
                }],
                ..BioTlsGroup::default()
            }],
            ..BioRefinementInfo::default()
        }];
        structure.metadata.software = vec![BioSoftwareItem {
            name: "REFMAC".to_string(),
            version: "5.8".to_string(),
            date: "2026-08-23".to_string(),
            classification: BioSoftwareClassification::Refinement,
            ..BioSoftwareItem::default()
        }];
        if let Some(crystal) = &mut structure.crystal {
            crystal.explicit_matrices = true;
        }
        structure.has_origx = true;
        structure.origx.vec = [1.0, 0.0, 0.0];
        structure.ncs_operators = vec![BioNcsOperator {
            id: "N1".to_string(),
            given: true,
            ..BioNcsOperator::default()
        }];
        let ala = BioAtomAddress {
            chain_name: "X".to_string(),
            seq_id: Some(PdbSeqId {
                seq_num: 10,
                ins_code: None,
            }),
            residue_name: "ALA".to_string(),
            atom_name: "C".to_string(),
            ..BioAtomAddress::default()
        };
        let gly = BioAtomAddress {
            chain_name: "X".to_string(),
            seq_id: Some(PdbSeqId {
                seq_num: 20,
                ins_code: None,
            }),
            residue_name: "GLY".to_string(),
            atom_name: "N".to_string(),
            ..BioAtomAddress::default()
        };
        structure.helices = vec![BioHelix {
            start: ala.clone(),
            end: gly.clone(),
            helix_class: BioHelixClass::RAlpha,
            length: 2,
        }];
        structure.sheets = vec![BioSheet {
            name: "S1".to_string(),
            strands: vec![
                BioSheetStrand {
                    start: ala.clone(),
                    end: ala.clone(),
                    ..BioSheetStrand::default()
                },
                BioSheetStrand {
                    start: gly.clone(),
                    end: gly.clone(),
                    hbond_atom1: ala.clone(),
                    hbond_atom2: gly.clone(),
                    sense: 1,
                    ..BioSheetStrand::default()
                },
            ],
        }];
        structure.mod_residues = vec![BioModRes {
            chain_name: "X".to_string(),
            res_id: PdbSeqId {
                seq_num: 20,
                ins_code: None,
            },
            residue_name: "MSE".to_string(),
            parent_comp_id: "MET".to_string(),
            mod_id: "MOD1".to_string(),
            details: "selenomethionine".to_string(),
        }];
        structure.assemblies = vec![BioAssembly {
            name: "1".to_string(),
            author_determined: true,
            generators: vec![BioAssemblyGenerator {
                subchains: vec!["A".to_string()],
                operators: vec![BioAssemblyOperator::default()],
                ..BioAssemblyGenerator::default()
            }],
            ..BioAssembly::default()
        }];
        structure.connections = vec![BioConnection {
            name: "conn1".to_string(),
            type_: BioConnectionType::Covale,
            partner1: ala.clone(),
            partner2: gly.clone(),
            ..BioConnection::default()
        }];
        structure.cispeps = vec![BioCisPep {
            partner_c: ala,
            partner_n: gly,
            model_num: 1,
            reported_angle: Some(0.0),
            ..BioCisPep::default()
        }];
        structure
    }

    pub(crate) fn check_mmcif_output_group_all_categories_emit_and_replace_prepopulated_categories()
    {
        let structure = full_output_structure();
        let mut block = CifBlock::default();
        for tag in [
            "_struct_asym.legacy",
            "_atom_type.legacy",
            "_software.legacy",
        ] {
            block.push_pair(
                tag.to_string(),
                CifToken {
                    value: "old".to_string(),
                    line_number: 0,
                },
            );
        }
        block.push_pair(
            "_entry.id".to_string(),
            CifToken {
                value: "OLD".to_string(),
                line_number: 0,
            },
        );

        let groups = MmcifOutputGroups::all(true);
        let entry_id = write_primary_mmcif_categories(&structure, &mut block, groups).unwrap();
        write_secondary_mmcif_categories(&structure, &mut block, groups, &entry_id).unwrap();

        for tag in [
            "_entry.id",
            "_pdbx_database_status.entry_id",
            "_audit_author.name",
            "_cell.length_a",
            "_symmetry.space_group_name_H-M",
            "_entity.id",
            "_entity_poly.entity_id",
            "_struct_ref.id",
            "_struct_ref_seq.align_id",
            "_chem_comp.id",
            "_exptl.entry_id",
            "_exptl_crystal.id",
            "_exptl_crystal_grow.crystal_id",
            "_diffrn.id",
            "_diffrn_detector.diffrn_id",
            "_diffrn_radiation.diffrn_id",
            "_diffrn_source.diffrn_id",
            "_reflns.entry_id",
            "_refine.entry_id",
            "_struct.title",
            "_struct_keywords.entry_id",
            "_struct_ncs_oper.id",
            "_struct_asym.id",
            "_database_PDB_matrix.entry_id",
            "_struct_conf.id",
            "_struct_conf_type.id",
            "_struct_sheet.id",
            "_struct_sheet_order.sheet_id",
            "_struct_sheet_range.sheet_id",
            "_pdbx_struct_sheet_hbond.sheet_id",
            "_struct_biol.id",
            "_pdbx_struct_assembly.id",
            "_pdbx_struct_assembly_gen.assembly_id",
            "_pdbx_struct_oper_list.id",
            "_struct_conn.id",
            "_struct_conn_type.id",
            "_struct_mon_prot_cis.pdbx_id",
            "_pdbx_struct_mod_residue.id",
            "_atom_sites.entry_id",
            "_atom_type.symbol",
            "_entity_poly_seq.entity_id",
            "_atom_site.id",
            "_pdbx_refine_tls.id",
            "_pdbx_refine_tls_group.id",
            "_software.pdbx_ordinal",
        ] {
            assert!(
                active_has_tag(&block, tag),
                "missing active writer tag {tag}"
            );
        }
        for tag in [
            "_struct_asym.legacy",
            "_atom_type.legacy",
            "_software.legacy",
        ] {
            assert!(!active_has_tag(&block, tag), "stale category item {tag}");
        }
        assert_eq!(
            loop_rows(loop_with_tag(&block, "_struct_conf.id"))[0][15],
            "2"
        );
        assert_eq!(
            loop_rows(loop_with_tag(
                &block,
                "_pdbx_struct_mod_residue.auth_comp_id"
            ))[0],
            [
                "1",
                "X",
                "20",
                "?",
                "MSE",
                "MSE",
                "MET",
                "selenomethionine",
                "MOD1"
            ]
        );
        assert_eq!(
            loop_rows(loop_with_tag(&block, "_pdbx_refine_tls_group.id"))[0],
            [
                "1",
                "1",
                "X-RAY",
                "X",
                "10",
                "?",
                "X",
                "20",
                "A",
                "'chain selection'"
            ]
        );
    }

    pub(crate) fn check_mmcif_output_group_all_disabled_preserves_prepopulated_block_exactly() {
        let structure = full_output_structure();
        let mut block = CifBlock::default();
        block.name = "existing".to_string();
        block.push_pair(
            "_entry.id".to_string(),
            CifToken {
                value: "EXISTING".to_string(),
                line_number: 7,
            },
        );
        block.push_pair(
            "_software.name".to_string(),
            CifToken {
                value: "legacy".to_string(),
                line_number: 8,
            },
        );
        let expected = block.clone();
        let groups = MmcifOutputGroups::all(false);

        let entry_id = write_primary_mmcif_categories(&structure, &mut block, groups).unwrap();
        write_secondary_mmcif_categories(&structure, &mut block, groups, &entry_id).unwrap();

        assert_eq!(entry_id, "EXISTING");
        assert_eq!(block, expected);
    }

    #[test]
    fn mmcif_cell_writer_replaces_existing_values_and_preserves_category_position() {
        let mut block = CifBlock::default();
        block.push_pair(
            "_entry.id".to_string(),
            CifToken {
                value: "X".to_string(),
                line_number: 0,
            },
        );
        block.push_pair(
            "_cell.length_a".to_string(),
            CifToken {
                value: "old".to_string(),
                line_number: 0,
            },
        );
        block.push_pair(
            "_struct.title".to_string(),
            CifToken {
                value: "title".to_string(),
                line_number: 0,
            },
        );

        write_cell_parameters(
            CrystalCell {
                a: 10.25,
                b: 20.5,
                c: 30.75,
                alpha: 90.0,
                beta: 91.5,
                gamma: 120.0,
            },
            &mut block,
        );

        let ordered_tags = block
            .entries
            .iter()
            .filter_map(|entry| match entry {
                CifEntry::Pair(index) => Some(block.items[*index].tag.as_str()),
                CifEntry::Loop(_) | CifEntry::Erased => None,
            })
            .collect::<Vec<_>>();
        assert_eq!(
            ordered_tags,
            [
                "_entry.id",
                "_cell.length_a",
                "_cell.length_b",
                "_cell.length_c",
                "_cell.angle_alpha",
                "_cell.angle_beta",
                "_cell.angle_gamma",
                "_struct.title",
            ]
        );
        assert_eq!(block.items[1].value.value, "10.25");
    }

    #[test]
    fn mmcif_atom_writer_emits_single_model_base_columns_and_group_policy() {
        let structure = parse_structure(
            r#"
data_atoms
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
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
ATOM 17 C CA . ALA A 7 ? 11.104 13.207 9.9 1 20 ? 70 ALA X CA 4
HETATM 99 Cl CL . LIG B . ? -2.5 0 8.25 0.5 31.125 -1 80 LIG Y CL 4
"#,
        );
        let mut block = CifBlock::default();

        add_cif_atoms(&structure, &mut block, true, false).unwrap();

        let atom_site = loop_with_tag(&block, "_atom_site.id");
        assert_eq!(
            atom_site.tags,
            [
                "_atom_site.group_PDB",
                "_atom_site.id",
                "_atom_site.type_symbol",
                "_atom_site.label_atom_id",
                "_atom_site.label_alt_id",
                "_atom_site.label_comp_id",
                "_atom_site.label_asym_id",
                "_atom_site.label_entity_id",
                "_atom_site.label_seq_id",
                "_atom_site.pdbx_PDB_ins_code",
                "_atom_site.Cartn_x",
                "_atom_site.Cartn_y",
                "_atom_site.Cartn_z",
                "_atom_site.occupancy",
                "_atom_site.B_iso_or_equiv",
                "_atom_site.pdbx_formal_charge",
                "_atom_site.auth_seq_id",
                "_atom_site.auth_asym_id",
                "_atom_site.pdbx_PDB_model_num",
            ]
        );
        assert_eq!(
            loop_rows(atom_site),
            [
                vec![
                    "ATOM", "1", "C", "CA", ".", "ALA", "A", ".", "7", "?", "11.104", "13.207",
                    "9.9", "1", "20", "?", "70", "X", "4",
                ],
                vec![
                    "HETATM", "2", "CL", "CL", ".", "LIG", "B", ".", ".", "?", "-2.5", "0", "8.25",
                    "0.5", "31.125", "-1", "80", "Y", "4",
                ],
            ]
        );
        assert!(block.loops.iter().all(|loop_| {
            !loop_
                .tags
                .iter()
                .any(|tag| tag.starts_with("_atom_site_anisotrop."))
        }));
    }

    #[test]
    fn mmcif_atom_writer_emits_multi_model_optional_and_anisotropic_columns() {
        let structure = parse_structure(
            r#"
data_atoms
loop_
_entity.id
_entity.type
1 polymer
2 water
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
B 2
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
1 ATOM C CA B ALA A 1 7 A 1.25 2.5 3.75 1 11 1 70 ALA X CA 3 calc 12 0.25
2 HETATM O O . HOH B 2 8 ? 4 5 6 0.5 22 ? 80 HOH Y O 9 dum ? 0.75
loop_
_atom_site_anisotrop.id
_atom_site_anisotrop.U[1][1]
_atom_site_anisotrop.U[2][2]
_atom_site_anisotrop.U[3][3]
_atom_site_anisotrop.U[1][2]
_atom_site_anisotrop.U[1][3]
_atom_site_anisotrop.U[2][3]
1 1.25 2.5 3.75 4.5 5.5 6.5
"#,
        );
        let mut block = CifBlock::default();

        add_cif_atoms(&structure, &mut block, false, true).unwrap();

        let atom_site = loop_with_tag(&block, "_atom_site.id");
        assert!(
            !atom_site
                .tags
                .iter()
                .any(|tag| tag == "_atom_site.group_PDB")
        );
        assert!(
            atom_site
                .tags
                .iter()
                .any(|tag| tag == "_atom_site.auth_atom_id")
        );
        assert!(
            atom_site
                .tags
                .iter()
                .any(|tag| tag == "_atom_site.auth_comp_id")
        );
        assert_eq!(
            &atom_site.tags[atom_site.tags.len() - 3..],
            [
                "_atom_site.calc_flag",
                "_atom_site.pdbx_tls_group_id",
                "_atom_site.ccp4_deuterium_fraction",
            ]
        );
        let atom_rows = loop_rows(atom_site);
        assert_eq!(atom_rows.len(), 2);
        assert_eq!(atom_rows[0][3], "B");
        assert_eq!(atom_rows[0][14], "1");
        assert_eq!(&atom_rows[0][15..17], ["CA", "ALA"]);
        assert_eq!(&atom_rows[0][19..], ["3", "c", "12", "0.25"]);
        assert_eq!(atom_rows[1][14], "?");
        assert_eq!(&atom_rows[1][19..], ["9", "dum", "?", "0.75"]);

        let anisotrop = loop_with_tag(&block, "_atom_site_anisotrop.id");
        assert_eq!(
            anisotrop.tags,
            [
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
        assert_eq!(
            loop_rows(anisotrop),
            [vec![
                "1", "C", "1.25", "2.5", "3.75", "4.5", "5.5", "6.5", "3"
            ]]
        );
    }

    #[test]
    fn mmcif_atom_writer_preserves_gemmi_null_sentinels_and_zero_trace_aniso_rule() {
        let mut structure = parse_structure(
            r#"
data_atoms
loop_
_entity.id
_entity.type
1 polymer
loop_
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
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
1 C CA . ALA A 1 1 1 2 3 1 A 1
"#,
        );
        structure.models[0].source_model_number = None;
        structure.chains[0].source.auth_chain_id = Some(PdbChainId([0; 4], 0));
        structure.residues[0].source.subchain_id = None;
        structure.residues[0].source.label_entity_id = None;
        structure.residues[0].source.label_seq_id = None;
        structure.residues[0].source.seq_id = None;
        structure.atoms[0].occupancy = None;
        structure.atoms[0].b_iso = None;
        structure.atoms[0].formal_charge = Some(0);
        structure.atoms[0].calc_flag = BioCalcFlag::NoHydrogen;
        structure.atoms[0].tls_group_id = Some(-1);
        structure.atoms[0].anisou = Some([1.0, -1.0, 0.0, 4.0, 5.0, 6.0]);

        let mut block = CifBlock::default();
        block.init_mmcif_loop("_atom_site_anisotrop.", &["id", "U[1][1]"]);
        add_cif_atoms(&structure, &mut block, false, false).unwrap();

        let atom_site = loop_with_tag(&block, "_atom_site.id");
        assert!(
            !atom_site
                .tags
                .iter()
                .any(|tag| tag == "_atom_site.calc_flag")
        );
        assert!(
            !atom_site
                .tags
                .iter()
                .any(|tag| tag == "_atom_site.pdbx_tls_group_id")
        );
        let rows = loop_rows(atom_site);
        assert_eq!(rows[0][5], ".");
        assert_eq!(rows[0][6], ".");
        assert_eq!(rows[0][7], ".");
        assert_eq!(rows[0][8], "?");
        assert_eq!(rows[0][12], "NaN");
        assert_eq!(rows[0][13], "NaN");
        assert_eq!(rows[0][14], "?");
        assert_eq!(rows[0][15], "?");
        assert_eq!(rows[0][16], "''");
        assert_eq!(rows[0][17], "0");
        assert!(block.entries.iter().all(|entry| {
            match entry {
                CifEntry::Loop(index) => !block.loops[*index]
                    .tags
                    .iter()
                    .any(|tag| tag.starts_with("_atom_site_anisotrop.")),
                CifEntry::Pair(_) | CifEntry::Erased => true,
            }
        }));
    }

    #[test]
    fn mmcif_structural_category_writer_assemblies_deduplicate_operators_and_expand_chains() {
        let mut structure = parse_structure(
            r#"
data_assembly
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
ATOM 1 C CA . ALA A 1 0 0 0 10 ALA X
"#,
        );
        let identity = BioAssemblyOperator {
            name: "source-identity".to_string(),
            ..BioAssemblyOperator::default()
        };
        let duplicate_identity = BioAssemblyOperator {
            name: "source-duplicate".to_string(),
            type_: "ignored duplicate type".to_string(),
            ..BioAssemblyOperator::default()
        };
        let translated = BioAssemblyOperator {
            name: "source-translation".to_string(),
            transform: BioTransform {
                mat: identity_transform().mat,
                vec: [10.0, 0.0, 0.0],
            },
            ..BioAssemblyOperator::default()
        };
        structure.assemblies = vec![BioAssembly {
            name: "1".to_string(),
            author_determined: true,
            software_determined: true,
            oligomeric_details: "DIMERIC".to_string(),
            software_name: "writer-tool".to_string(),
            absa: Some(12.5),
            ssa: Some(f64::NAN),
            generators: vec![
                BioAssemblyGenerator {
                    subchains: vec!["A".to_string(), "B".to_string()],
                    operators: vec![identity, duplicate_identity],
                    ..BioAssemblyGenerator::default()
                },
                BioAssemblyGenerator {
                    chains: vec!["X".to_string()],
                    operators: vec![translated],
                    ..BioAssemblyGenerator::default()
                },
            ],
            ..BioAssembly::default()
        }];
        let mut block = CifBlock::default();

        write_assemblies(&structure, &mut block).unwrap();

        assert_eq!(
            loop_rows(loop_with_tag(&block, "_pdbx_struct_assembly.id")),
            [vec![
                "1",
                "author_and_software_defined_assembly",
                "writer-tool",
                "dimeric",
                "2",
            ]]
        );
        assert_eq!(
            loop_rows(loop_with_tag(&block, "_pdbx_struct_assembly_prop.biol_id")),
            [vec!["1", "'ABSA (A^2)'", "12.5"]]
        );
        assert_eq!(
            loop_rows(loop_with_tag(
                &block,
                "_pdbx_struct_assembly_gen.assembly_id"
            )),
            [vec!["1", "1,1", "A,B"], vec!["1", "2", "A"]]
        );
        let operator_rows = loop_rows(loop_with_tag(&block, "_pdbx_struct_oper_list.id"));
        assert_eq!(operator_rows.len(), 2);
        assert_eq!(&operator_rows[0][..2], ["1", "'identity operation'"]);
        assert_eq!(
            &operator_rows[1][..2],
            ["2", "'crystal symmetry operation'"]
        );
        assert_eq!(operator_rows[1][5], "10");
    }

    #[test]
    fn mmcif_structural_category_writer_ncs_inserts_only_missing_identity() {
        let mut structure = BioStructure::default();
        structure.ncs_oper_identity_id = Some("I".to_string());
        structure.ncs_operators = vec![BioNcsOperator {
            id: "N".to_string(),
            given: false,
            transform: BioTransform {
                mat: identity_transform().mat,
                vec: [1.0, 2.0, 3.0],
            },
        }];
        let mut block = CifBlock::default();

        write_ncs_oper(&structure, &mut block);

        let rows = loop_rows(loop_with_tag(&block, "_struct_ncs_oper.id"));
        assert_eq!(rows.len(), 2);
        assert_eq!(&rows[0][..2], ["I", "given"]);
        assert_eq!(&rows[1][..2], ["N", "generate"]);
        assert_eq!(&rows[1][5..=5], ["1"]);
        assert_eq!(&rows[1][9..=9], ["2"]);
        assert_eq!(&rows[1][13..=13], ["3"]);

        structure.ncs_operators.insert(
            0,
            BioNcsOperator {
                id: "I".to_string(),
                given: true,
                transform: identity_transform(),
            },
        );
        let mut duplicate_block = CifBlock::default();
        write_ncs_oper(&structure, &mut duplicate_block);
        assert_eq!(
            loop_rows(loop_with_tag(&duplicate_block, "_struct_ncs_oper.id"))
                .iter()
                .filter(|row| row[0] == "I")
                .count(),
            1
        );

        structure.ncs_operators.clear();
        let mut empty_block = CifBlock::default();
        write_ncs_oper(&structure, &mut empty_block);
        assert!(empty_block.loops.is_empty());
    }

    #[test]
    fn mmcif_structural_category_writer_connections_use_source_lookup_and_nearest_image() {
        let mut structure = parse_structure(
            r#"
data_conn
_cell.length_a 10
_cell.length_b 10
_cell.length_c 10
_cell.angle_alpha 90
_cell.angle_beta 90
_cell.angle_gamma 90
_symmetry.space_group_name_H-M 'P 1'
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
ATOM 1 C C . ALA A 1 0 0 0 1 ALA X
ATOM 2 N N A GLY A 2 9 0 0 2 GLY X
"#,
        );
        let connection = BioConnection {
            name: "hydrogen-1".to_string(),
            type_: BioConnectionType::Hydrog,
            partner1: BioAtomAddress {
                chain_name: "X".to_string(),
                seq_id: Some(PdbSeqId {
                    seq_num: 1,
                    ins_code: None,
                }),
                residue_name: "ALA".to_string(),
                atom_name: "C".to_string(),
                altloc: Some(AltLocLabel(b'B')),
            },
            partner2: BioAtomAddress {
                chain_name: "X".to_string(),
                seq_id: Some(PdbSeqId {
                    seq_num: 2,
                    ins_code: None,
                }),
                residue_name: "GLY".to_string(),
                atom_name: "N".to_string(),
                altloc: Some(AltLocLabel(b'A')),
            },
            asu: BioAsu::Any,
            link_id: "LINK-H".to_string(),
            ..BioConnection::default()
        };
        let mut missing = connection.clone();
        missing.name = "omitted".to_string();
        missing.type_ = BioConnectionType::Covale;
        missing.partner2.residue_name = "MISSING".to_string();
        structure.connections = vec![connection, missing];
        let mut block = CifBlock::default();

        write_struct_conn(&structure, &mut block).unwrap();

        let connection_loop = loop_with_tag(&block, "_struct_conn.id");
        assert_eq!(connection_loop.tags.len(), 23);
        assert_eq!(
            loop_rows(connection_loop),
            [vec![
                "hydrogen-1",
                "hydrog",
                "A",
                "ALA",
                "1",
                "C",
                "?",
                "X",
                "1",
                "?",
                "1_555",
                "A",
                "GLY",
                "2",
                "N",
                "A",
                "X",
                "2",
                "?",
                "1_455",
                "?",
                "1.0000",
                "LINK-H",
            ]]
        );
        assert_eq!(
            loop_rows(loop_with_tag(&block, "_struct_conn_type.id")),
            [vec!["hydrog"]]
        );
    }

    #[test]
    fn mmcif_structural_category_writer_cispeps_select_models_and_omit_missing_partners() {
        let mut structure = parse_structure(
            r#"
data_cis
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
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA . ALA A 1 0 0 0 1 ALA X 1
ATOM 2 C C  . ALA A 10 0 0 0 10 ALA X 2
ATOM 3 N N  . GLY A 11 1 0 0 11 GLY X 2
"#,
        );
        let valid = BioCisPep {
            partner_c: BioAtomAddress {
                chain_name: "X".to_string(),
                seq_id: Some(PdbSeqId {
                    seq_num: 10,
                    ins_code: None,
                }),
                residue_name: "ALA".to_string(),
                ..BioAtomAddress::default()
            },
            partner_n: BioAtomAddress {
                chain_name: "X".to_string(),
                seq_id: Some(PdbSeqId {
                    seq_num: 11,
                    ins_code: None,
                }),
                residue_name: "GLY".to_string(),
                ..BioAtomAddress::default()
            },
            model_num: 2,
            only_altloc: Some(AltLocLabel(b'B')),
            reported_angle: Some(f64::NAN),
        };
        let mut missing_model = valid.clone();
        missing_model.model_num = 9;
        let mut missing_partner = valid.clone();
        missing_partner.partner_n.residue_name = "MISSING".to_string();
        structure.cispeps = vec![missing_model, valid, missing_partner];
        let mut block = CifBlock::default();

        write_cispeps(&structure, &mut block).unwrap();

        assert_eq!(
            loop_rows(loop_with_tag(&block, "_struct_mon_prot_cis.pdbx_id")),
            [vec![
                "1", "2", "A", "10", "ALA", "X", "10", "?", "A", "11", "GLY", "X", "11", "?", "B",
                "?",
            ]]
        );

        structure.cispeps.retain(|cispep| cispep.model_num == 9);
        let mut omitted_block = CifBlock::default();
        write_cispeps(&structure, &mut omitted_block).unwrap();
        assert!(omitted_block.loops.is_empty());
    }
}
