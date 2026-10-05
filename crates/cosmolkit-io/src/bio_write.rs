//! IO-owned full Gemmi mmCIF writer over detached BIO values.

pub mod atoms;
mod author;
pub mod categories;
mod cell;
mod cell_category;
mod database_status;
mod entity;
mod full;
mod ncs;
mod origx;
pub mod tags;
mod title_keywords;
pub mod value;

use std::fmt;

use cosmolkit_bio::{BioStructureData, BioStructureError};

use crate::cif::CifReadError;

/// Full source-defined mmCIF category switches and canonical CIF layout controls.
///
/// The default enables every category and `group_pdb`; `auth_all` is false.
/// Switches are flattened into this single IO-owned parameter value.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BioMmcifWriteParams {
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

    /// Emit the `_atom_site.group_PDB` column (ATOM/HETATM records).
    pub group_pdb: bool,
    /// Emit `_atom_site.auth_atom_id`/`auth_comp_id` in addition to label ids.
    pub auth_all: bool,
    /// Write single-row loops as tag/value pairs.
    pub prefer_pairs: bool,
    /// Omit blank lines between categories (blocks still separated).
    pub compact: bool,
    /// Put `#` empty comments before/after categories.
    pub misuse_hash: bool,
    /// Column at which pair values start (0 = single space).
    pub align_pairs: u16,
    /// Max loop column width used for value alignment (0 = no alignment).
    pub align_loops: u16,
}

impl Default for BioMmcifWriteParams {
    fn default() -> Self {
        // Gemmi❗✔️:   explicit MmcifOutputGroups(bool all)
        // Gemmi❗✔️:     : atoms(all), block_name(all), entry(all), database_status(all),
        // Gemmi❗✔️:       author(all), cell(all), symmetry(all), entity(all), entity_poly(all),
        // Gemmi❗✔️:       struct_ref(all), chem_comp(all), exptl(all), diffrn(all),
        // Gemmi❗✔️:       reflns(all), refine(all), title_keywords(all), ncs(all),
        // Gemmi❗✔️:       struct_asym(all), origx(all), struct_conf(all), struct_sheet(all),
        // Gemmi❗✔️:       struct_biol(all), assembly(all), conn(all), cis(all), modres(all),
        // Gemmi❗✔️:       scale(all), atom_type(all), entity_poly_seq(all), tls(all),
        // Gemmi❗✔️:       software(all), group_pdb(all), auth_all(false) {}
        // Cost: fixed-size scalar construction, no allocation.
        Self {
            atoms: true,
            block_name: true,
            entry: true,
            database_status: true,
            author: true,
            cell: true,
            symmetry: true,
            entity: true,
            entity_poly: true,
            struct_ref: true,
            chem_comp: true,
            exptl: true,
            diffrn: true,
            reflns: true,
            refine: true,
            title_keywords: true,
            ncs: true,
            struct_asym: true,
            origx: true,
            struct_conf: true,
            struct_sheet: true,
            struct_biol: true,
            assembly: true,
            conn: true,
            cis: true,
            modres: true,
            scale: true,
            atom_type: true,
            entity_poly_seq: true,
            tls: true,
            software: true,
            group_pdb: true,
            auth_all: false,
            prefer_pairs: false,
            compact: false,
            misuse_hash: false,
            align_pairs: 0,
            align_loops: 0,
        }
    }
}

/// A structural invariant, CIF mutation, invalid text or file IO failure.
#[derive(Debug)]
pub enum BioMmcifWriteError {
    Structure(BioStructureError),
    Cif(CifReadError),
    InvalidText {
        field: &'static str,
    },
    Io(std::io::Error),
    /// A destination file could not be written; the path and the
    /// underlying `std::io::Error` are both retained.
    FileWrite {
        path: std::path::PathBuf,
        source: std::io::Error,
    },
}

impl fmt::Display for BioMmcifWriteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Structure(error) => write!(f, "BIO structure invalid for mmCIF writing: {error}"),
            Self::Cif(error) => write!(f, "mmCIF document error: {error}"),
            Self::InvalidText { field } => write!(f, "invalid UTF-8 in {field}"),
            Self::Io(error) => write!(f, "mmCIF IO error: {error}"),
            Self::FileWrite { path, source } => {
                write!(f, "failed to write mmCIF file {}: {source}", path.display())
            }
        }
    }
}

impl std::error::Error for BioMmcifWriteError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Structure(error) => Some(error),
            Self::Cif(error) => Some(error),
            Self::Io(error) => Some(error),
            Self::FileWrite { source, .. } => Some(source),
            Self::InvalidText { .. } => None,
        }
    }
}

impl From<BioStructureError> for BioMmcifWriteError {
    fn from(error: BioStructureError) -> Self {
        Self::Structure(error)
    }
}
impl From<CifReadError> for BioMmcifWriteError {
    fn from(error: CifReadError) -> Self {
        Self::Cif(error)
    }
}
impl From<std::io::Error> for BioMmcifWriteError {
    fn from(error: std::io::Error) -> Self {
        Self::Io(error)
    }
}

/// Build the complete selected-group document using the single canonical CIF owner.
pub(crate) fn make_mmcif_document(
    data: &BioStructureData,
    params: &BioMmcifWriteParams,
) -> Result<crate::cif::CifDocument, BioMmcifWriteError> {
    // Gemmi❗✔️: cif::Document make_mmcif_document(const Structure& st, MmcifOutputGroups groups) {
    // Gemmi❗✔️:   cif::Document doc;
    // Gemmi❗✔️:   doc.blocks.resize(1);
    // Gemmi❗✔️:   update_mmcif_block(st, doc.blocks[0], groups);
    // Gemmi❗✔️:   return doc;
    // Gemmi❗✔️: }
    // Behavior: preserve the source whole-function zero-model early return.
    // Cost: one document and category output storage; no reparse or cloned BIO data.
    let mut document = crate::cif::CifDocument::with_single_block("");
    if !data.models().is_empty() {
        let block = document.sole_block_mut()?;
        let id = full::write_primary_mmcif_categories(data, block, params)?;
        full::write_secondary_mmcif_categories(data, block, *params, &id)?;
    }
    Ok(document)
}

/// Serialize a detached BIO structure with all selected source categories.
///
pub fn bio_structure_to_mmcif_text(
    data: &BioStructureData,
    params: &BioMmcifWriteParams,
) -> Result<String, BioMmcifWriteError> {
    // Gemmi✔️✔️: cif::Document make_mmcif_document(const Structure& st, MmcifOutputGroups groups) {
    // Gemmi✔️✔️:   cif::Document doc;
    // Gemmi✔️✔️:   doc.blocks.resize(1);
    // Gemmi✔️✔️:   update_mmcif_block(st, doc.blocks[0], groups);
    // Gemmi✔️✔️:   return doc;
    // Gemmi✔️✔️: }
    // Gemmi✔️✔️: inline void write_cif_to_stream(std::ostream& os, const Document& doc,
    // Gemmi✔️✔️:                                 WriteOptions options=WriteOptions()) {
    // Gemmi✔️✔️:   bool first = true;
    // Gemmi✔️✔️:   for (const Block& block : doc.blocks) {
    // Gemmi✔️✔️:     if (!first)
    // Gemmi✔️✔️:       os.put('\n'); // extra blank line for readability
    // Gemmi✔️✔️:     write_cif_block_to_stream(os, block, options);
    // Gemmi✔️✔️:     first = false;
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️: }
    // Behavior: compose the profile document, then serialize it with the
    // c08/c09 pipeline driven by the params' five layout controls; the
    // returned String holds exactly the emitted bytes. The input is
    // borrowed and never mutated; group dispatch belongs to the composition owner.
    // Complexity: one document construction plus one serialization pass;
    // no reparse, document copy or double buffering beyond the returned
    // String.
    let document = make_mmcif_document(data, params)?;
    crate::cif::bio_coordinate_document_to_string(&document, params).map_err(BioMmcifWriteError::Io)
}

/// Write a detached BIO structure to a selected-group mmCIF file.
///
/// The full document is serialized in memory before the destination is
/// opened, so a serialization failure never creates or truncates the file.
/// A failure during the physical write itself surfaces as
/// [`BioMmcifWriteError::FileWrite`] retaining the path and the underlying
/// `std::io::Error`; no atomicity is promised for that late stage.
pub fn write_bio_structure_mmcif_file(
    data: &BioStructureData,
    params: &BioMmcifWriteParams,
    path: &std::path::Path,
) -> Result<(), BioMmcifWriteError> {
    // Behavior: serialize completely (composition + text) before any
    // filesystem effect, then write the bytes through std::fs::write which
    // creates or truncates the destination; IO failures keep the path and
    // the original error.
    // Complexity: one serialization pass plus one file write; no partial
    // buffered handle is held open across serialization.
    let text = bio_structure_to_mmcif_text(data, params)?;
    std::fs::write(path, text.as_bytes()).map_err(|source| BioMmcifWriteError::FileWrite {
        path: path.to_path_buf(),
        source,
    })
}

#[cfg(test)]
fn coordinate_test_params() -> BioMmcifWriteParams {
    // Proposal: retain the original coordinate-profile expected bytes with
    // explicit pinned group controls. Full defaults have separate native tests.
    BioMmcifWriteParams {
        atoms: true,
        block_name: true,
        entry: true,
        database_status: false,
        author: false,
        cell: false,
        symmetry: false,
        entity: false,
        entity_poly: false,
        struct_ref: false,
        chem_comp: false,
        exptl: false,
        diffrn: false,
        reflns: false,
        refine: false,
        title_keywords: false,
        ncs: false,
        struct_asym: false,
        origx: false,
        struct_conf: false,
        struct_sheet: false,
        struct_biol: false,
        assembly: false,
        conn: false,
        cis: false,
        modres: false,
        scale: false,
        atom_type: false,
        entity_poly_seq: false,
        tls: false,
        software: false,
        group_pdb: true,
        auth_all: false,
        prefer_pairs: false,
        compact: false,
        misuse_hash: false,
        align_pairs: 0,
        align_loops: 0,
    }
}

#[cfg(test)]
mod compose_tests {
    use super::{BioMmcifWriteParams, make_mmcif_document};
    use crate::cif::CifItem;
    use cosmolkit_bio::{
        AtomName, AtomSourceIds, BioAtomRow, BioCalcFlag, BioChainId, BioChainRow, BioConnection,
        BioCoordinateBlock, BioCoordinateFormat, BioModelId, BioModelRow, BioNcsOperator,
        BioResidueId, BioResidueRow, BioRowSpan, BioSiftsUnpResidue, BioStructureData,
        BioStructureParts, BioStructureSourceState, BioTransform, ChainKind, ChainSourceIds,
        EntityKind, ResidueInfoKind, ResidueName, ResidueSourceIds,
    };
    use cosmolkit_types::Element;
    use std::collections::BTreeMap;

    fn populated_parts() -> BioStructureParts {
        let mut parts = BioStructureParts {
            input_format: BioCoordinateFormat::Mmcif,
            ..empty_base()
        };
        parts.source_state = BioStructureSourceState {
            name: "1abc".to_string(),
            info: BTreeMap::new(),
            ..BioStructureSourceState::default()
        };
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
            BioRowSpan::new(0, 1).unwrap(),
            ResidueName::from_ascii(b"ALA").unwrap(),
            ResidueInfoKind::Aa,
            EntityKind::Polymer,
            None,
            None,
            ResidueSourceIds::default(),
            BioSiftsUnpResidue::default(),
        )];
        parts.atoms = vec![BioAtomRow::new(
            BioResidueId::new(0),
            AtomName::from_ascii(b"CA").unwrap(),
            Element::C,
            None,
            None,
            0,
            BioCalcFlag::NotSet,
            1.0,
            20.0,
            [0.05, 0.0, 0.0, 0.0, 0.0, 0.0],
            -1,
            0.0,
            AtomSourceIds::default(),
        )];
        parts.coordinates = BioCoordinateBlock::new(vec![[1.0, 2.0, 3.0]]);
        parts
    }

    fn empty_base() -> BioStructureParts {
        BioStructureParts {
            input_format: BioCoordinateFormat::Mmcif,
            models: Vec::new(),
            chains: Vec::new(),
            residues: Vec::new(),
            atoms: Vec::new(),
            entities: Vec::new(),
            connections: Vec::new(),
            cispeps: Vec::new(),
            mod_residues: Vec::new(),
            helices: Vec::new(),
            sheets: Vec::new(),
            metadata: cosmolkit_bio::BioMetadata::default(),
            source_state: BioStructureSourceState::default(),
            coordinates: BioCoordinateBlock::new(Vec::new()),
            crystal: None,
            ncs_operators: Vec::new(),
            assemblies: Vec::new(),
        }
    }

    #[test]
    fn bio_pdbscope_compose_empty_models_yields_single_empty_block() {
        // make_mmcif_document + the preamble early return: the resized block
        // keeps Gemmi's default empty name and receives no items at all.
        let data = BioStructureData::from_parts(empty_base()).unwrap();
        let document = make_mmcif_document(&data, &super::coordinate_test_params()).unwrap();
        assert_eq!(document.blocks().len(), 1);
        let block = document.sole_block().unwrap();
        assert_eq!(block.name(), "");
        assert!(block.items().is_empty());
    }

    #[test]
    fn bio_pdbscope_compose_exact_category_order_populated() {
        // Exactly: _entry.id pair first (preamble), then the _atom_site loop,
        // then the conditional _atom_site_anisotrop loop for the nonzero
        // tensor. Nothing else exists in the document.
        let data = BioStructureData::from_parts(populated_parts()).unwrap();
        let document = make_mmcif_document(&data, &super::coordinate_test_params()).unwrap();
        let block = document.sole_block().unwrap();
        assert_eq!(block.name(), "1abc");
        let tags: Vec<&str> = block
            .items()
            .iter()
            .map(|item| match item {
                CifItem::Pair(pair) => pair.tag(),
                CifItem::Loop(loop_) => loop_.tags().first().map(String::as_str).unwrap_or(""),
                CifItem::Frame(_) => "<frame>",
            })
            .collect();
        assert_eq!(
            tags,
            vec![
                "_entry.id",
                "_atom_site.group_PDB",
                "_atom_site_anisotrop.id"
            ]
        );
        assert_eq!(block.find_value("_entry.id").unwrap().raw(), "1abc");
        assert!(
            block
                .find_mmcif_category("_atom_site.")
                .unwrap()
                .is_present()
        );
        assert_eq!(
            block.find_mmcif_category("_atom_site.").unwrap().len(),
            1,
            "one atom row"
        );
    }

    #[test]
    fn bio_pdbscope_compose_metadata_rich_input_stays_coordinate_only() {
        // NCS/connection/refinement-heavy input must not change the emitted
        // category set, and the borrowed input is not mutated.
        let mut parts = populated_parts();
        parts.ncs_operators = vec![BioNcsOperator::new(
            "1".to_string(),
            true,
            BioTransform::default(),
        )];
        parts.connections = vec![BioConnection::default()];
        parts.metadata = cosmolkit_bio::BioMetadata {
            authors: vec!["Author, A.".to_string()],
            ..cosmolkit_bio::BioMetadata::default()
        };
        let data = BioStructureData::from_parts(parts).unwrap();
        let atoms_before = data.atoms().len();
        let coords_before = data.coordinates().positions().len();
        let document = make_mmcif_document(&data, &super::coordinate_test_params()).unwrap();
        let block = document.sole_block().unwrap();
        for forbidden in [
            "_struct_ncs_oper.",
            "_struct_conn.",
            "_audit_author.",
            "_cell.",
            "_symmetry.",
            "_refine.",
            "_pdbx_struct_assembly.",
        ] {
            assert!(
                !block.find_mmcif_category(forbidden).unwrap().is_present(),
                "{forbidden} must not be serialized"
            );
        }
        let categories: Vec<&str> = block
            .items()
            .iter()
            .map(|item| match item {
                CifItem::Pair(pair) => pair.tag(),
                CifItem::Loop(loop_) => loop_.tags().first().map(String::as_str).unwrap_or(""),
                CifItem::Frame(_) => "<frame>",
            })
            .collect();
        assert_eq!(
            categories,
            vec![
                "_entry.id",
                "_atom_site.group_PDB",
                "_atom_site_anisotrop.id"
            ]
        );
        assert_eq!(data.atoms().len(), atoms_before);
        assert_eq!(data.coordinates().positions().len(), coords_before);
        assert_eq!(data.source_state().name, "1abc");
    }
}

#[cfg(test)]
mod text_tests {
    use super::{BioMmcifWriteParams, bio_structure_to_mmcif_text};
    use cosmolkit_bio::{
        AtomName, AtomSourceIds, BioAtomRow, BioCalcFlag, BioChainId, BioChainRow,
        BioCoordinateBlock, BioCoordinateFormat, BioModelId, BioModelRow, BioResidueId,
        BioResidueRow, BioRowSpan, BioSiftsUnpResidue, BioStructureData, BioStructureParts,
        BioStructureSourceState, ChainKind, ChainSourceIds, EntityKind, PdbSeqId, ResidueInfoKind,
        ResidueName, ResidueSourceIds,
    };
    use cosmolkit_types::Element;
    use std::collections::BTreeMap;

    fn one_atom_parts(name: &str, entry_id: Option<&str>) -> BioStructureParts {
        let mut parts = BioStructureParts {
            input_format: BioCoordinateFormat::Mmcif,
            ..empty_base()
        };
        let mut info = BTreeMap::new();
        if let Some(id) = entry_id {
            info.insert("_entry.id".to_string(), id.to_string());
        }
        parts.source_state = BioStructureSourceState {
            name: name.to_string(),
            info,
            ..BioStructureSourceState::default()
        };
        parts.models = vec![BioModelRow::new(BioRowSpan::new(0, 1).unwrap(), Some(1))];
        parts.chains = vec![BioChainRow::new(
            BioModelId::new(0),
            None,
            BioRowSpan::new(0, 1).unwrap(),
            ChainKind::Protein,
            ChainSourceIds::new(
                Some(cosmolkit_bio::PdbChainId::from_ascii(b"A").unwrap()),
                Some("A".to_string()),
            ),
        )];
        parts.residues = vec![BioResidueRow::new(
            BioChainId::new(0),
            BioRowSpan::new(0, 1).unwrap(),
            ResidueName::from_ascii(b"ALA").unwrap(),
            ResidueInfoKind::Aa,
            EntityKind::Polymer,
            None,
            None,
            ResidueSourceIds::new(
                Some(PdbSeqId::new(1, None)),
                Some(1),
                None,
                Some("A".to_string()),
                None,
            )
            .unwrap(),
            BioSiftsUnpResidue::default(),
        )];
        parts.atoms = vec![BioAtomRow::new(
            BioResidueId::new(0),
            AtomName::from_ascii(b"CA").unwrap(),
            Element::C,
            None,
            None,
            0,
            BioCalcFlag::NotSet,
            1.0,
            20.0,
            [0.05, 0.0, 0.0, 0.0, 0.0, 0.0],
            -1,
            0.0,
            AtomSourceIds::default(),
        )];
        parts.coordinates = BioCoordinateBlock::new(vec![[1.0, 2.0, 3.0]]);
        parts
    }

    fn empty_base() -> BioStructureParts {
        BioStructureParts {
            input_format: BioCoordinateFormat::Mmcif,
            models: Vec::new(),
            chains: Vec::new(),
            residues: Vec::new(),
            atoms: Vec::new(),
            entities: Vec::new(),
            connections: Vec::new(),
            cispeps: Vec::new(),
            mod_residues: Vec::new(),
            helices: Vec::new(),
            sheets: Vec::new(),
            metadata: cosmolkit_bio::BioMetadata::default(),
            source_state: BioStructureSourceState::default(),
            coordinates: BioCoordinateBlock::new(Vec::new()),
            crystal: None,
            ncs_operators: Vec::new(),
            assemblies: Vec::new(),
        }
    }

    #[test]
    fn bio_pdbscope_text_default_exact_bytes() {
        // Every line below derives from the pinned write_out_pair /
        // write_out_loop / write_cif_block_to_stream layout rules and the
        // add_cif_atoms column profile, not from CK output.
        let data = BioStructureData::from_parts(one_atom_parts("1abc", None)).unwrap();
        let text = bio_structure_to_mmcif_text(&data, &super::coordinate_test_params()).unwrap();
        let expected = concat!(
            "data_1abc\n",
            "_entry.id 1abc\n",
            "\n",
            "loop_\n",
            "_atom_site.group_PDB\n",
            "_atom_site.id\n",
            "_atom_site.type_symbol\n",
            "_atom_site.label_atom_id\n",
            "_atom_site.label_alt_id\n",
            "_atom_site.label_comp_id\n",
            "_atom_site.label_asym_id\n",
            "_atom_site.label_entity_id\n",
            "_atom_site.label_seq_id\n",
            "_atom_site.pdbx_PDB_ins_code\n",
            "_atom_site.Cartn_x\n",
            "_atom_site.Cartn_y\n",
            "_atom_site.Cartn_z\n",
            "_atom_site.occupancy\n",
            "_atom_site.B_iso_or_equiv\n",
            "_atom_site.pdbx_formal_charge\n",
            "_atom_site.auth_seq_id\n",
            "_atom_site.auth_asym_id\n",
            "_atom_site.pdbx_PDB_model_num\n",
            "ATOM 1 C CA . ALA A . 1 ? 1 2 3 1 20 ? 1 A 1\n",
            "\n",
            "loop_\n",
            "_atom_site_anisotrop.id\n",
            "_atom_site_anisotrop.type_symbol\n",
            "_atom_site_anisotrop.U[1][1]\n",
            "_atom_site_anisotrop.U[2][2]\n",
            "_atom_site_anisotrop.U[3][3]\n",
            "_atom_site_anisotrop.U[1][2]\n",
            "_atom_site_anisotrop.U[1][3]\n",
            "_atom_site_anisotrop.U[2][3]\n",
            "1 C 0.05 0 0 0 0 0\n",
        );
        assert_eq!(text, expected);
    }

    #[test]
    fn bio_pdbscope_text_four_atom_column_options() {
        let data = BioStructureData::from_parts(one_atom_parts("1abc", None)).unwrap();
        for (group_pdb, auth_all, first_tag, has_auth, width) in [
            (true, false, "_atom_site.group_PDB", false, 19usize),
            (false, false, "_atom_site.id", false, 18),
            (true, true, "_atom_site.group_PDB", true, 21),
            (false, true, "_atom_site.id", true, 20),
        ] {
            let params = BioMmcifWriteParams {
                group_pdb,
                auth_all,
                ..super::coordinate_test_params()
            };
            let text = bio_structure_to_mmcif_text(&data, &params).unwrap();
            let lines: Vec<&str> = text.lines().collect();
            let loop_start = lines
                .iter()
                .position(|line| *line == "loop_")
                .expect("atom loop present");
            assert_eq!(lines[loop_start + 1], first_tag);
            assert_eq!(
                lines
                    .iter()
                    .filter(|line| **line == "_atom_site.auth_atom_id")
                    .count(),
                usize::from(has_auth)
            );
            assert_eq!(
                lines
                    .iter()
                    .filter(|line| **line == "_atom_site.auth_comp_id")
                    .count(),
                usize::from(has_auth)
            );
            // One atom row: exactly `width` space-separated values on the
            // single row line right after the tag block.
            let row_line = lines[loop_start + 1 + width];
            let tokens: Vec<&str> = row_line.split_whitespace().collect();
            assert_eq!(tokens.len(), width);
            assert_eq!(tokens.first(), Some(&if group_pdb { "ATOM" } else { "1" }));
            assert_eq!(tokens.last(), Some(&"1"));
            // Auth duplication duplicates the label atom/comp strings.
            if has_auth {
                assert_eq!((tokens[width - 5], tokens[width - 4]), ("CA", "ALA"));
            }
        }
    }

    #[test]
    fn bio_pdbscope_text_zero_model_and_quoted_entry_id() {
        // Zero-model early return: the composed document holds one fresh
        // empty block, so the text is exactly the mandatory data_ header.
        let empty = BioStructureData::from_parts(empty_base()).unwrap();
        assert_eq!(
            bio_structure_to_mmcif_text(&empty, &super::coordinate_test_params()).unwrap(),
            "data_\n"
        );
        // An entry id containing a newline is quoted as a semicolon text
        // field (cif::quote \n branch + write_out_pair is_text_field path).
        let data = BioStructureData::from_parts(one_atom_parts("n", Some("a\nb"))).unwrap();
        let text = bio_structure_to_mmcif_text(&data, &super::coordinate_test_params()).unwrap();
        assert!(text.starts_with("data_n\n_entry.id\n;a\nb\n;\n\nloop_\n"));
    }

    #[test]
    fn bio_pdbscope_file_exact_bytes_and_replacement() {
        let data = BioStructureData::from_parts(one_atom_parts("1abc", None)).unwrap();
        let directory = std::env::temp_dir();
        let path = directory.join(format!("ck-bio-pdb-file-{}.cif", std::process::id()));
        // A longer pre-existing sentinel proves full truncation, not append.
        std::fs::write(&path, b"sentinel content that is longer than the output\n").unwrap();
        super::super::write_bio_structure_mmcif_file(
            &data,
            &super::coordinate_test_params(),
            &path,
        )
        .unwrap();
        let written = std::fs::read(&path).unwrap();
        let text = String::from_utf8(written).unwrap();
        // The same independently derived default-profile bytes as the TEXT
        // unit (source layout rules, not CK output).
        let expected = concat!(
            "data_1abc\n",
            "_entry.id 1abc\n",
            "\n",
            "loop_\n",
            "_atom_site.group_PDB\n",
            "_atom_site.id\n",
            "_atom_site.type_symbol\n",
            "_atom_site.label_atom_id\n",
            "_atom_site.label_alt_id\n",
            "_atom_site.label_comp_id\n",
            "_atom_site.label_asym_id\n",
            "_atom_site.label_entity_id\n",
            "_atom_site.label_seq_id\n",
            "_atom_site.pdbx_PDB_ins_code\n",
            "_atom_site.Cartn_x\n",
            "_atom_site.Cartn_y\n",
            "_atom_site.Cartn_z\n",
            "_atom_site.occupancy\n",
            "_atom_site.B_iso_or_equiv\n",
            "_atom_site.pdbx_formal_charge\n",
            "_atom_site.auth_seq_id\n",
            "_atom_site.auth_asym_id\n",
            "_atom_site.pdbx_PDB_model_num\n",
            "ATOM 1 C CA . ALA A . 1 ? 1 2 3 1 20 ? 1 A 1\n",
            "\n",
            "loop_\n",
            "_atom_site_anisotrop.id\n",
            "_atom_site_anisotrop.type_symbol\n",
            "_atom_site_anisotrop.U[1][1]\n",
            "_atom_site_anisotrop.U[2][2]\n",
            "_atom_site_anisotrop.U[3][3]\n",
            "_atom_site_anisotrop.U[1][2]\n",
            "_atom_site_anisotrop.U[1][3]\n",
            "_atom_site_anisotrop.U[2][3]\n",
            "1 C 0.05 0 0 0 0 0\n",
        );
        assert_eq!(text, expected);
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn bio_pdbscope_file_errors_retain_path_and_source() {
        let data = BioStructureData::from_parts(one_atom_parts("1abc", None)).unwrap();
        // Nonexistent parent directory: no file is created and the error
        // keeps both the path and the underlying io::Error.
        let missing = std::path::PathBuf::from("/nonexistent-ck-bio-pdb-dir/out.cif");
        let error = super::super::write_bio_structure_mmcif_file(
            &data,
            &super::coordinate_test_params(),
            &missing,
        )
        .unwrap_err();
        match &error {
            super::BioMmcifWriteError::FileWrite { path, source } => {
                assert_eq!(path, &missing);
                assert_eq!(source.kind(), std::io::ErrorKind::NotFound);
            }
            other => panic!("expected FileWrite, got {other:?}"),
        }
        assert!(
            error
                .to_string()
                .contains("/nonexistent-ck-bio-pdb-dir/out.cif")
        );
        use std::error::Error as _;
        assert!(error.source().is_some());
        // A directory destination is also a retained FileWrite failure.
        let directory = std::env::temp_dir();
        let error = super::super::write_bio_structure_mmcif_file(
            &data,
            &super::coordinate_test_params(),
            &directory,
        )
        .unwrap_err();
        assert!(matches!(error, super::BioMmcifWriteError::FileWrite { .. }));
    }
}

#[cfg(test)]
mod scope_tests {
    use super::{BioMmcifWriteError as Error, BioMmcifWriteParams as Params};
    use crate::cif::{CifCheckLevel, read_cif_document};
    use cosmolkit_bio::BioStructureError;
    use std::error::Error as _;

    #[test]
    fn bio_pdbscope_scope_exact_seven_defaults_and_independence() {
        let defaults = Params::default();
        assert!(defaults.group_pdb);
        assert!(!defaults.auth_all);
        assert!(!defaults.prefer_pairs);
        assert!(!defaults.compact);
        assert!(!defaults.misuse_hash);
        assert_eq!(defaults.align_pairs, 0);
        assert_eq!(defaults.align_loops, 0);
        // Preserve the original seven layout/atom controls; category defaults
        // now follow the pinned full MmcifOutputGroups constructor.
        let all_set = Params {
            group_pdb: false,
            auth_all: true,
            prefer_pairs: true,
            compact: true,
            misuse_hash: true,
            align_pairs: 33,
            align_loops: 30,
            ..Params::default()
        };
        assert_ne!(all_set, defaults);
        // Each control flips independently of the others.
        let mut one = defaults;
        one.group_pdb = false;
        assert_eq!(
            (one.group_pdb, one.auth_all, one.prefer_pairs),
            (false, false, false)
        );
        let mut two = defaults;
        two.auth_all = true;
        assert_eq!(
            (two.group_pdb, two.auth_all, two.align_pairs),
            (true, true, 0)
        );
        let mut three = defaults;
        three.align_loops = 30;
        assert_eq!(
            (three.compact, three.misuse_hash, three.align_loops),
            (false, false, 30)
        );
    }

    #[test]
    fn bio_pdbscope_scope_error_chains_retained() {
        let io = Error::from(std::io::Error::new(
            std::io::ErrorKind::PermissionDenied,
            "source path",
        ));
        assert_eq!(io.source().unwrap().to_string(), "source path");
        assert!(matches!(io, Error::Io(_)));
        let cif = Error::from(
            read_cif_document("invalid", "original.cif", CifCheckLevel::Syntax).unwrap_err(),
        );
        assert!(cif.source().unwrap().to_string().contains("original.cif"));
        assert!(matches!(cif, Error::Cif(_)));
        let structure = Error::from(BioStructureError::RowIndexTooLarge { value: usize::MAX });
        assert!(structure.source().is_some());
        assert!(matches!(structure, Error::Structure(_)));
        assert!(Error::InvalidText { field: "metadata" }.source().is_none());
    }
}
