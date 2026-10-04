//! Detached mmCIF-to-BIO parsing helpers.

use super::bio_pdb::{
    PdbCcdAlias, PdbResidueKeyError, mapped_ccd_atom_address, mapped_ccd_residue_address,
    rename_ccd_sequence_tokens,
};
use crate::cif::{
    CifBlock, CifCheckLevel, CifDocument, CifReadError, CifReadErrorKind, CifRow, CifTable,
    CifValue, cif_as_char, cif_as_f64, cif_as_i32, cif_as_string, cif_is_null, read_cif_document,
};
use cosmolkit_bio::{
    AltLocLabel, AtomAddress, AtomName, AtomSourceIds, BioAssembly, BioAssemblyGenerator,
    BioAssemblyOperator, BioAssemblySpecialKind, BioAsu, BioAtomRow, BioCalcFlag, BioChainId,
    BioChainRow, BioCisPep, BioConnection, BioConnectionKind, BioCoordinateBlock,
    BioCoordinateFormat, BioCrystalCell, BioCrystalInfo, BioDiffractionInfo, BioEntityDbRef,
    BioEntityId, BioEntityRow, BioExperimentInfo, BioExperimentalCrystalInfo, BioHelix,
    BioMetadata, BioModRes, BioModelId, BioModelRow, BioNcsOperator, BioRefinementInfo,
    BioResidueId, BioResidueRow, BioRowSpan, BioSheet, BioSiftsUnpResidue,
    BioSoftwareClassification, BioSoftwareItem, BioStrand, BioStructureData, BioStructureError,
    BioStructureParts, BioStructureSourceState, BioTlsGroup, BioTlsSelection, BioTransform,
    ChainKind, ChainSourceIds, EntityKind, EntitySourceIds, PdbAtomSerial, PdbChainId, PdbSeqId,
    PolymerKind, ResidueAddress, ResidueName, ResidueSourceIds, find_residue_info,
    set_crystal_cell, set_crystal_fractional_transform, set_crystal_space_group_hm,
    setup_cell_images,
};
use cosmolkit_types::Element;
use std::collections::HashMap;
use std::error::Error;
use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BioMmcifReadStage {
    CifDocument,
    CoordinateBlock,
    CrystalCell,
    Refinement,
    Tls,
    Experimental,
    Reflections,
    Software,
    Ncs,
    FractionalTransform,
    Origx,
    AnisotropicU,
    AtomSites,
    EntitySequence,
    Helices,
    Sheets,
    Connections,
    CisPeptides,
    ModifiedResidues,
    Assemblies,
    SiftsUnp,
    CcdRestoration,
    Materialization,
    StructureValidation,
}

impl fmt::Display for BioMmcifReadStage {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let name = match self {
            Self::CifDocument => "CIF document",
            Self::CoordinateBlock => "coordinate block selection",
            Self::CrystalCell => "crystal cell",
            Self::Refinement => "refinement metadata",
            Self::Tls => "TLS metadata",
            Self::Experimental => "experimental metadata",
            Self::Reflections => "reflection metadata",
            Self::Software => "software metadata",
            Self::Ncs => "NCS operators",
            Self::FractionalTransform => "fractional transform",
            Self::Origx => "ORIGX transform",
            Self::AnisotropicU => "anisotropic displacement",
            Self::AtomSites => "atom sites",
            Self::EntitySequence => "entity and sequence metadata",
            Self::Helices => "helices",
            Self::Sheets => "sheets",
            Self::Connections => "connections",
            Self::CisPeptides => "cis peptides",
            Self::ModifiedResidues => "modified residues",
            Self::Assemblies => "assemblies",
            Self::SiftsUnp => "SIFTS UNP mapping",
            Self::CcdRestoration => "CCD code restoration",
            Self::Materialization => "BIO hierarchy materialization",
            Self::StructureValidation => "BIO structure validation",
        };
        f.write_str(name)
    }
}

#[derive(Debug)]
pub struct BioMmcifReadError {
    stage: BioMmcifReadStage,
    cause: Box<dyn Error>,
}

impl BioMmcifReadError {
    fn new(stage: BioMmcifReadStage, cause: impl Error + 'static) -> Self {
        Self {
            stage,
            cause: Box::new(cause),
        }
    }

    #[must_use]
    pub const fn stage(&self) -> BioMmcifReadStage {
        self.stage
    }
}

impl fmt::Display for BioMmcifReadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "mmCIF {} failed: {}", self.stage, self.cause)
    }
}

impl Error for BioMmcifReadError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        Some(self.cause.as_ref())
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum CoordinateBlockSelectionError {
    EmptyDocument { source: String },
    LaterCoordinateBlock { block_number: usize, source: String },
}

impl fmt::Display for CoordinateBlockSelectionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyDocument { source } => {
                write!(f, "coordinate mmCIF document has no first block: {source}")
            }
            Self::LaterCoordinateBlock {
                block_number,
                source,
            } => write!(
                f,
                "2+ blocks are ok if only the first one has coordinates;\n_atom_site in block #{block_number}: {source}"
            ),
        }
    }
}

impl std::error::Error for CoordinateBlockSelectionError {}

fn select_coordinate_block<'a>(
    blocks: &'a [CifBlock],
    source: &str,
) -> Result<&'a CifBlock, CoordinateBlockSelectionError> {
    // Gemmi✔️✔️: for (size_t i = 1; i < doc.blocks.size(); ++i)
    // Gemmi✔️✔️:   if (doc.blocks[i].has_tag("_atom_site.id"))
    // Gemmi✔️✔️:     fail("2+ blocks are ok if only the first one has coordinates;\n"
    // Gemmi✔️✔️:          "_atom_site in block #" + std::to_string(i+1) + ": " + doc.source);
    // Gemmi✔️✔️: Structure st = make_structure_from_block(doc.blocks.at(0));
    // Behavior review: scan later blocks in source order using canonical CIF
    // tag-presence semantics before selecting block zero. Empty direct
    // documents remain a distinct out-of-range failure; parser-produced empty
    // text fails earlier in the CIF grammar. Complexity review: O(B + I),
    // where B is block count and I is the total item scan performed by
    // CifBlock::has_tag; this follows Gemmi's ordered block/tag scan and does
    // not clone blocks or their values.
    for (block_index, block) in blocks.iter().enumerate().skip(1) {
        if block.has_tag("_atom_site.id") {
            return Err(CoordinateBlockSelectionError::LaterCoordinateBlock {
                block_number: block_index + 1,
                source: source.to_owned(),
            });
        }
    }

    blocks
        .first()
        .ok_or_else(|| CoordinateBlockSelectionError::EmptyDocument {
            source: source.to_owned(),
        })
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum MmcifCrystalInfoError {
    Cif(CifReadError),
    Bio(BioStructureError),
    StagingInvariant { stage: &'static str },
}

impl From<CifReadError> for MmcifCrystalInfoError {
    fn from(error: CifReadError) -> Self {
        Self::Cif(error)
    }
}

impl From<BioStructureError> for MmcifCrystalInfoError {
    fn from(error: BioStructureError) -> Self {
        Self::Bio(error)
    }
}

impl fmt::Display for MmcifCrystalInfoError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Cif(error) => fmt::Display::fmt(error, f),
            Self::Bio(error) => fmt::Display::fmt(error, f),
            Self::StagingInvariant { stage } => {
                write!(
                    f,
                    "internal mmCIF crystal staging invariant failed: {stage}"
                )
            }
        }
    }
}

impl std::error::Error for MmcifCrystalInfoError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Cif(error) => Some(error),
            Self::Bio(error) => Some(error),
            Self::StagingInvariant { .. } => None,
        }
    }
}

fn parse_mmcif_crystal_info(block: &CifBlock) -> Result<BioCrystalInfo, MmcifCrystalInfoError> {
    let mut crystal = parse_mmcif_crystal_cell_info(block)?;
    apply_mmcif_fractional_transform(block, &mut crystal)?;
    Ok(crystal)
}

fn parse_mmcif_crystal_cell_info(
    block: &CifBlock,
) -> Result<BioCrystalInfo, MmcifCrystalInfoError> {
    // Gemmi❗✔️: inline void set_cell_from_mmcif(cif::Block& block, UnitCell& cell,
    // Gemmi❗✔️:                                 bool mmcif=true) {
    // Gemmi❗✔️:   cif::Table tab = block.find((mmcif ? "_cell." : "_cell_"),
    // Gemmi❗✔️:                               {"length_a", "length_b", "length_c",
    // Gemmi❗✔️:                                "angle_alpha", "angle_beta", "angle_gamma"});
    // Gemmi❗✔️:   if (tab.ok()) {
    // Gemmi❗✔️:     auto c = tab.one();
    // Gemmi❗✔️:     if (!cif::is_null(c[0]) && !cif::is_null(c[1]) && !cif::is_null(c[2]))
    // Gemmi❗✔️:       cell.set(cif::as_number(c[0]), cif::as_number(c[1]), cif::as_number(c[2]),
    // Gemmi❗✔️:                cif::as_number(c[3]), cif::as_number(c[4]), cif::as_number(c[5]));
    // Gemmi❗✔️:   }
    // Gemmi❗✔️: }
    // Gemmi❗✔️: inline const std::string* find_spacegroup_hm_value(const cif::Block& block) {
    // Gemmi❗✔️:   const char* hm_tag = "_symmetry.space_group_name_H-M";
    // Gemmi❗✔️:   return block.find_value(hm_tag);
    // Gemmi❗✔️: }
    // Behavior review: fixed `_cell.` and H-M keys retain Gemmi's field
    // precedence; alternate names are not consulted. The default BIO crystal
    // begins with the same unit-cell and identity-transform values as Gemmi's
    // default UnitCell. Pair/loop tag selection delegates to the canonical CIF
    // owner. Ordinary numeric bit parity remains bounded by the CIF owner's
    // documented floating-conversion scope.
    // Complexity review: six fixed cell reads and one H-M lookup; no model
    // clone, row scan, or fractional-transform conversion is introduced.
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

    let cell_table = block.find(
        "_cell.",
        &[
            "length_a",
            "length_b",
            "length_c",
            "angle_alpha",
            "angle_beta",
            "angle_gamma",
        ],
    )?;
    if cell_table.is_present() {
        let row = cell_table.one()?;
        let lengths = [
            row.get(0).map(CifValue::raw).unwrap_or(""),
            row.get(1).map(CifValue::raw).unwrap_or(""),
            row.get(2).map(CifValue::raw).unwrap_or(""),
        ];
        if lengths.iter().all(|value| !cif_is_null(value)) {
            let mut numbers = [f64::NAN; 6];
            for (column, number) in numbers.iter_mut().enumerate() {
                let raw = row.get(column).map(CifValue::raw).unwrap_or("");
                *number = cif_as_f64(raw, f64::NAN)?;
            }
            set_crystal_cell(
                &mut crystal,
                BioCrystalCell {
                    a: numbers[0],
                    b: numbers[1],
                    c: numbers[2],
                    alpha: numbers[3],
                    beta: numbers[4],
                    gamma: numbers[5],
                },
            )?;
        }
    }

    let space_group_hm = block
        .find_value("_symmetry.space_group_name_H-M")
        .map(CifValue::raw)
        .unwrap_or("");
    set_crystal_space_group_hm(&mut crystal, cif_as_string(space_group_hm));

    Ok(crystal)
}

fn apply_mmcif_fractional_transform(
    block: &CifBlock,
    crystal: &mut BioCrystalInfo,
) -> Result<(), MmcifCrystalInfoError> {
    // Gemmi❗✔️:   cif::Table fract_tv = block.find("_atom_sites.fract_transf_",
    // Gemmi❗✔️:                                    transform_tags("matrix", "vector"));
    // Gemmi❗✔️:   if (fract_tv.length() > 0) {
    // Gemmi❗✔️:     Transform fract = get_transform_matrix(fract_tv[0]);
    // Gemmi❗✔️:     st.cell.set_matrices_from_fract(fract);
    // Gemmi❗✔️:   }
    // Behavior review: only the first row of a nonempty complete transform
    // table overrides the current fractional transform, after NCS processing
    // in the full reader. The canonical crystal transition retains Gemmi's
    // stored matrices and recalculates its inverse representation.
    // Complexity review: twelve fixed scalar conversions and one bounded
    // crystal update; no extra scan or structure clone.
    let fractional_tags = mmcif_transform_tags("matrix", "vector");
    let fractional_tag_refs: [&str; 12] =
        std::array::from_fn(|index| fractional_tags[index].as_str());
    let fractional_table = block.find("_atom_sites.fract_transf_", &fractional_tag_refs)?;
    if fractional_table.len() > 0 {
        let row = fractional_table
            .row(0)
            .ok_or(MmcifCrystalInfoError::StagingInvariant {
                stage: "nonempty fractional transform table lacks row zero",
            })?;
        let fractional = mmcif_transform_from_row(row)?;
        set_crystal_fractional_transform(crystal, fractional);
    }
    Ok(())
}

fn add_mmcif_info_values(block: &CifBlock, source_state: &mut BioStructureSourceState, tag: &str) {
    // Gemmi✔️✔️: auto add_info = [&](const std::string& tag) {
    // Gemmi✔️✔️:     bool first = true;
    // Gemmi✔️✔️:     for (const std::string& v : block.find_values(tag))
    // Gemmi✔️✔️:         if (!cif::is_null(v)) {
    // Gemmi✔️✔️:             if (first)
    // Gemmi✔️✔️:                 st.info[tag] = cif::as_string(v);
    // Gemmi✔️✔️:             else
    // Gemmi✔️✔️:                 st.info[tag] += "; " + cif::as_string(v);
    // Gemmi✔️✔️:             first = false;
    // Gemmi✔️✔️:         }
    // Gemmi✔️✔️: };
    // Behavior review: use the first matching CIF column in source item order,
    // skip only the raw `.` and `?` null values, decode each retained token,
    // and join source rows with the exact `; ` delimiter. No map entry is
    // created or overwritten when the column is missing or all-null; a
    // decoded empty string is still a retained, present value.
    // Complexity review: one ordered CIF-item lookup, one pass over the
    // selected column, amortized string appends, and at most one BTreeMap
    // insertion; this matches the source's scan and asymptotic map costs.
    let Some(values) = block.find_values(tag) else {
        return;
    };

    let mut combined = None::<String>;
    for value in values.iter() {
        if value.is_null() {
            continue;
        }
        let decoded = value.decoded();
        match &mut combined {
            Some(previous) => {
                previous.push_str("; ");
                previous.push_str(&decoded);
            }
            None => combined = Some(decoded),
        }
    }
    if let Some(combined) = combined {
        source_state.info.insert(tag.to_owned(), combined);
    }
}

fn parse_mmcif_entry_info(block: &CifBlock, source_state: &mut BioStructureSourceState) {
    // Gemmi✔️✔️: void read_entry_info(cif::Block& block, gemmi::Structure& st) {
    // Gemmi✔️✔️:     auto add_info = [&](const std::string& tag) {
    // Gemmi✔️✔️:         bool first = true;
    // Gemmi✔️✔️:         for (const std::string& v : block.find_values(tag))
    // Gemmi✔️✔️:             if (!cif::is_null(v)) {
    // Gemmi✔️✔️:                 if (first)
    // Gemmi✔️✔️:                     st.info[tag] = cif::as_string(v);
    // Gemmi✔️✔️:                 else
    // Gemmi✔️✔️:                     st.info[tag] += "; " + cif::as_string(v);
    // Gemmi✔️✔️:                 first = false;
    // Gemmi✔️✔️:             }
    // Gemmi✔️✔️:     };
    // Gemmi✔️✔️:     add_info("_entry.id");
    // Gemmi✔️✔️:     add_info("_cell.Z_PDB");
    // Gemmi✔️✔️:     add_info("_exptl.method");
    // Gemmi✔️✔️:     add_info("_struct.title");
    // Gemmi✔️✔️:     // in pdbx/mmcif v5 date_original was replaced with a much longer tag
    // Gemmi✔️✔️:     std::string old_date_tag = "_database_PDB_rev.date_original";
    // Gemmi✔️✔️:     std::string new_date_tag = "_pdbx_database_status.recvd_initial_deposition_date";
    // Gemmi✔️✔️:     add_info(old_date_tag);
    // Gemmi✔️✔️:     add_info(new_date_tag);
    // Gemmi✔️✔️:     if (st.info.count(old_date_tag) == 1 && st.info.count(new_date_tag) == 0)
    // Gemmi✔️✔️:         st.info[new_date_tag] = st.info[old_date_tag];
    // Gemmi✔️✔️:     add_info("_struct_keywords.pdbx_keywords");
    // Gemmi✔️✔️:     add_info("_struct_keywords.text");
    // Gemmi✔️✔️: }
    // Behavior review: retain Gemmi's exact eight-key call order and date
    // fallback condition. A new-date key already present in the source state,
    // including a decoded empty value, suppresses copying the old date.
    // Complexity review: exactly eight CIF lookups, each with the source
    // ordered item scan; only selected column values are traversed.
    add_mmcif_info_values(block, source_state, "_entry.id");
    add_mmcif_info_values(block, source_state, "_cell.Z_PDB");
    add_mmcif_info_values(block, source_state, "_exptl.method");
    add_mmcif_info_values(block, source_state, "_struct.title");

    let old_date_tag = "_database_PDB_rev.date_original";
    let new_date_tag = "_pdbx_database_status.recvd_initial_deposition_date";
    add_mmcif_info_values(block, source_state, old_date_tag);
    add_mmcif_info_values(block, source_state, new_date_tag);
    if let Some(old_date) = source_state.info.get(old_date_tag).cloned()
        && !source_state.info.contains_key(new_date_tag)
    {
        source_state.info.insert(new_date_tag.to_owned(), old_date);
    }

    add_mmcif_info_values(block, source_state, "_struct_keywords.pdbx_keywords");
    add_mmcif_info_values(block, source_state, "_struct_keywords.text");
}

fn append_mmcif_audit_authors(block: &CifBlock, metadata: &mut BioMetadata) {
    // Gemmi✔️✔️: void read_audit_author(cif::Block& block, Structure& st) {
    // Gemmi✔️✔️:     for (const std::string& v : block.find_values("_audit_author.name"))
    // Gemmi✔️✔️:         if (!cif::is_null(v))
    // Gemmi✔️✔️:             st.meta.authors.push_back(cif::as_string(v));
    // Gemmi✔️✔️: }
    // Behavior review: use the first source-ordered matching pair/loop
    // column, skip only raw CIF null tokens, decode each remaining value, and
    // append without trimming, deduplication, or replacing prior authors.
    // Complexity review: one ordered item scan plus one pass over the selected
    // values, with amortized Vec pushes; this matches the source's linear scan
    // and append costs without cloning the author collection.
    let Some(values) = block.find_values("_audit_author.name") else {
        return;
    };
    for value in values.iter() {
        if !value.is_null() {
            metadata.authors.push(value.decoded());
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum MmcifTlsInfoError {
    Cif(CifReadError),
    SequenceId(MmcifSequenceIdError),
    NumericIdOutsideSourceDefinedRange { id: String },
    SelectionChainNotRepresentable { chain: String },
}

impl From<CifReadError> for MmcifTlsInfoError {
    fn from(error: CifReadError) -> Self {
        Self::Cif(error)
    }
}

impl From<MmcifSequenceIdError> for MmcifTlsInfoError {
    fn from(error: MmcifSequenceIdError) -> Self {
        Self::SequenceId(error)
    }
}

impl fmt::Display for MmcifTlsInfoError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Cif(error) => fmt::Display::fmt(error, f),
            Self::SequenceId(error) => fmt::Display::fmt(error, f),
            Self::NumericIdOutsideSourceDefinedRange { id } => write!(
                f,
                "TLS group ID is outside the source-defined signed-int range: {id:?}"
            ),
            Self::SelectionChainNotRepresentable { chain } => write!(
                f,
                "TLS selection chain is outside the approved BIO identifier representation: {chain:?}"
            ),
        }
    }
}

impl std::error::Error for MmcifTlsInfoError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Cif(error) => Some(error),
            Self::SequenceId(error) => Some(error),
            Self::NumericIdOutsideSourceDefinedRange { .. }
            | Self::SelectionChainNotRepresentable { .. } => None,
        }
    }
}

fn parse_mmcif_refinement_info(
    block: &CifBlock,
    source_state: &mut BioStructureSourceState,
    metadata: &mut BioMetadata,
) -> Result<(), CifReadError> {
    // Gemmi❗✔️: void read_refinement_info(cif::Block& block, Structure& st) {
    // Gemmi❗✔️:     for (auto row : block.find("_refine.", {"pdbx_refine_id",           // 0
    // Gemmi❗✔️:                                             "?ls_d_res_high",           // 1
    // Gemmi❗✔️:                                             "?ls_d_res_low",            // 2
    // Gemmi❗✔️:                                             "?ls_percent_reflns_obs",   // 3
    // Gemmi❗✔️:                                             "?ls_number_reflns_obs",    // 4
    // Gemmi❗✔️:                                             "?ls_number_reflns_R_work", // 5
    // Gemmi❗✔️:                                             "?ls_number_reflns_R_free", // 6
    // Gemmi❗✔️:                                             "?ls_R_factor_obs",         // 7
    // Gemmi❗✔️:                                             "?ls_R_factor_R_work",      // 8
    // Gemmi❗✔️:                                             "?ls_R_factor_R_free"})) { // 9
    // Gemmi❗✔️:         st.meta.refinement.emplace_back();
    // Gemmi❗✔️:         RefinementInfo& ref = st.meta.refinement.back();
    // Gemmi❗✔️:         ref.id = row.str(0);
    // Gemmi❗✔️:         if (row.has(1)) {
    // Gemmi❗✔️:             ref.resolution_high = cif::as_number(row[1]);
    // Gemmi❗✔️:             if (ref.resolution_high > 0 &&
    // Gemmi❗✔️:                 (st.resolution == 0 || ref.resolution_high < st.resolution))
    // Gemmi❗✔️:                 st.resolution = ref.resolution_high;
    // Gemmi❗✔️:         }
    // Gemmi❗✔️:         copy_double(row, 2, ref.resolution_low);
    // Gemmi❗✔️:         copy_double(row, 3, ref.completeness);
    // Gemmi❗✔️:         copy_int(row, 4, ref.reflection_count);
    // Gemmi❗✔️:         copy_int(row, 5, ref.work_set_count);
    // Gemmi❗✔️:         copy_int(row, 6, ref.rfree_set_count);
    // Gemmi❗✔️:         copy_double(row, 7, ref.r_all);
    // Gemmi❗✔️:         copy_double(row, 8, ref.r_work);
    // Gemmi❗✔️:         copy_double(row, 9, ref.r_free);
    // Gemmi❗✔️:     }
    // Gemmi❗✔️:     if (st.resolution == 0.) {
    // Gemmi❗✔️:       const std::string* em_res = block.find_value("_em_3d_reconstruction.resolution");
    // Gemmi❗✔️:       if (em_res)
    // Gemmi❗✔️:         st.resolution = cif::as_number(*em_res);
    // Gemmi❗✔️:
    // Gemmi❗✔️:     }
    // Gemmi❗✔️: }
    // Gemmi❗✔️: void copy_int(const cif::Table::Row& row, int n, int& dest) {
    // Gemmi❗✔️:   if (row.has2(n))
    // Gemmi❗✔️:     dest = cif::as_int(row[n]);
    // Gemmi❗✔️: }
    // Gemmi❗✔️: void copy_double(const cif::Table::Row& row, int n, double& dest) {
    // Gemmi❗✔️:   if (row.has2(n))
    // Gemmi❗✔️:     dest = cif::as_number(row[n]);
    // Gemmi❗✔️: }
    // Gemmi❗✔️: bool has(size_t n) const { return tab.positions.at(n) >= 0; }
    // Gemmi❗✔️: bool has2(size_t n) const { return has(n) && !cif::is_null(operator[](n)); }
    // Gemmi❗✔️: std::string str(int n) const { return as_string(at(n)); }
    // Gemmi❗✔️: inline double as_number(const std::string& s, double nan=NAN) {
    // Gemmi❗✔️:   const char* start = s.data();
    // Gemmi❗✔️:   const char* end = s.data() + s.size();
    // Gemmi❗✔️:   if (*start == '+')
    // Gemmi❗✔️:     ++start;
    // Gemmi❗✔️:   char f = start[int(*start == '-')] | 0x20;
    // Gemmi❗✔️:   if (f == 'i' || f == 'n')
    // Gemmi❗✔️:     return nan;
    // Gemmi❗✔️:   double d;
    // Gemmi❗✔️:   auto result = fast_float::from_chars(start, end, d);
    // Gemmi❗✔️:   if (result.ec != std::errc())
    // Gemmi❗✔️:     return nan;
    // Gemmi❗✔️:   if (*result.ptr == '(') {
    // Gemmi❗✔️:     const char* p = result.ptr + 1;
    // Gemmi❗✔️:     while (*p >= '0' && *p <= '9')
    // Gemmi❗✔️:       ++p;
    // Gemmi❗✔️:     if (*p == ')')
    // Gemmi❗✔️:       result.ptr = p + 1;
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   return result.ptr == end ? d : nan;
    // Gemmi❗✔️: }
    // Gemmi❗✔️: inline int as_int(const std::string& str) {
    // Gemmi❗✔️:   return string_to_int(str, true);
    // Gemmi❗✔️: }
    // Gemmi❗✔️: inline bool is_null(const std::string& value) {
    // Gemmi❗✔️:   return value.size() == 1 && (value[0] == '?' || value[0] == '.');
    // Gemmi❗✔️: }
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
    // Gemmi❗✔️:   for (; (length == 0 || i < length) && is_digit(p[i]); ++i) {
    // Gemmi❗✔️:     n = n * 10 - (p[i] - '0');
    // Gemmi❗✔️:     has_digits = true;
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   if (checked) {
    // Gemmi❗✔️:     while ((length == 0 || i < length) && is_space(p[i]))
    // Gemmi❗✔️:       ++i;
    // Gemmi❗✔️:     if (!has_digits || p[i] != '\0')
    // Gemmi❗✔️:       throw std::invalid_argument("not an integer: " +
    // Gemmi❗✔️:                                   std::string(p, length ? length : i+1));
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   return mult * n;
    // Gemmi❗✔️: }
    // Behavior review: preserve source row order/IDs, all nine parsed scalar
    // fields, the distinction between column-presence `has` and non-null
    // `has2`, and the positive-minimum then exact-zero EM fallback. Gemmi's
    // reader contains no `_refine_ls_shell` lookup; `bins` and other modeled
    // refinement members therefore retain their source defaults. Float parse
    // failures use the source NaN fallback; checked integer syntax failures
    // are returned as structured CIF errors. Signed C++ int overflow is
    // undefined and is not claimed equivalent outside the `i32` domain.
    // Complexity review: one canonical table selection and a single row pass;
    // each row performs a fixed number of indexed value checks/conversions and
    // one amortized vector append, matching Gemmi's source-shaped scan.
    const REFINE_TAGS: [&str; 10] = [
        "pdbx_refine_id",
        "?ls_d_res_high",
        "?ls_d_res_low",
        "?ls_percent_reflns_obs",
        "?ls_number_reflns_obs",
        "?ls_number_reflns_R_work",
        "?ls_number_reflns_R_free",
        "?ls_R_factor_obs",
        "?ls_R_factor_R_work",
        "?ls_R_factor_R_free",
    ];
    let table = block.find("_refine.", &REFINE_TAGS)?;
    for row in table.iter() {
        metadata.refinement.push(BioRefinementInfo::default());
        let refinement_index = metadata.refinement.len() - 1;
        let refinement = &mut metadata.refinement[refinement_index];
        refinement.id = row.decoded(0).unwrap_or_default();

        if row.has(1) {
            let resolution_high = row
                .get(1)
                .map(|value| cif_as_f64(value.raw(), f64::NAN).unwrap_or(f64::NAN))
                .unwrap_or(f64::NAN);
            refinement.basic.resolution_high = resolution_high;
            if resolution_high > 0.0
                && (source_state.resolution == 0.0 || resolution_high < source_state.resolution)
            {
                source_state.resolution = resolution_high;
            }
        }

        if let Some(value) = row.get(2).filter(|value| !value.is_null()) {
            refinement.basic.resolution_low = cif_as_f64(value.raw(), f64::NAN).unwrap_or(f64::NAN);
        }
        if let Some(value) = row.get(3).filter(|value| !value.is_null()) {
            refinement.basic.completeness = cif_as_f64(value.raw(), f64::NAN).unwrap_or(f64::NAN);
        }
        if let Some(value) = row.get(4).filter(|value| !value.is_null()) {
            refinement.basic.reflection_count = cif_as_i32(value.raw())?;
        }
        if let Some(value) = row.get(5).filter(|value| !value.is_null()) {
            refinement.basic.work_set_count = cif_as_i32(value.raw())?;
        }
        if let Some(value) = row.get(6).filter(|value| !value.is_null()) {
            refinement.basic.rfree_set_count = cif_as_i32(value.raw())?;
        }
        if let Some(value) = row.get(7).filter(|value| !value.is_null()) {
            refinement.basic.r_all = cif_as_f64(value.raw(), f64::NAN).unwrap_or(f64::NAN);
        }
        if let Some(value) = row.get(8).filter(|value| !value.is_null()) {
            refinement.basic.r_work = cif_as_f64(value.raw(), f64::NAN).unwrap_or(f64::NAN);
        }
        if let Some(value) = row.get(9).filter(|value| !value.is_null()) {
            refinement.basic.r_free = cif_as_f64(value.raw(), f64::NAN).unwrap_or(f64::NAN);
        }
    }

    if source_state.resolution == 0.0
        && let Some(em_resolution) = block.find_value("_em_3d_reconstruction.resolution")
    {
        source_state.resolution = cif_as_f64(em_resolution.raw(), f64::NAN).unwrap_or(f64::NAN);
    }

    Ok(())
}

fn parse_mmcif_tls_info(
    block: &CifBlock,
    metadata: &mut BioMetadata,
) -> Result<(), MmcifTlsInfoError> {
    // Gemmi❗✔️: void read_tls_info(cif::Block& block, Structure& st) {
    // Gemmi❗✔️:     for (auto row : block.find("_pdbx_refine_tls.", {
    // Gemmi❗✔️:         "id", "?pdbx_refine_id",
    // Gemmi❗✔️:         "T[1][1]", "T[2][2]", "T[3][3]", "T[1][2]", "T[1][3]", "T[2][3]",
    // Gemmi❗✔️:         "L[1][1]", "L[2][2]", "L[3][3]", "L[1][2]", "L[1][3]", "L[2][3]",
    // Gemmi❗✔️:         "S[1][1]", "S[1][2]", "S[1][3]",
    // Gemmi❗✔️:         "S[2][1]", "S[2][2]", "S[2][3]",
    // Gemmi❗✔️:         "S[3][1]", "S[3][2]", "S[3][3]",
    // Gemmi❗✔️:         "origin_x", "origin_y", "origin_z"})) {
    // Gemmi❗✔️:         if (st.meta.refinement.empty())
    // Gemmi❗✔️:             break;
    // Gemmi❗✔️:         RefinementInfo* ref = nullptr;
    // Gemmi❗✔️:         if (row.has(1))
    // Gemmi❗✔️:             ref = get_by_id(st.meta.refinement, row.str(1));
    // Gemmi❗✔️:         if (!ref)
    // Gemmi❗✔️:             ref = &st.meta.refinement[0];
    // Gemmi❗✔️:         ref->tls_groups.emplace_back();
    // Gemmi❗✔️:         TlsGroup& tls = ref->tls_groups.back();
    // Gemmi❗✔️:         tls.id = row.str(0);
    // Gemmi❗✔️:         tls.num_id = (short) no_sign_atoi(tls.id.c_str());
    // Gemmi❗✔️:         tls.T = get_smat33<double>(row, 2);
    // Gemmi❗✔️:         tls.L = get_smat33<double>(row, 8);
    // Gemmi❗✔️:         for (int i = 0; i < 3; ++i)
    // Gemmi❗✔️:             for (int j = 0; j < 3; ++j)
    // Gemmi❗✔️:                 tls.S[i][j] = cif::as_number(row[14+3*i+j]);
    // Gemmi❗✔️:         tls.origin.x = cif::as_number(row[23]);
    // Gemmi❗✔️:         tls.origin.y = cif::as_number(row[24]);
    // Gemmi❗✔️:         tls.origin.z = cif::as_number(row[25]);
    // Gemmi❗✔️:     }
    // Gemmi❗✔️:     for (auto row : block.find("_pdbx_refine_tls_group.", {
    // Gemmi❗✔️:         "refine_tls_id", "?beg_auth_asym_id", "?beg_auth_seq_id", "?beg_PDB_ins_code",
    // Gemmi❗✔️:         "?end_auth_seq_id", "?end_PDB_ins_code", "?selection_details"})) {
    // Gemmi❗✔️:         for (RefinementInfo& ref : st.meta.refinement)
    // Gemmi❗✔️:             if (TlsGroup* tls = get_by_id(ref.tls_groups, row.str(0))) {
    // Gemmi❗✔️:                 tls->selections.emplace_back();
    // Gemmi❗✔️:                 TlsGroup::Selection& sel = tls->selections.back();
    // Gemmi❗✔️:                 if (row.has(1))
    // Gemmi❗✔️:                     sel.chain = row.str(1);
    // Gemmi❗✔️:                 if (row.has(2))
    // Gemmi❗✔️:                     sel.res_begin = make_seqid(row.str(2), row.ptr_at(3));
    // Gemmi❗✔️:                 if (row.has(4))
    // Gemmi❗✔️:                     sel.res_end = make_seqid(row.str(4), row.ptr_at(5));
    // Gemmi❗✔️:                 if (row.has(6))
    // Gemmi❗✔️:                     sel.details = row.str(6);
    // Gemmi❗✔️:                 break;
    // Gemmi❗✔️:             }
    // Gemmi❗✔️:     }
    // Gemmi❗✔️: }
    // Gemmi❗✔️: template<typename T>
    // Gemmi❗✔️: SMat33<T> get_smat33(cif::Table::Row& row, int n) {
    // Gemmi❗✔️:   return SMat33<T>{(T) cif::as_number(row[n+0]),
    // Gemmi❗✔️:                    (T) cif::as_number(row[n+1]),
    // Gemmi❗✔️:                    (T) cif::as_number(row[n+2]),
    // Gemmi❗✔️:                    (T) cif::as_number(row[n+3]),
    // Gemmi❗✔️:                    (T) cif::as_number(row[n+4]),
    // Gemmi❗✔️:                    (T) cif::as_number(row[n+5])};
    // Gemmi❗✔️: }
    // Gemmi❗✔️: template<class T>
    // Gemmi❗✔️: T* get_by_id(std::vector<T>& vec, const std::string& id) {
    // Gemmi❗✔️:   for (T& item : vec)
    // Gemmi❗✔️:     if (item.id == id)
    // Gemmi❗✔️:       return &item;
    // Gemmi❗✔️:   return nullptr;
    // Gemmi❗✔️: }
    // Behavior review: preserve the no-refinement early break, first exact
    // refinement-ID match with fallback to the first refinement, unsigned
    // decimal-prefix TLS numeric IDs and the fixed-profile i16 narrowing,
    // matrix/component order, and ordered first matching group attachment.
    // Selection fields are assigned only under source column-presence gates;
    // sequence conversion reuses the canonical Gemmi-shaped helper. The
    // approved four-byte ASCII chain representation and source-undefined
    // signed-int overflow remain explicit boundaries, not source parity claims.
    // Complexity review: the two table traversals and nested ordered
    // first-match scans have the same linear/nested-scan shape as the source.
    // The requested column maps are fixed-size; group/selection vectors append
    // once per accepted row, and strings are decoded only where source `str`
    // would create them. No whole-metadata clone or lookup index is added.
    const TLS_TAGS: [&str; 26] = [
        "id",
        "?pdbx_refine_id",
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
        "origin_x",
        "origin_y",
        "origin_z",
    ];
    let tls_table = block.find("_pdbx_refine_tls.", &TLS_TAGS)?;
    for row in tls_table.iter() {
        if metadata.refinement.is_empty() {
            break;
        }

        let refinement_index = if row.has(1) {
            let refinement_id = row.decoded(1).unwrap_or_default();
            metadata
                .refinement
                .iter()
                .position(|refinement| refinement.id == refinement_id)
        } else {
            None
        }
        .unwrap_or(0);

        let refinement = &mut metadata.refinement[refinement_index];
        refinement.tls_groups.push(BioTlsGroup::default());
        let group = refinement
            .tls_groups
            .last_mut()
            .expect("the source row just appended one TLS group");
        group.id = row.decoded(0).unwrap_or_default();
        let (numeric_id, _) =
            super::bio_pdb::gemmi_no_sign_atoi(group.id.as_bytes()).ok_or_else(|| {
                MmcifTlsInfoError::NumericIdOutsideSourceDefinedRange {
                    id: group.id.clone(),
                }
            })?;
        group.num_id = numeric_id as i16;

        let source_number = |column: usize| {
            let value = row.get(column).map(CifValue::raw).unwrap_or_default();
            cif_as_f64(value, f64::NAN).unwrap_or(f64::NAN)
        };
        let symmetric_tensor =
            |first_column| std::array::from_fn(|offset| source_number(first_column + offset));
        group.t = symmetric_tensor(2);
        group.l = symmetric_tensor(8);
        for i in 0..3 {
            for j in 0..3 {
                group.s[i][j] = source_number(14 + 3 * i + j);
            }
        }
        group.origin = [source_number(23), source_number(24), source_number(25)];
    }

    const TLS_SELECTION_TAGS: [&str; 7] = [
        "refine_tls_id",
        "?beg_auth_asym_id",
        "?beg_auth_seq_id",
        "?beg_PDB_ins_code",
        "?end_auth_seq_id",
        "?end_PDB_ins_code",
        "?selection_details",
    ];
    let selection_table = block.find("_pdbx_refine_tls_group.", &TLS_SELECTION_TAGS)?;
    for row in selection_table.iter() {
        let group_id = row.decoded(0).unwrap_or_default();
        let mut match_location = None;
        'refinements: for (refinement_index, refinement) in metadata.refinement.iter().enumerate() {
            for (group_index, group) in refinement.tls_groups.iter().enumerate() {
                if group.id == group_id {
                    match_location = Some((refinement_index, group_index));
                    break 'refinements;
                }
            }
        }
        let Some((refinement_index, group_index)) = match_location else {
            continue;
        };

        let group = &mut metadata.refinement[refinement_index].tls_groups[group_index];
        group.selections.push(BioTlsSelection::default());
        let selection = group
            .selections
            .last_mut()
            .expect("the source row just appended one TLS selection");
        if row.has(1) {
            let chain = row.decoded(1).unwrap_or_default();
            selection.chain = PdbChainId::from_ascii(chain.as_bytes())
                .ok_or_else(|| MmcifTlsInfoError::SelectionChainNotRepresentable { chain })?;
        }
        if row.has(2) {
            selection.res_begin = make_mmcif_seq_id(
                row.decoded(2).unwrap_or_default(),
                row.get(3).map(CifValue::raw),
            )?;
        }
        if row.has(4) {
            selection.res_end = make_mmcif_seq_id(
                row.decoded(4).unwrap_or_default(),
                row.get(5).map(CifValue::raw),
            )?;
        }
        if row.has(6) {
            selection.details = row.decoded(6).unwrap_or_default();
        }
    }

    Ok(())
}

fn find_mmcif_diffraction_mut<'a>(
    metadata: &'a mut BioMetadata,
    diffraction_id: &str,
) -> Option<&'a mut BioDiffractionInfo> {
    // Gemmi❗✔️: DiffractionInfo* find_diffrn(Metadata& meta, const std::string& diffrn_id) {
    // Gemmi❗✔️:   for (CrystalInfo& crystal_info : meta.crystals)
    // Gemmi❗✔️:     for (DiffractionInfo& diffr_info : crystal_info.diffractions)
    // Gemmi❗✔️:       if (diffr_info.id == diffrn_id)
    // Gemmi❗✔️:         return &diffr_info;
    // Gemmi❗✔️:   return nullptr;
    // Gemmi❗✔️: }
    // Behavior review: nested source-order scanning returns the first exact ID
    // match across crystals and their diffraction rows; duplicates are not
    // merged or skipped. Complexity review: the same O(total diffractions)
    // scan is used for each detector/radiation/source row, with no index or
    // collection clone.
    for crystal in &mut metadata.crystals {
        for diffraction in &mut crystal.diffractions {
            if diffraction.id == diffraction_id {
                return Some(diffraction);
            }
        }
    }
    None
}

fn parse_mmcif_experimental_info(
    block: &CifBlock,
    metadata: &mut BioMetadata,
) -> Result<(), CifReadError> {
    // Gemmi❗✔️: void read_experimental_info(cif::Block& block, Structure& st) {
    // Gemmi❗✔️:     for (auto row : block.find("_exptl.", {"method", "?crystals_number"})) {
    // Gemmi❗✔️:         st.meta.experiments.emplace_back();
    // Gemmi❗✔️:         st.meta.experiments.back().method = row.str(0);
    // Gemmi❗✔️:         copy_int(row, 1, st.meta.experiments.back().number_of_crystals);
    // Gemmi❗✔️:     }
    // Gemmi❗✔️:     for (auto row : block.find("_exptl_crystal.", {"id", "?description"})) {
    // Gemmi❗✔️:         st.meta.crystals.emplace_back();
    // Gemmi❗✔️:         st.meta.crystals.back().id = row.str(0);
    // Gemmi❗✔️:         copy_string(row, 1, st.meta.crystals.back().description);
    // Gemmi❗✔️:     }
    // Gemmi❗✔️:     for (auto row : block.find("_diffrn.",
    // Gemmi❗✔️:                                {"id", "crystal_id", "?ambient_temp"})) {
    // Gemmi❗✔️:         std::string id = row.str(1);
    // Gemmi❗✔️:         auto cryst = std::find_if(st.meta.crystals.begin(), st.meta.crystals.end(),
    // Gemmi❗✔️:                                   [&](const CrystalInfo& c) { return c.id == id; });
    // Gemmi❗✔️:         if (cryst != st.meta.crystals.end()) {
    // Gemmi❗✔️:             cryst->diffractions.emplace_back();
    // Gemmi❗✔️:             cryst->diffractions.back().id = row.str(0);
    // Gemmi❗✔️:             copy_double(row, 2, cryst->diffractions.back().temperature);
    // Gemmi❗✔️:         }
    // Gemmi❗✔️:     }
    // Gemmi❗✔️:     for (auto row : block.find("_diffrn_detector.", {"diffrn_id",
    // Gemmi❗✔️:                                                      "?pdbx_collection_date",
    // Gemmi❗✔️:                                                      "?detector",
    // Gemmi❗✔️:                                                      "?type",
    // Gemmi❗✔️:                                                      "?details"}))
    // Gemmi❗✔️:         if (DiffractionInfo* di = find_diffrn(st.meta, row.str(0))) {
    // Gemmi❗✔️:             copy_string(row, 1, di->collection_date);
    // Gemmi❗✔️:             copy_string(row, 2, di->detector);
    // Gemmi❗✔️:             copy_string(row, 3, di->detector_make);
    // Gemmi❗✔️:             copy_string(row, 4, di->optics);
    // Gemmi❗✔️:         }
    // Gemmi❗✔️:     for (auto row : block.find("_diffrn_radiation.",
    // Gemmi❗✔️:                                {"diffrn_id",
    // Gemmi❗✔️:                                 "?pdbx_scattering_type",
    // Gemmi❗✔️:                                 "?pdbx_monochromatic_or_laue_m_l",
    // Gemmi❗✔️:                                 "?monochromator"}))
    // Gemmi❗✔️:         if (DiffractionInfo* di = find_diffrn(st.meta, row.str(0))) {
    // Gemmi❗✔️:             copy_string(row, 1, di->scattering_type);
    // Gemmi❗✔️:             if (row.has2(2))
    // Gemmi❗✔️:                 di->mono_or_laue = row.str(2)[0];
    // Gemmi❗✔️:             copy_string(row, 3, di->monochromator);
    // Gemmi❗✔️:         }
    // Gemmi❗✔️:     for (auto row : block.find("_diffrn_source.", {"diffrn_id",
    // Gemmi❗✔️:                                                    "?source",
    // Gemmi❗✔️:                                                    "?type",
    // Gemmi❗✔️:                                                    "?pdbx_synchrotron_site",
    // Gemmi❗✔️:                                                    "?pdbx_synchrotron_beamline",
    // Gemmi❗✔️:                                                    "?pdbx_wavelength_list"}))
    // Gemmi❗✔️:         if (DiffractionInfo* di = find_diffrn(st.meta, row.str(0))) {
    // Gemmi❗✔️:             copy_string(row, 1, di->source);
    // Gemmi❗✔️:             copy_string(row, 2, di->source_type);
    // Gemmi❗✔️:             copy_string(row, 3, di->synchrotron);
    // Gemmi❗✔️:             copy_string(row, 4, di->beamline);
    // Gemmi❗✔️:             copy_string(row, 5, di->wavelengths);
    // Gemmi❗✔️:         }
    // Gemmi❗✔️: }
    // Gemmi❗✔️: void copy_int(const cif::Table::Row& row, int n, int& dest) {
    // Gemmi❗✔️:   if (row.has2(n))
    // Gemmi❗✔️:     dest = cif::as_int(row[n]);
    // Gemmi❗✔️: }
    // Gemmi❗✔️: void copy_double(const cif::Table::Row& row, int n, double& dest) {
    // Gemmi❗✔️:   if (row.has2(n))
    // Gemmi❗✔️:     dest = cif::as_number(row[n]);
    // Gemmi❗✔️: }
    // Gemmi❗✔️: void copy_string(const cif::Table::Row& row, int n, std::string& dest) {
    // Gemmi❗✔️:   if (row.has2(n))
    // Gemmi❗✔️:     dest = cif::as_string(row[n]);
    // Gemmi❗✔️: }
    // Gemmi❗✔️: bool has(size_t n) const { return tab.positions.at(n) >= 0; }
    // Gemmi❗✔️: bool has2(size_t n) const { return has(n) && !cif::is_null(operator[](n)); }
    // Gemmi❗✔️: std::string str(int n) const { return as_string(at(n)); }
    // Behavior review: preserve all six table passes and their order, appending
    // experiment/crystal/diffraction records at the same points as Gemmi.
    // Required identifiers and method strings use decoded source `str`
    // behavior. Optional copies require both column presence and a non-null
    // raw value; quoted empty strings still assign. Crystal and diffraction
    // joins are first exact source-order matches; unmatched association rows
    // do nothing. The radiation code retains the first decoded byte, with an
    // empty decoded string mapping to C++'s NUL character at `string[0]`.
    // Checked integer conversion uses the existing CIF owner; signed overflow
    // remains the separately documented source-undefined boundary. This
    // function does not import `_exptl_crystal_grow` or `_reflns` data, whose
    // source readers are distinct.
    // Complexity review: table/row traversal remains linear per category;
    // crystal association scans crystals once per `_diffrn` row, and metadata
    // association scans all diffractions once per detector/radiation/source
    // row, matching the source's first-match scans. No lookup map or cloned
    // metadata is introduced; decoded strings allocate only when assigned.
    const EXPERIMENT_TAGS: [&str; 2] = ["method", "?crystals_number"];
    let experiments = block.find("_exptl.", &EXPERIMENT_TAGS)?;
    for row in experiments.iter() {
        metadata.experiments.push(BioExperimentInfo::default());
        let experiment = metadata
            .experiments
            .last_mut()
            .expect("the source row just appended one experiment");
        experiment.method = row.decoded(0).unwrap_or_default();
        if let Some(value) = row.get(1).filter(|value| !value.is_null()) {
            experiment.number_of_crystals = cif_as_i32(value.raw())?;
        }
    }

    const CRYSTAL_TAGS: [&str; 2] = ["id", "?description"];
    let crystals = block.find("_exptl_crystal.", &CRYSTAL_TAGS)?;
    for row in crystals.iter() {
        metadata
            .crystals
            .push(BioExperimentalCrystalInfo::default());
        let crystal = metadata
            .crystals
            .last_mut()
            .expect("the source row just appended one crystal");
        crystal.id = row.decoded(0).unwrap_or_default();
        if let Some(value) = row.get(1).filter(|value| !value.is_null()) {
            crystal.description = value.decoded();
        }
    }

    const DIFFRACTION_TAGS: [&str; 3] = ["id", "crystal_id", "?ambient_temp"];
    let diffractions = block.find("_diffrn.", &DIFFRACTION_TAGS)?;
    for row in diffractions.iter() {
        let crystal_id = row.decoded(1).unwrap_or_default();
        let Some(crystal_index) = metadata
            .crystals
            .iter()
            .position(|crystal| crystal.id == crystal_id)
        else {
            continue;
        };
        let crystal = &mut metadata.crystals[crystal_index];
        crystal.diffractions.push(BioDiffractionInfo::default());
        let diffraction = crystal
            .diffractions
            .last_mut()
            .expect("the source row just appended one diffraction");
        diffraction.id = row.decoded(0).unwrap_or_default();
        if let Some(value) = row.get(2).filter(|value| !value.is_null()) {
            diffraction.temperature = cif_as_f64(value.raw(), f64::NAN)?;
        }
    }

    const DETECTOR_TAGS: [&str; 5] = [
        "diffrn_id",
        "?pdbx_collection_date",
        "?detector",
        "?type",
        "?details",
    ];
    let detectors = block.find("_diffrn_detector.", &DETECTOR_TAGS)?;
    for row in detectors.iter() {
        let diffraction_id = row.decoded(0).unwrap_or_default();
        if let Some(diffraction) = find_mmcif_diffraction_mut(metadata, &diffraction_id) {
            if let Some(value) = row.get(1).filter(|value| !value.is_null()) {
                diffraction.collection_date = value.decoded();
            }
            if let Some(value) = row.get(2).filter(|value| !value.is_null()) {
                diffraction.detector = value.decoded();
            }
            if let Some(value) = row.get(3).filter(|value| !value.is_null()) {
                diffraction.detector_make = value.decoded();
            }
            if let Some(value) = row.get(4).filter(|value| !value.is_null()) {
                diffraction.optics = value.decoded();
            }
        }
    }

    const RADIATION_TAGS: [&str; 4] = [
        "diffrn_id",
        "?pdbx_scattering_type",
        "?pdbx_monochromatic_or_laue_m_l",
        "?monochromator",
    ];
    let radiation = block.find("_diffrn_radiation.", &RADIATION_TAGS)?;
    for row in radiation.iter() {
        let diffraction_id = row.decoded(0).unwrap_or_default();
        if let Some(diffraction) = find_mmcif_diffraction_mut(metadata, &diffraction_id) {
            if let Some(value) = row.get(1).filter(|value| !value.is_null()) {
                diffraction.scattering_type = value.decoded();
            }
            if let Some(value) = row.get(2).filter(|value| !value.is_null()) {
                let code = value.decoded();
                diffraction.mono_or_laue = code.as_bytes().first().copied().unwrap_or_default();
            }
            if let Some(value) = row.get(3).filter(|value| !value.is_null()) {
                diffraction.monochromator = value.decoded();
            }
        }
    }

    const SOURCE_TAGS: [&str; 6] = [
        "diffrn_id",
        "?source",
        "?type",
        "?pdbx_synchrotron_site",
        "?pdbx_synchrotron_beamline",
        "?pdbx_wavelength_list",
    ];
    let sources = block.find("_diffrn_source.", &SOURCE_TAGS)?;
    for row in sources.iter() {
        let diffraction_id = row.decoded(0).unwrap_or_default();
        if let Some(diffraction) = find_mmcif_diffraction_mut(metadata, &diffraction_id) {
            if let Some(value) = row.get(1).filter(|value| !value.is_null()) {
                diffraction.source = value.decoded();
            }
            if let Some(value) = row.get(2).filter(|value| !value.is_null()) {
                diffraction.source_type = value.decoded();
            }
            if let Some(value) = row.get(3).filter(|value| !value.is_null()) {
                diffraction.synchrotron = value.decoded();
            }
            if let Some(value) = row.get(4).filter(|value| !value.is_null()) {
                diffraction.beamline = value.decoded();
            }
            if let Some(value) = row.get(5).filter(|value| !value.is_null()) {
                diffraction.wavelengths = value.decoded();
            }
        }
    }

    Ok(())
}

fn parse_mmcif_reflns_info(
    block: &CifBlock,
    metadata: &mut BioMetadata,
) -> Result<(), CifReadError> {
    // Gemmi❗✔️: void read_reflns_info(cif::Block& block, Structure& st) {
    // Gemmi❗✔️:     size_t n = 0;
    // Gemmi❗✔️:     for (auto row : block.find("_reflns.", {"pdbx_diffrn_id",        // 0
    // Gemmi❗✔️:                                             "?number_obs",           // 1
    // Gemmi❗✔️:                                             "?d_resolution_high",    // 2
    // Gemmi❗✔️:                                             "?d_resolution_low",     // 3
    // Gemmi❗✔️:                                             "?percent_possible_obs", // 4
    // Gemmi❗✔️:                                             "?pdbx_redundancy",      // 5
    // Gemmi❗✔️:                                             "?pdbx_Rmerge_I_obs",    // 6
    // Gemmi❗✔️:                                             "?pdbx_Rsym_value",      // 7
    // Gemmi❗✔️:                                             "?pdbx_netI_over_sigmaI"})) {
    // Gemmi❗✔️:         if (n >= st.meta.experiments.size())
    // Gemmi❗✔️:             break;
    // Gemmi❗✔️:         ExperimentInfo& exper = st.meta.experiments[n++];
    // Gemmi❗✔️:         split_str_into(row.str(0), ',', exper.diffraction_ids);
    // Gemmi❗✔️:         copy_int(row, 1, exper.unique_reflections);
    // Gemmi❗✔️:         copy_double(row, 2, exper.reflections.resolution_high);
    // Gemmi❗✔️:         copy_double(row, 3, exper.reflections.resolution_low);
    // Gemmi❗✔️:         copy_double(row, 4, exper.reflections.completeness);
    // Gemmi❗✔️:         copy_double(row, 5, exper.reflections.redundancy);
    // Gemmi❗✔️:         copy_double(row, 6, exper.reflections.r_merge);
    // Gemmi❗✔️:         copy_double(row, 7, exper.reflections.r_sym);
    // Gemmi❗✔️:         copy_double(row, 8, exper.reflections.mean_I_over_sigma);
    // Gemmi❗✔️:     }
    // Gemmi❗✔️: }
    // Gemmi❗✔️: void copy_int(const cif::Table::Row& row, int n, int& dest) {
    // Gemmi❗✔️:   if (row.has2(n))
    // Gemmi❗✔️:     dest = cif::as_int(row[n]);
    // Gemmi❗✔️: }
    // Gemmi❗✔️: void copy_double(const cif::Table::Row& row, int n, double& dest) {
    // Gemmi❗✔️:   if (row.has2(n))
    // Gemmi❗✔️:     dest = cif::as_number(row[n]);
    // Gemmi❗✔️: }
    // Gemmi❗✔️: bool has2(size_t n) const { return has(n) && !cif::is_null(operator[](n)); }
    // Gemmi❗✔️: std::string str(int n) const { return as_string(at(n)); }
    // Gemmi❗✔️: inline size_t length(char) { return 1; }
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
    // Behavior review: pair each `_reflns` row to the next experiment by
    // source order, stopping before any excess row is decoded or converted.
    // The comma split appends all fields, including empty leading, adjacent,
    // trailing, or sole-empty tokens. Optional numeric columns update only
    // when present and non-null; defaults and existing values otherwise stay.
    // `_reflns_shell` and `_reflns.B_iso_Wilson_estimate` are not read here.
    // Source signed overflow remains undefined; floating numeric limits remain
    // those of the existing canonical CIF conversion owner.
    // Complexity review: one row pass, one linear split per consumed row, and
    // a fixed number of scalar copies; no experiment search, sorting, or clone.
    const REFLECTION_TAGS: [&str; 9] = [
        "pdbx_diffrn_id",
        "?number_obs",
        "?d_resolution_high",
        "?d_resolution_low",
        "?percent_possible_obs",
        "?pdbx_redundancy",
        "?pdbx_Rmerge_I_obs",
        "?pdbx_Rsym_value",
        "?pdbx_netI_over_sigmaI",
    ];
    let rows = block.find("_reflns.", &REFLECTION_TAGS)?;
    for (index, row) in rows.iter().enumerate() {
        if index >= metadata.experiments.len() {
            break;
        }
        let experiment = &mut metadata.experiments[index];
        let diffraction_ids = row.decoded(0).unwrap_or_default();
        experiment
            .diffraction_ids
            .extend(diffraction_ids.split(',').map(str::to_owned));

        if let Some(value) = row.get(1).filter(|value| !value.is_null()) {
            experiment.unique_reflections = cif_as_i32(value.raw())?;
        }
        if let Some(value) = row.get(2).filter(|value| !value.is_null()) {
            experiment.reflections.resolution_high = cif_as_f64(value.raw(), f64::NAN)?;
        }
        if let Some(value) = row.get(3).filter(|value| !value.is_null()) {
            experiment.reflections.resolution_low = cif_as_f64(value.raw(), f64::NAN)?;
        }
        if let Some(value) = row.get(4).filter(|value| !value.is_null()) {
            experiment.reflections.completeness = cif_as_f64(value.raw(), f64::NAN)?;
        }
        if let Some(value) = row.get(5).filter(|value| !value.is_null()) {
            experiment.reflections.redundancy = cif_as_f64(value.raw(), f64::NAN)?;
        }
        if let Some(value) = row.get(6).filter(|value| !value.is_null()) {
            experiment.reflections.r_merge = cif_as_f64(value.raw(), f64::NAN)?;
        }
        if let Some(value) = row.get(7).filter(|value| !value.is_null()) {
            experiment.reflections.r_sym = cif_as_f64(value.raw(), f64::NAN)?;
        }
        if let Some(value) = row.get(8).filter(|value| !value.is_null()) {
            experiment.reflections.mean_i_over_sigma = cif_as_f64(value.raw(), f64::NAN)?;
        }
    }

    Ok(())
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum MmcifHelixError {
    Cif(CifReadError),
    ChainNameNotRepresentable {
        endpoint: &'static str,
        value: String,
    },
    ResidueNameNotRepresentable {
        endpoint: &'static str,
        value: String,
    },
    ResidueAddressNotRepresentable {
        endpoint: &'static str,
    },
    SequenceId {
        endpoint: &'static str,
        error: MmcifSequenceIdError,
    },
    HelixClass(CifReadError),
    HelixLength(CifReadError),
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum MmcifSheetError {
    Cif(CifReadError),
    ChainNameNotRepresentable {
        address: &'static str,
        value: String,
    },
    ResidueNameNotRepresentable {
        address: &'static str,
        value: String,
    },
    ResidueAddressNotRepresentable {
        address: &'static str,
    },
    SequenceId {
        address: &'static str,
        error: MmcifSequenceIdError,
    },
}

impl From<CifReadError> for MmcifSheetError {
    fn from(error: CifReadError) -> Self {
        Self::Cif(error)
    }
}

impl fmt::Display for MmcifSheetError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Cif(error) => fmt::Display::fmt(error, f),
            Self::ChainNameNotRepresentable { address, value } => write!(
                f,
                "{address} sheet chain name exceeds the approved BIO representation: {value:?}"
            ),
            Self::ResidueNameNotRepresentable { address, value } => write!(
                f,
                "{address} sheet residue name exceeds the approved BIO representation: {value:?}"
            ),
            Self::ResidueAddressNotRepresentable { address } => {
                write!(f, "{address} sheet residue address is not representable")
            }
            Self::SequenceId { address, error } => {
                write!(f, "invalid {address} sheet sequence identifier: {error}")
            }
        }
    }
}

impl std::error::Error for MmcifSheetError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Cif(error) => Some(error),
            Self::SequenceId { error, .. } => Some(error),
            Self::ChainNameNotRepresentable { .. }
            | Self::ResidueNameNotRepresentable { .. }
            | Self::ResidueAddressNotRepresentable { .. } => None,
        }
    }
}

impl From<CifReadError> for MmcifHelixError {
    fn from(error: CifReadError) -> Self {
        Self::Cif(error)
    }
}

impl fmt::Display for MmcifHelixError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Cif(error) => fmt::Display::fmt(error, f),
            Self::ChainNameNotRepresentable { endpoint, value } => write!(
                f,
                "{endpoint} helix chain name exceeds the approved BIO representation: {value:?}"
            ),
            Self::ResidueNameNotRepresentable { endpoint, value } => write!(
                f,
                "{endpoint} helix residue name exceeds the approved BIO representation: {value:?}"
            ),
            Self::ResidueAddressNotRepresentable { endpoint } => {
                write!(f, "{endpoint} helix residue address is not representable")
            }
            Self::SequenceId { endpoint, error } => {
                write!(f, "invalid {endpoint} helix sequence identifier: {error}")
            }
            Self::HelixClass(error) => write!(f, "invalid mmCIF helix class: {error}"),
            Self::HelixLength(error) => write!(f, "invalid mmCIF helix length: {error}"),
        }
    }
}

impl std::error::Error for MmcifHelixError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Cif(error) => Some(error),
            Self::SequenceId { error, .. } => Some(error),
            Self::HelixClass(error) | Self::HelixLength(error) => Some(error),
            Self::ChainNameNotRepresentable { .. }
            | Self::ResidueNameNotRepresentable { .. }
            | Self::ResidueAddressNotRepresentable { .. } => None,
        }
    }
}

fn make_mmcif_helix_endpoint(
    row: &CifRow<'_>,
    chain_column: usize,
    name_column: usize,
    sequence_column: usize,
    insertion_column: usize,
    endpoint: &'static str,
) -> Result<AtomAddress, MmcifHelixError> {
    // Gemmi❗✔️: inline ResidueId make_resid(const std::string& name,
    // Gemmi❗✔️:                             const std::string& seqid,
    // Gemmi❗✔️:                             const std::string* icode) {
    // Gemmi❗✔️:   return ResidueId{make_seqid(seqid, icode), {}, name};
    // Gemmi❗✔️: }
    // Behavior review: the canonical IO sequence helper implements the full
    // `make_seqid` closure; the endpoint name and chain are retained exactly
    // within their already-approved four-byte ASCII BIO representation.
    // A present nullable insertion column passes its raw value, while an
    // absent column passes no pointer, matching `row.ptr_at`.
    // Complexity review: a fixed number of row lookups and at most one linear
    // numeric conversion; owned decoded strings mirror source `row.str`.
    let chain_text = row.decoded(chain_column).unwrap_or_default();
    let chain = PdbChainId::from_ascii(chain_text.as_bytes()).ok_or_else(|| {
        MmcifHelixError::ChainNameNotRepresentable {
            endpoint,
            value: chain_text,
        }
    })?;
    let residue_text = row.decoded(name_column).unwrap_or_default();
    let residue_name = ResidueName::from_ascii(residue_text.as_bytes()).ok_or_else(|| {
        MmcifHelixError::ResidueNameNotRepresentable {
            endpoint,
            value: residue_text,
        }
    })?;
    let sequence_text = row.decoded(sequence_column).unwrap_or_default();
    let insertion_text = row
        .has(insertion_column)
        .then(|| row.get(insertion_column).map(CifValue::raw).unwrap_or(""));
    let sequence_id = make_mmcif_seq_id(sequence_text, insertion_text)
        .map_err(|error| MmcifHelixError::SequenceId { endpoint, error })?;
    let residue = make_mmcif_residue_address(residue_name, sequence_id)
        .ok_or(MmcifHelixError::ResidueAddressNotRepresentable { endpoint })?;

    Ok(AtomAddress::new(chain, residue, String::new(), None))
}

fn parse_mmcif_helices(block: &CifBlock) -> Result<Vec<BioHelix>, MmcifHelixError> {
    // Gemmi❗✔️: std::vector<Helix> read_helices(cif::Block& block) {
    // Gemmi❗✔️:   std::vector<Helix> helices;
    // Gemmi❗✔️:   for (const auto row : block.find("_struct_conf.", {
    // Gemmi❗✔️:         "conf_type_id",                                          // 0
    // Gemmi❗✔️:         "beg_auth_asym_id", "beg_label_comp_id",                 // 1-2
    // Gemmi❗✔️:         "beg_auth_seq_id", "?pdbx_beg_PDB_ins_code",             // 3-4
    // Gemmi❗✔️:         "end_auth_asym_id", "end_label_comp_id",                 // 5-6
    // Gemmi❗✔️:         "end_auth_seq_id", "?pdbx_end_PDB_ins_code",             // 7-8
    // Gemmi❗✔️:         "?pdbx_PDB_helix_class", "?pdbx_PDB_helix_length"})) {  // 9-10
    // Gemmi❗✔️:     if (alpha_up(row.str(0)[0]) != 'H')
    // Gemmi❗✔️:       continue;
    // Gemmi❗✔️:     Helix h;
    // Gemmi❗✔️:     h.start.chain_name = row.str(1);
    // Gemmi❗✔️:     h.start.res_id = make_resid(row.str(2), row.str(3), row.ptr_at(4));
    // Gemmi❗✔️:     h.end.chain_name = row.str(5);
    // Gemmi❗✔️:     h.end.res_id = make_resid(row.str(6), row.str(7), row.ptr_at(8));
    // Gemmi❗✔️:     if (row.has(9))
    // Gemmi❗✔️:       h.set_helix_class_as_int(cif::as_int(row[9], -1));
    // Gemmi❗✔️:     if (row.has(10))
    // Gemmi❗✔️:       h.length = cif::as_int(row[10], -1);
    // Gemmi❗✔️:     helices.push_back(h);
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   return helices;
    // Gemmi❗✔️: }
    // Gemmi❗✔️: inline char alpha_up(char c) { return c & ~0x20; }
    // Gemmi❗✔️: inline int as_int(const std::string& str, int null) {
    // Gemmi❗✔️:   return is_null(str) ? null : as_int(str);
    // Gemmi❗✔️: }
    // Behavior review: source row order and first-byte filtering are retained.
    // Missing/null insertion and optional integer columns preserve the source
    // pointer/presence distinction; class values outside 1..=10 retain the
    // canonical Unknown class. Unsupported source string widths fail with a
    // structured boundary error; signed source integer overflow is undefined.
    // Complexity review: one row traversal, fixed per-row endpoint accesses,
    // and source-shaped output allocations; no search, sort, or row clone.
    const HELIX_TAGS: [&str; 11] = [
        "conf_type_id",
        "beg_auth_asym_id",
        "beg_label_comp_id",
        "beg_auth_seq_id",
        "?pdbx_beg_PDB_ins_code",
        "end_auth_asym_id",
        "end_label_comp_id",
        "end_auth_seq_id",
        "?pdbx_end_PDB_ins_code",
        "?pdbx_PDB_helix_class",
        "?pdbx_PDB_helix_length",
    ];
    let rows = block.find("_struct_conf.", &HELIX_TAGS)?;
    let mut helices = Vec::new();
    for row in rows.iter() {
        let conf_type = row.decoded(0).unwrap_or_default();
        if conf_type.as_bytes().first().copied().unwrap_or_default() & !0x20 != b'H' {
            continue;
        }
        let start = make_mmcif_helix_endpoint(&row, 1, 2, 3, 4, "start")?;
        let end = make_mmcif_helix_endpoint(&row, 5, 6, 7, 8, "end")?;
        let mut helix = BioHelix {
            start,
            end,
            ..BioHelix::default()
        };
        if row.has(9) {
            let raw = row.get(9).map(CifValue::raw).unwrap_or("");
            let value = if cif_is_null(raw) {
                -1
            } else {
                cif_as_i32(raw).map_err(MmcifHelixError::HelixClass)?
            };
            helix.set_helix_class_as_int(value);
        }
        if row.has(10) {
            let raw = row.get(10).map(CifValue::raw).unwrap_or("");
            helix.length = if cif_is_null(raw) {
                -1
            } else {
                cif_as_i32(raw).map_err(MmcifHelixError::HelixLength)?
            };
        }
        helices.push(helix);
    }
    Ok(helices)
}

fn make_mmcif_sheet_address(
    row: &CifRow<'_>,
    chain_column: usize,
    residue_name_column: usize,
    sequence_column: usize,
    insertion_column: usize,
    atom_name_column: Option<usize>,
    address: &'static str,
) -> Result<AtomAddress, MmcifSheetError> {
    let chain_text = row.decoded(chain_column).unwrap_or_default();
    let chain = PdbChainId::from_ascii(chain_text.as_bytes()).ok_or_else(|| {
        MmcifSheetError::ChainNameNotRepresentable {
            address,
            value: chain_text,
        }
    })?;
    let residue_text = row.decoded(residue_name_column).unwrap_or_default();
    let residue_name = ResidueName::from_ascii(residue_text.as_bytes()).ok_or_else(|| {
        MmcifSheetError::ResidueNameNotRepresentable {
            address,
            value: residue_text,
        }
    })?;
    let sequence_text = row.decoded(sequence_column).unwrap_or_default();
    let insertion_text = row
        .has(insertion_column)
        .then(|| row.get(insertion_column).map(CifValue::raw).unwrap_or(""));
    let sequence_id = make_mmcif_seq_id(sequence_text, insertion_text)
        .map_err(|error| MmcifSheetError::SequenceId { address, error })?;
    let residue = make_mmcif_residue_address(residue_name, sequence_id)
        .ok_or(MmcifSheetError::ResidueAddressNotRepresentable { address })?;
    let atom_name = atom_name_column
        .and_then(|column| row.decoded(column))
        .unwrap_or_default();

    Ok(AtomAddress::new(chain, residue, atom_name, None))
}

fn parse_mmcif_sheets(block: &CifBlock) -> Result<Vec<BioSheet>, MmcifSheetError> {
    // Gemmi❗✔️: std::vector<Sheet> read_sheets(cif::Block& block) {
    // Gemmi❗✔️:   std::vector<Sheet> sheets;
    // Gemmi❗✔️:   for (const std::string& sheet_id : block.find_values("_struct_sheet.id"))
    // Gemmi❗✔️:     sheets.emplace_back(sheet_id);
    // Gemmi❗✔️:   for (const auto row : block.find("_struct_sheet_range.", {
    // Gemmi❗✔️:         "sheet_id", "id",                                        // 0-1
    // Gemmi❗✔️:         "beg_auth_asym_id", "beg_label_comp_id",                 // 2-3
    // Gemmi❗✔️:         "beg_auth_seq_id", "?pdbx_beg_PDB_ins_code",             // 4-5
    // Gemmi❗✔️:         "end_auth_asym_id", "end_label_comp_id",                 // 6-7
    // Gemmi❗✔️:         "end_auth_seq_id", "?pdbx_end_PDB_ins_code"})) {        // 8-9
    // Gemmi❗✔️:     Sheet& sheet = impl::find_or_add(sheets, row.str(0));
    // Gemmi❗✔️:     sheet.strands.emplace_back();
    // Gemmi❗✔️:     Sheet::Strand& strand = sheet.strands.back();
    // Gemmi❗✔️:     strand.name = row.str(1);
    // Gemmi❗✔️:     strand.start.chain_name = row.str(2);
    // Gemmi❗✔️:     strand.start.res_id = make_resid(row.str(3), row.str(4), row.ptr_at(5));
    // Gemmi❗✔️:     strand.end.chain_name = row.str(6);
    // Gemmi❗✔️:     strand.end.res_id = make_resid(row.str(7), row.str(8), row.ptr_at(9));
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   // below we assume that range_id_1 is the strand preceding range_id_2
    // Gemmi❗✔️:   for (const auto row : block.find("_struct_sheet_order.", {
    // Gemmi❗✔️:         "sheet_id", "range_id_2", "sense"}))
    // Gemmi❗✔️:     if (Sheet* sheet = impl::find_or_null(sheets, row.str(0)))
    // Gemmi❗✔️:       if (Sheet::Strand* ss = impl::find_or_null(sheet->strands, row.str(1)))
    // Gemmi❗✔️:         switch (alpha_up(row.str(2)[0])) {
    // Gemmi❗✔️:           case 'P': ss->sense = 1; break; // parallel
    // Gemmi❗✔️:           case 'A': ss->sense = -1; break; // anti-parallel
    // Gemmi❗✔️:         }
    // Gemmi❗✔️:
    // Gemmi❗✔️:   for (const auto row : block.find("_pdbx_struct_sheet_hbond.", {
    // Gemmi❗✔️:         "sheet_id", "range_id_2",
    // Gemmi❗✔️:         "range_1_auth_asym_id", "range_1_label_comp_id",
    // Gemmi❗✔️:         "range_1_auth_seq_id", "?range_1_PDB_ins_code",
    // Gemmi❗✔️:         "range_1_label_atom_id",
    // Gemmi❗✔️:         "range_2_auth_asym_id", "range_2_label_comp_id",
    // Gemmi❗✔️:         "range_2_auth_seq_id", "?range_2_PDB_ins_code",
    // Gemmi❗✔️:         "range_2_label_atom_id"}))
    // Gemmi❗✔️:     if (Sheet* sheet = impl::find_or_null(sheets, row.str(0)))
    // Gemmi❗✔️:       if (Sheet::Strand* ss = impl::find_or_null(sheet->strands, row.str(1))) {
    // Gemmi❗✔️:         ss->hbond_atom1.chain_name = row.str(2);
    // Gemmi❗✔️:         ss->hbond_atom1.res_id = make_resid(row.str(3), row.str(4),
    // Gemmi❗✔️:                                             row.ptr_at(5));
    // Gemmi❗✔️:         ss->hbond_atom1.atom_name = row.str(6);
    // Gemmi❗✔️:         ss->hbond_atom2.chain_name = row.str(7);
    // Gemmi❗✔️:         ss->hbond_atom2.res_id = make_resid(row.str(8), row.str(9),
    // Gemmi❗✔️:                                             row.ptr_at(10));
    // Gemmi❗✔️:         ss->hbond_atom2.atom_name = row.str(11);
    // Gemmi❗✔️:       }
    // Gemmi❗✔️:
    // Gemmi❗✔️:   return sheets;
    // Gemmi❗✔️: }
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
    // Gemmi❗✔️: inline Column Block::find_values(const std::string& tag) {
    // Gemmi❗✔️:   std::string lctag = gemmi::to_lower(tag);
    // Gemmi❗✔️:   for (Item& i : items)
    // Gemmi❗✔️:     if (i.type == ItemType::Loop) {
    // Gemmi❗✔️:       int pos = i.loop.find_tag_lc(lctag);
    // Gemmi❗✔️:       if (pos != -1)
    // Gemmi❗✔️:         return Column{&i, static_cast<size_t>(pos)};
    // Gemmi❗✔️:     } else if (i.type == ItemType::Pair) {
    // Gemmi❗✔️:       if (gemmi::iequal(i.pair[0], lctag))
    // Gemmi❗✔️:         return Column{&i, 0};
    // Gemmi❗✔️:     }
    // Gemmi❗✔️:   return Column{nullptr, 0};
    // Gemmi❗✔️: }
    // Gemmi❗✔️: inline char alpha_up(char c) { return c & ~0x20; }
    // Behavior review: this preserves source row order, duplicate IDs, first
    // exact-name matching, append-on-range creation, foreign-key ignoring,
    // P/A sense updates and both full source-address projections. Strand
    // `sense` starts at zero because the source's empty-argument
    // `vector::emplace_back()` value-initializes the aggregate field; this was
    // verified against the pinned header/compiler. Errors at the approved
    // fixed-width BIO boundary remain typed and do not truncate names.
    // Complexity review: each source phase traverses its own rows once and
    // performs the same linear first-match vector searches as Gemmi; no sort,
    // deduplication, index table, or extra row materialization is introduced.
    const SHEET_RANGE_TAGS: [&str; 10] = [
        "sheet_id",
        "id",
        "beg_auth_asym_id",
        "beg_label_comp_id",
        "beg_auth_seq_id",
        "?pdbx_beg_PDB_ins_code",
        "end_auth_asym_id",
        "end_label_comp_id",
        "end_auth_seq_id",
        "?pdbx_end_PDB_ins_code",
    ];
    const SHEET_ORDER_TAGS: [&str; 3] = ["sheet_id", "range_id_2", "sense"];
    const SHEET_HBOND_TAGS: [&str; 12] = [
        "sheet_id",
        "range_id_2",
        "range_1_auth_asym_id",
        "range_1_label_comp_id",
        "range_1_auth_seq_id",
        "?range_1_PDB_ins_code",
        "range_1_label_atom_id",
        "range_2_auth_asym_id",
        "range_2_label_comp_id",
        "range_2_auth_seq_id",
        "?range_2_PDB_ins_code",
        "range_2_label_atom_id",
    ];

    let mut sheets = Vec::new();
    if let Some(sheet_ids) = block.find_values("_struct_sheet.id") {
        for sheet_id in sheet_ids.iter() {
            sheets.push(BioSheet {
                name: sheet_id.decoded(),
                strands: Vec::new(),
            });
        }
    }

    for row in block
        .find("_struct_sheet_range.", &SHEET_RANGE_TAGS)?
        .iter()
    {
        let sheet_id = row.decoded(0).unwrap_or_default();
        let sheet_index =
            if let Some(index) = sheets.iter().position(|sheet| sheet.name == sheet_id) {
                index
            } else {
                sheets.push(BioSheet {
                    name: sheet_id,
                    strands: Vec::new(),
                });
                sheets.len() - 1
            };
        let strand_name = row.decoded(1).unwrap_or_default();
        let start = make_mmcif_sheet_address(&row, 2, 3, 4, 5, None, "start")?;
        let end = make_mmcif_sheet_address(&row, 6, 7, 8, 9, None, "end")?;
        sheets[sheet_index].strands.push(BioStrand::new(
            start,
            end,
            AtomAddress::default(),
            AtomAddress::default(),
            0,
            strand_name,
        ));
    }

    for row in block
        .find("_struct_sheet_order.", &SHEET_ORDER_TAGS)?
        .iter()
    {
        let sheet_id = row.decoded(0).unwrap_or_default();
        let strand_id = row.decoded(1).unwrap_or_default();
        if let Some(sheet) = sheets.iter_mut().find(|sheet| sheet.name == sheet_id)
            && let Some(strand) = sheet
                .strands
                .iter_mut()
                .find(|strand| strand.name == strand_id)
        {
            let sense = row
                .decoded(2)
                .unwrap_or_default()
                .as_bytes()
                .first()
                .copied()
                .unwrap_or_default()
                & !0x20;
            match sense {
                b'P' => strand.sense = 1,
                b'A' => strand.sense = -1,
                _ => {}
            }
        }
    }

    for row in block
        .find("_pdbx_struct_sheet_hbond.", &SHEET_HBOND_TAGS)?
        .iter()
    {
        let sheet_id = row.decoded(0).unwrap_or_default();
        let strand_id = row.decoded(1).unwrap_or_default();
        if let Some(sheet) = sheets.iter_mut().find(|sheet| sheet.name == sheet_id)
            && let Some(strand) = sheet
                .strands
                .iter_mut()
                .find(|strand| strand.name == strand_id)
        {
            let hbond_atom1 = make_mmcif_sheet_address(&row, 2, 3, 4, 5, Some(6), "hbond_atom1")?;
            let hbond_atom2 = make_mmcif_sheet_address(&row, 7, 8, 9, 10, Some(11), "hbond_atom2")?;
            strand.hbond_atom1 = hbond_atom1;
            strand.hbond_atom2 = hbond_atom2;
        }
    }

    Ok(sheets)
}

fn parse_mmcif_software_info(
    block: &CifBlock,
    metadata: &mut BioMetadata,
) -> Result<(), CifReadError> {
    // Gemmi❗✔️: void read_software_info(cif::Block& block, Structure& st) {
    // Gemmi❗✔️:     for (auto row : block.find("_software.", {"name",
    // Gemmi❗✔️:                                               "?classification",
    // Gemmi❗✔️:                                               "?version",
    // Gemmi❗✔️:                                               "?date",
    // Gemmi❗✔️:                                               "?description",
    // Gemmi❗✔️:                                               "?contact_author",
    // Gemmi❗✔️:                                               "?contact_author_email"})) {
    // Gemmi❗✔️:         st.meta.software.emplace_back();
    // Gemmi❗✔️:         SoftwareItem& item = st.meta.software.back();
    // Gemmi❗✔️:         item.name = row.str(0);
    // Gemmi❗✔️:         if (row.has2(1))
    // Gemmi❗✔️:             item.classification = software_classification_from_string(row.str(1));
    // Gemmi❗✔️:         copy_string(row, 2, item.version);
    // Gemmi❗✔️:         copy_string(row, 3, item.date);
    // Gemmi❗✔️:         copy_string(row, 4, item.description);
    // Gemmi❗✔️:         copy_string(row, 5, item.contact_author);
    // Gemmi❗✔️:         copy_string(row, 6, item.contact_author_email);
    // Gemmi❗✔️:     }
    // Gemmi❗✔️: }
    // Gemmi❗✔️: void copy_string(const cif::Table::Row& row, int n, std::string& dest) {
    // Gemmi❗✔️:   if (row.has2(n))
    // Gemmi❗✔️:     dest = cif::as_string(row[n]);
    // Gemmi❗✔️: }
    // Gemmi❗✔️: bool has(size_t n) const { return tab.positions.at(n) >= 0; }
    // Gemmi❗✔️: bool has2(size_t n) const { return has(n) && !cif::is_null(operator[](n)); }
    // Gemmi❗✔️: std::string str(int n) const { return as_string(at(n)); }
    // Gemmi❗✔️: inline SoftwareItem::Classification
    // Gemmi❗✔️: software_classification_from_string(const std::string& str) {
    // Gemmi❗✔️:   if (iequal(str, "data collection")) return SoftwareItem::DataCollection;
    // Gemmi❗✔️:   if (iequal(str, "data extraction")) return SoftwareItem::DataExtraction;
    // Gemmi❗✔️:   if (iequal(str, "data processing")) return SoftwareItem::DataProcessing;
    // Gemmi❗✔️:   if (iequal(str, "data reduction"))  return SoftwareItem::DataReduction;
    // Gemmi❗✔️:   if (iequal(str, "data scaling"))    return SoftwareItem::DataScaling;
    // Gemmi❗✔️:   if (iequal(str, "model building"))  return SoftwareItem::ModelBuilding;
    // Gemmi❗✔️:   if (iequal(str, "phasing"))         return SoftwareItem::Phasing;
    // Gemmi❗✔️:   if (iequal(str, "refinement"))      return SoftwareItem::Refinement;
    // Gemmi❗✔️:   return SoftwareItem::Unspecified;
    // Gemmi❗✔️: }
    // Gemmi❗✔️: inline bool iequal_from(const std::string& str, size_t offset, const std::string& low) {
    // Gemmi❗✔️:   return str.length() == low.length() + offset &&
    // Gemmi❗✔️:          std::equal(std::begin(low), std::end(low), str.begin() + offset,
    // Gemmi❗✔️:                     [](char c1, char c2) { return c1 == lower(c2); });
    // Gemmi❗✔️: }
    // Gemmi❗✔️: inline bool iequal(const std::string& str, const std::string& low) {
    // Gemmi❗✔️:   return iequal_from(str, 0, low);
    // Gemmi❗✔️: }
    // Gemmi❗✔️: inline char lower(char c) {
    // Gemmi❗✔️:   if (c >= 'A' && c <= 'Z')
    // Gemmi❗✔️:     return c | 0x20;
    // Gemmi❗✔️:   return c;
    // Gemmi❗✔️: }
    // Behavior review: use Gemmi's exact required/optional table columns,
    // preserve table row order and append one default BIO value per row.
    // Missing/null classification leaves Unspecified; recognized names use
    // Gemmi's ASCII case-insensitive mapping and unknown non-null names also
    // map to Unspecified. Optional strings are assigned only when present and
    // non-null, with source CIF decoding; absent/null values retain defaults.
    // Complexity review: one canonical table-selection scan and one pass over
    // rows, with a fixed classification comparison set and at most five
    // decoded optional strings per row; this matches source lookup/row costs.
    const SOFTWARE_TAGS: [&str; 7] = [
        "name",
        "?classification",
        "?version",
        "?date",
        "?description",
        "?contact_author",
        "?contact_author_email",
    ];
    let table = block.find("_software.", &SOFTWARE_TAGS)?;
    for row in table.iter() {
        let mut item = BioSoftwareItem::default();
        item.name = row.decoded(0).unwrap_or_default();

        if let Some(value) = row.get(1).filter(|value| !value.is_null()) {
            let classification = value.decoded();
            item.classification = if classification.eq_ignore_ascii_case("data collection") {
                BioSoftwareClassification::DataCollection
            } else if classification.eq_ignore_ascii_case("data extraction") {
                BioSoftwareClassification::DataExtraction
            } else if classification.eq_ignore_ascii_case("data processing") {
                BioSoftwareClassification::DataProcessing
            } else if classification.eq_ignore_ascii_case("data reduction") {
                BioSoftwareClassification::DataReduction
            } else if classification.eq_ignore_ascii_case("data scaling") {
                BioSoftwareClassification::DataScaling
            } else if classification.eq_ignore_ascii_case("model building") {
                BioSoftwareClassification::ModelBuilding
            } else if classification.eq_ignore_ascii_case("phasing") {
                BioSoftwareClassification::Phasing
            } else if classification.eq_ignore_ascii_case("refinement") {
                BioSoftwareClassification::Refinement
            } else {
                BioSoftwareClassification::Unspecified
            };
        }

        if let Some(value) = row.get(2).filter(|value| !value.is_null()) {
            item.version = value.decoded();
        }
        if let Some(value) = row.get(3).filter(|value| !value.is_null()) {
            item.date = value.decoded();
        }
        if let Some(value) = row.get(4).filter(|value| !value.is_null()) {
            item.description = value.decoded();
        }
        if let Some(value) = row.get(5).filter(|value| !value.is_null()) {
            item.contact_author = value.decoded();
        }
        if let Some(value) = row.get(6).filter(|value| !value.is_null()) {
            item.contact_author_email = value.decoded();
        }

        metadata.software.push(item);
    }
    Ok(())
}

#[derive(Debug, Clone, PartialEq)]
struct MmcifNcsOrigxInfo {
    source_state: BioStructureSourceState,
    ncs_operators: Vec<BioNcsOperator>,
}

fn parse_mmcif_ncs_origx_info(
    block: &CifBlock,
) -> Result<MmcifNcsOrigxInfo, MmcifCrystalInfoError> {
    let mut source_state = BioStructureSourceState::default();
    let ncs_operators = parse_mmcif_ncs_operators(block, &mut source_state)?;
    parse_mmcif_origx(block, &mut source_state)?;

    Ok(MmcifNcsOrigxInfo {
        source_state,
        ncs_operators,
    })
}

fn parse_mmcif_ncs_operators(
    block: &CifBlock,
    source_state: &mut BioStructureSourceState,
) -> Result<Vec<BioNcsOperator>, MmcifCrystalInfoError> {
    // Gemmi❗✔️: void read_ncs_info(cif::Block& block, Structure& st) {
    // Gemmi❗✔️:     std::vector<std::string> ncs_oper_tags = transform_tags("matrix", "vector");
    // Gemmi❗✔️:     ncs_oper_tags.emplace_back("id");  // 12
    // Gemmi❗✔️:     ncs_oper_tags.emplace_back("?code");  // 13
    // Gemmi❗✔️:     cif::Table ncs_oper = block.find("_struct_ncs_oper.", ncs_oper_tags);
    // Gemmi❗✔️:     for (auto op : ncs_oper) {
    // Gemmi❗✔️:         bool given = op.has(13) && op.str(13) == "given";
    // Gemmi❗✔️:         Transform tr = get_transform_matrix(op);
    // Gemmi❗✔️:         if (tr.is_identity())
    // Gemmi❗✔️:             // ignore identity, but store its id so we can write it back to mmCIF
    // Gemmi❗✔️:             st.info["_struct_ncs_oper.id"] = op.str(12);
    // Gemmi❗✔️:         else if (tr.has_nan())
    // Gemmi❗✔️:             // As of 2022 some entries (7qb5, 6tsd) have incomplete _struct_ncs_oper.
    // Gemmi❗✔️:             // It is safer to skip them.
    // Gemmi❗✔️:             continue;
    // Gemmi❗✔️:         else
    // Gemmi❗✔️:             st.ncs.push_back({op.str(12), given, tr});
    // Gemmi❗✔️:     }
    // Gemmi❗✔️: }
    // Behavior review: NCS identity and NaN checks retain source order and
    // exact component semantics; identity ids overwrite the same metadata key,
    // while non-identity finite/non-NaN transforms append in table order.
    // Complexity review: one linear pass over NCS rows, a fixed 12-value
    // conversion per row, amortized vector appends and source-equivalent
    // ordered-map insertion; ORIGX is intentionally applied by the following
    // source-ordered helper.
    let mut ncs_operators = Vec::new();

    let mut ncs_tags = mmcif_transform_tags("matrix", "vector");
    ncs_tags.push("id".to_owned());
    ncs_tags.push("?code".to_owned());
    let ncs_tag_refs: Vec<&str> = ncs_tags.iter().map(String::as_str).collect();
    let ncs_table = block.find("_struct_ncs_oper.", &ncs_tag_refs)?;
    for row in ncs_table.iter() {
        let given = row.has(13) && row.decoded(13).as_deref() == Some("given");
        let transform = mmcif_transform_from_row(row)?;
        let id = row.decoded(12).unwrap_or_default();
        if mmcif_transform_is_identity(&transform) {
            source_state
                .info
                .insert("_struct_ncs_oper.id".to_owned(), id);
        } else if mmcif_transform_has_nan(&transform) {
            continue;
        } else {
            ncs_operators.push(BioNcsOperator::new(id, given, transform));
        }
    }

    Ok(ncs_operators)
}

fn parse_mmcif_origx(
    block: &CifBlock,
    source_state: &mut BioStructureSourceState,
) -> Result<(), MmcifCrystalInfoError> {
    // Gemmi❗✔️:   cif::Table origx_tv = block.find("_database_PDB_matrix.",
    // Gemmi❗✔️:                                    transform_tags("origx", "origx_vector"));
    // Gemmi❗✔️:   if (origx_tv.length() > 0) {
    // Gemmi❗✔️:     st.has_origx = true;
    // Gemmi❗✔️:     st.origx = get_transform_matrix(origx_tv[0]);
    // Gemmi❗✔️:   }
    // Behavior review: ORIGX is read only after the fractional transform in
    // the full reader, requires the complete required-tag table and at least
    // one row, and retains only row zero.
    // Complexity review: one fixed table lookup and at most one 12-field row
    // conversion; no additional model traversal or temporary collection.
    let origx_tags = mmcif_transform_tags("origx", "origx_vector");
    let origx_tag_refs: [&str; 12] = std::array::from_fn(|index| origx_tags[index].as_str());
    let origx_table = block.find("_database_PDB_matrix.", &origx_tag_refs)?;
    if origx_table.len() > 0 {
        let row = origx_table
            .row(0)
            .ok_or(MmcifCrystalInfoError::StagingInvariant {
                stage: "nonempty ORIGX table lacks row zero",
            })?;
        source_state.has_origx = true;
        source_state.origx = mmcif_transform_from_row(row)?;
    }
    Ok(())
}

fn mmcif_transform_is_identity(transform: &BioTransform) -> bool {
    // Gemmi❗✔️: bool is_identity() const {
    // Gemmi❗✔️:   return mat.is_identity() && vec.x == 0. && vec.y == 0. && vec.z == 0.;
    // Gemmi❗✔️: }
    // Gemmi❗✔️: bool is_identity() const {
    // Gemmi❗✔️:   return a[0][0] == 1 && a[0][1] == 0 && a[0][2] == 0 &&
    // Gemmi❗✔️:          a[1][0] == 0 && a[1][1] == 1 && a[1][2] == 0 &&
    // Gemmi❗✔️:          a[2][0] == 0 && a[2][1] == 0 && a[2][2] == 1;
    // Gemmi❗✔️: }
    // Behavior review: exact f64 equality, including IEEE signed-zero equality,
    // matches Gemmi's matrix and translation comparisons without a tolerance.
    // Complexity review: twelve fixed scalar comparisons, with no allocation.
    *transform.matrix() == *BioTransform::identity().matrix()
        && *transform.translation() == [0.0; 3]
}

fn mmcif_transform_has_nan(transform: &BioTransform) -> bool {
    // Gemmi❗✔️: bool has_nan() const {
    // Gemmi❗✔️:   for (int i = 0; i < 3; ++i)
    // Gemmi❗✔️:     for (int j = 0; j < 3; ++j)
    // Gemmi❗✔️:       if (std::isnan(a[i][j]))
    // Gemmi❗✔️:           return true;
    // Gemmi❗✔️:     return false;
    // Gemmi❗✔️: }
    // Gemmi❗✔️: bool has_nan() const {
    // Gemmi❗✔️:   return std::isnan(x) || std::isnan(y) || std::isnan(z);
    // Gemmi❗✔️: }
    // Gemmi❗✔️: bool has_nan() const {
    // Gemmi❗✔️:   return mat.has_nan() || vec.has_nan();
    // Gemmi❗✔️: }
    // Behavior review: scans row-major matrix then translation and stops at
    // the first NaN, matching Gemmi; infinities are not classified as NaN.
    // Complexity review: at most nine matrix and three vector checks, no heap
    // allocation or temporary collection.
    transform
        .matrix()
        .iter()
        .flatten()
        .any(|component| component.is_nan())
        || transform
            .translation()
            .iter()
            .any(|component| component.is_nan())
}

fn mmcif_transform_tags(matrix: &str, vector: &str) -> Vec<String> {
    // Gemmi❗✔️: std::vector<std::string> transform_tags(const std::string& mstr, const std::string& vstr) {
    // Gemmi❗✔️:   return {mstr + "[1][1]", mstr + "[1][2]", mstr + "[1][3]", vstr + "[1]",
    // Gemmi❗✔️:           mstr + "[2][1]", mstr + "[2][2]", mstr + "[2][3]", vstr + "[2]",
    // Gemmi❗✔️:           mstr + "[3][1]", mstr + "[3][2]", mstr + "[3][3]", vstr + "[3]"};
    // Gemmi❗✔️: }
    // Behavior review: the 12 matrix/vector tags retain Gemmi's row-major
    // order; later transform readers share this one fixed mapping.
    // Complexity review: constructs the same fixed 12 owned tag strings as
    // Gemmi's returned vector.
    vec![
        format!("{matrix}[1][1]"),
        format!("{matrix}[1][2]"),
        format!("{matrix}[1][3]"),
        format!("{vector}[1]"),
        format!("{matrix}[2][1]"),
        format!("{matrix}[2][2]"),
        format!("{matrix}[2][3]"),
        format!("{vector}[2]"),
        format!("{matrix}[3][1]"),
        format!("{matrix}[3][2]"),
        format!("{matrix}[3][3]"),
        format!("{vector}[3]"),
    ]
}

fn mmcif_transform_from_row(row: CifRow<'_>) -> Result<BioTransform, MmcifCrystalInfoError> {
    // Gemmi❗✔️: Transform get_transform_matrix(const cif::Table::Row& r) {
    // Gemmi❗✔️:   Transform t;
    // Gemmi❗✔️:   for (int i = 0; i < 3; ++i) {
    // Gemmi❗✔️:     for (int j = 0; j < 3; ++j)
    // Gemmi❗✔️:       t.mat[i][j] = cif::as_number(r[4*i+j]);
    // Gemmi❗✔️:     t.vec.at(i) = cif::as_number(r[4*i+3]);
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   return t;
    // Gemmi❗✔️: }
    // Behavior review: all 12 values use Gemmi's row-major matrix/vector
    // indexing; absent pair payloads produce the source empty string and the
    // same CIF NaN fallback. Numeric edge parity remains governed by the CIF
    // owner's documented floating-conversion boundary.
    // Complexity review: reads and converts exactly 12 scalar fields with no
    // temporary row/vector collection.
    let mut values = [f64::NAN; 12];
    for (column, value) in values.iter_mut().enumerate() {
        let raw = row.get(column).map(CifValue::raw).unwrap_or("");
        *value = cif_as_f64(raw, f64::NAN)?;
    }
    Ok(BioTransform::new(
        [
            [values[0], values[1], values[2]],
            [values[4], values[5], values[6]],
            [values[8], values[9], values[10]],
        ],
        [values[3], values[7], values[11]],
    ))
}

struct MmcifEntityBuildState {
    source_entity_id: String,
    kind: EntityKind,
    polymer_kind: PolymerKind,
    full_sequence: Vec<String>,
    dbrefs: Vec<BioEntityDbRef>,
    subchains: Vec<String>,
    sifts_unp_accessions: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum MmcifSiftsUnpError {
    Cif(CifReadError),
    EntityNotFound { entity_id: String },
    SubchainNotFound { asym_id: String },
    SequenceNotFound { seq_id: String },
    NumberOutsideU16 { value: String },
    StagingInvariant { stage: &'static str },
}

impl From<CifReadError> for MmcifSiftsUnpError {
    fn from(error: CifReadError) -> Self {
        Self::Cif(error)
    }
}

impl fmt::Display for MmcifSiftsUnpError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Cif(error) => fmt::Display::fmt(error, f),
            Self::EntityNotFound { entity_id } => {
                write!(f, "_pdbx_sifts_xref_db: entity_id not found: {entity_id}")
            }
            Self::SubchainNotFound { asym_id } => {
                write!(f, "_pdbx_sifts_xref_db: asym_id not found: {asym_id}")
            }
            Self::SequenceNotFound { seq_id } => {
                write!(f, "_pdbx_sifts_xref_db: seq_id not found: {seq_id}")
            }
            Self::NumberOutsideU16 { value } => {
                write!(f, "_pdbx_sifts_xref_db.unp_num: {value}")
            }
            Self::StagingInvariant { stage } => {
                write!(f, "internal mmCIF SIFTS staging invariant failed: {stage}")
            }
        }
    }
}

impl std::error::Error for MmcifSiftsUnpError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Cif(error) => Some(error),
            Self::EntityNotFound { .. }
            | Self::SubchainNotFound { .. }
            | Self::SequenceNotFound { .. }
            | Self::NumberOutsideU16 { .. }
            | Self::StagingInvariant { .. } => None,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum MmcifEntityInfoError {
    Cif(CifReadError),
    SequenceId(MmcifSequenceIdError),
    Sifts(MmcifSiftsUnpError),
}

impl MmcifEntityInfoError {
    fn kind(&self) -> CifReadErrorKind {
        match self {
            Self::Cif(error) => error.kind(),
            Self::SequenceId(MmcifSequenceIdError::InvalidInsertionCode(error))
            | Self::SequenceId(MmcifSequenceIdError::InvalidSequenceNumber(error)) => error.kind(),
            Self::SequenceId(MmcifSequenceIdError::InconsistentInsertionCode { .. }) => {
                CifReadErrorKind::InvalidValue
            }
            Self::Sifts(MmcifSiftsUnpError::Cif(error)) => error.kind(),
            Self::Sifts(_) => CifReadErrorKind::InvalidValue,
        }
    }
}

impl From<CifReadError> for MmcifEntityInfoError {
    fn from(error: CifReadError) -> Self {
        Self::Cif(error)
    }
}

impl From<MmcifSequenceIdError> for MmcifEntityInfoError {
    fn from(error: MmcifSequenceIdError) -> Self {
        Self::SequenceId(error)
    }
}

impl From<MmcifSiftsUnpError> for MmcifEntityInfoError {
    fn from(error: MmcifSiftsUnpError) -> Self {
        Self::Sifts(error)
    }
}

impl fmt::Display for MmcifEntityInfoError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Cif(error) => fmt::Display::fmt(error, f),
            Self::SequenceId(error) => fmt::Display::fmt(error, f),
            Self::Sifts(error) => fmt::Display::fmt(error, f),
        }
    }
}

impl std::error::Error for MmcifEntityInfoError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Cif(error) => Some(error),
            Self::SequenceId(error) => Some(error),
            Self::Sifts(error) => Some(error),
        }
    }
}

fn parse_mmcif_entity_polymer_values(
    block: &CifBlock,
    mut models: Option<&mut MmcifAtomSiteGrouping>,
) -> Result<Vec<BioEntityRow>, MmcifEntityInfoError> {
    let mut entities = parse_mmcif_entity_build_state(block, models.as_deref())?;
    if let Some(grouping) = models.as_deref_mut() {
        parse_mmcif_sifts_unp(block, &mut entities, &mut grouping.models)?;
    }
    Ok(mmcif_entity_rows(entities))
}

fn parse_mmcif_entity_build_state(
    block: &CifBlock,
    models: Option<&MmcifAtomSiteGrouping>,
) -> Result<Vec<MmcifEntityBuildState>, MmcifEntityInfoError> {
    // Gemmi❗✔️: void read_entity_and_sequence_info(cif::Block& block, Structure& st) {
    // Gemmi❗✔️:     cif::Table polymer_types = block.find("_entity_poly.", {"entity_id", "type"});
    // Gemmi❗✔️:     for (auto row : block.find("_entity.", {"id", "?type"})) {
    // Gemmi❗✔️:         Entity ent(row.str(0));
    // Gemmi❗✔️:         if (row.has(1))
    // Gemmi❗✔️:             ent.entity_type = entity_type_from_string(row.str(1));
    // Gemmi❗✔️:         ent.polymer_type = PolymerType::Unknown;
    // Gemmi❗✔️:         if (polymer_types.ok()) {
    // Gemmi❗✔️:             try {
    // Gemmi❗✔️:                 std::string poly_type = polymer_types.find_row(ent.name).str(1);
    // Gemmi❗✔️:                 if (ent.entity_type == EntityType::Unknown)
    // Gemmi❗✔️:                     ent.entity_type = EntityType::Polymer;
    // Gemmi❗✔️:                 ent.polymer_type = polymer_type_from_string(poly_type);
    // Gemmi❗✔️:             } catch (std::runtime_error&) {}
    // Gemmi❗✔️:         }
    // Gemmi❗✔️:         // _entity_poly_seq is supposed to reflect heterogeneities in _atom_site.
    // Gemmi❗✔️:         ent.reflects_microhetero = true;
    // Gemmi❗✔️:         st.entities.push_back(ent);
    // Gemmi❗✔️:     }
    // Gemmi❗✔️: for (auto row : block.find("_entity_poly_seq.",
    // Gemmi❗✔️:                            {"entity_id", "num", "mon_id"}))
    // Gemmi❗✔️:   if (Entity* ent = st.get_entity(row.str(0))) {
    // Gemmi❗✔️:     int pos = cif::as_int(row[1], 0) - 1;
    // Gemmi❗✔️:     if (pos == (int) ent->full_sequence.size())
    // Gemmi❗✔️:       ent->full_sequence.push_back(row.str(2));
    // Gemmi❗✔️:     else if (pos >= 0 && pos < (int) ent->full_sequence.size())
    // Gemmi❗✔️:       cat_to(ent->full_sequence[pos], ',', row.str(2));
    // Gemmi❗✔️:   }
    // Gemmi❗✔️: Entity* get_entity(const std::string& ent_id) {
    // Gemmi❗✔️:   return impl::find_or_null(entities, ent_id);
    // Gemmi❗✔️: }
    // Gemmi❗✔️: template<typename Vec, typename S>
    // Gemmi❗✔️: auto find_iter_(Vec& vec, const S& name) {
    // Gemmi❗✔️:   return std::find_if(vec.begin(), vec.end(), [&name](const auto& m) { return get_id(m) == name; });
    // Gemmi❗✔️: }
    // Gemmi❗✔️: template<typename T, typename S>
    // Gemmi❗✔️: T* find_or_null(std::vector<T>& vec, const S& name) {
    // Gemmi❗✔️:   auto it = find_iter_(vec, name);
    // Gemmi❗✔️:   return it != vec.end() ? &*it : nullptr;
    // Gemmi❗✔️: }
    // Gemmi❗✔️: template<typename T, typename... Args>
    // Gemmi❗✔️: void cat_to(std::string& out, const T& value, Args const&... args) {
    // Gemmi❗✔️:   append_to_str(out, value);
    // Gemmi❗✔️:   cat_to(out, args...);
    // Gemmi❗✔️: }
    // Gemmi❗✔️: template<typename T>
    // Gemmi❗✔️: void append_to_str(std::string& out, const T& v) { out += v; }
    // Gemmi❗✔️: std::string str(int n) const { return as_string(at(n)); }
    // Gemmi❗✔️: inline int as_int(const std::string& str, int null) {
    // Gemmi❗✔️:   return is_null(str) ? null : as_int(str);
    // Gemmi❗✔️: }
    // Gemmi❗✔️: cif::Table struct_ref = block.find("_struct_ref.",
    // Gemmi❗✔️:     {"id", "entity_id", "db_name", "db_code",
    // Gemmi❗✔️:      "?pdbx_db_accession", "?pdbx_db_isoform"});
    // Gemmi❗✔️: cif::Table struct_ref_seq = block.find("_struct_ref_seq.",
    // Gemmi❗✔️:     {"ref_id", "seq_align_beg", "seq_align_end",                   // 0-2
    // Gemmi❗✔️:      "db_align_beg", "db_align_end",                               // 3-4
    // Gemmi❗✔️:      "?pdbx_auth_seq_align_beg", "?pdbx_seq_align_beg_ins_code",   // 5-6
    // Gemmi❗✔️:      "?pdbx_auth_seq_align_end", "?pdbx_seq_align_end_ins_code"}); // 7-8
    // Gemmi❗✔️: // DbRef doesn't correspond 1:1 to the mmCIF tables; we need to remove
    // Gemmi❗✔️: // duplicates from _struct_ref_seq to make it work.
    // Gemmi❗✔️: std::vector<std::string> seen;
    // Gemmi❗✔️: for (cif::Table::Row seq : struct_ref_seq) {
    // Gemmi❗✔️:   std::string str = seq[0];
    // Gemmi❗✔️:   for (int i = 1; i < 5; ++i) {
    // Gemmi❗✔️:     str += '\t';
    // Gemmi❗✔️:     str += seq[i];
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   if (in_vector(str, seen))
    // Gemmi❗✔️:     continue;
    // Gemmi❗✔️:   seen.push_back(str);
    // Gemmi❗✔️:   cif::Table::Row row = struct_ref.find_row(seq.str(0));
    // Gemmi❗✔️:   if (Entity* ent = st.get_entity(row.str(1))) {
    // Gemmi❗✔️:     ent->dbrefs.emplace_back();
    // Gemmi❗✔️:     Entity::DbRef& dbref = ent->dbrefs.back();
    // Gemmi❗✔️:     dbref.db_name = row.str(2);
    // Gemmi❗✔️:     dbref.id_code = row.str(3);
    // Gemmi❗✔️:     if (row.has(4))
    // Gemmi❗✔️:       dbref.accession_code = row.str(4);
    // Gemmi❗✔️:     if (row.has(5))
    // Gemmi❗✔️:       dbref.isoform = row.str(5);
    // Gemmi❗✔️:     constexpr int None = SeqId::OptionalNum::None;
    // Gemmi❗✔️:     dbref.label_seq_begin = cif::as_int(seq[1], None);
    // Gemmi❗✔️:     dbref.label_seq_end = cif::as_int(seq[2], None);
    // Gemmi❗✔️:     dbref.db_begin.num = cif::as_int(seq[3], None);
    // Gemmi❗✔️:     dbref.db_end.num = cif::as_int(seq[4], None);
    // Gemmi❗✔️:     if (seq.has(5))
    // Gemmi❗✔️:       dbref.seq_begin = make_seqid(seq.str(5), seq.ptr_at(6));
    // Gemmi❗✔️:     if (seq.has(7))
    // Gemmi❗✔️:       dbref.seq_end = make_seqid(seq.str(7), seq.ptr_at(8));
    // Gemmi❗✔️:   }
    // Gemmi❗✔️: }
    // Gemmi❗✔️: template <class T>
    // Gemmi❗✔️: bool in_vector(const T& x, const std::vector<T>& v) {
    // Gemmi❗✔️:   return std::find(v.begin(), v.end(), x) != v.end();
    // Gemmi❗✔️: }
    // Gemmi❗✔️: inline std::string& Table::Row::operator[](size_t n) {
    // Gemmi❗✔️:   int pos = tab.positions[n];
    // Gemmi❗✔️:   if (Loop* loop = tab.get_loop()) {
    // Gemmi❗✔️:     if (row_index == -1) // tags
    // Gemmi❗✔️:       return loop->tags[pos];
    // Gemmi❗✔️:     return loop->values[loop->width() * row_index + pos];
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   return tab.bloc.items[pos].pair[row_index == -1 ? 0 : 1];
    // Gemmi❗✔️: }
    // Gemmi❗✔️: bool has(size_t n) const { return tab.positions.at(n) >= 0; }
    // Gemmi❗✔️: std::string* ptr_at(int n) {
    // Gemmi❗✔️:   int pos = tab.positions.at(n < 0 ? n + size() : n);
    // Gemmi❗✔️:   return pos >= 0 ? &value_at(pos) : nullptr;
    // Gemmi❗✔️: }
    // Gemmi❗✔️: std::string str(int n) const { return as_string(at(n)); }
    // Gemmi❗✔️: inline Table::Row Table::find_row(const std::string& s) {
    // Gemmi❗✔️:   int pos = positions.at(0);
    // Gemmi❗✔️:   if (const Loop* loop = get_loop()) {
    // Gemmi❗✔️:     for (size_t i = 0; i < loop->values.size(); i += loop->width())
    // Gemmi❗✔️:       if (as_string(loop->values[i + pos]) == s)
    // Gemmi❗✔️:         return Row{*this, static_cast<int>(i / loop->width())};
    // Gemmi❗✔️:   } else if (as_string(bloc.items[pos].pair[1]) == s) {
    // Gemmi❗✔️:     return Row{*this, 0};
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   fail("Not found in " + *column_at_pos(pos).get_tag() + ": " + s);
    // Gemmi❗✔️: }
    // Gemmi❗✔️: Entity* get_entity(const std::string& ent_id) {
    // Gemmi❗✔️:   return impl::find_or_null(entities, ent_id);
    // Gemmi❗✔️: }
    // Behavior review: `_struct_ref_seq` order is preserved; its raw first-five
    // cell strings form the global duplicate key exactly as in the source.
    // The first decoded matching reference row is used; an absent reference
    // remains an error, while an unresolved entity is skipped. Nullable IDs
    // map to Gemmi's INT_MIN sentinel. The existing `make_mmcif_seq_id` owner
    // handles auth sequence/insertion semantics and errors.
    // Complexity review: for D sequence-reference rows, source-equivalent
    // vector scans are O(D^2) for duplicate detection, O(D*R) for reference
    // lookup, and O(D*E) for entity lookup; copied strings are retained once
    // per accepted DbRef plus once per seen raw key, with no hash conversion.
    // Behavior review: C08 retains each `_entity` row in source order, uses
    // the optional type column only when present, and suppresses only a
    // missing `_entity_poly` row. C09 applies `_entity_poly_seq` rows in source
    // order, resolves the first matching source entity, appends a monomer only
    // at the current sequence end, joins an earlier-position alternative with
    // a comma, and ignores negative positions, forward holes, and missing
    // entity IDs. Null sequence numbers become zero before the source's
    // one-based-to-zero-based subtraction. `cif_as_i32` supplies the already
    // source-anchored C-locale integer conversion; its checked overflow result
    // is explicit because upstream signed overflow is undefined. Position
    // subtraction is widened to i64, so `INT_MIN` remains a deterministic
    // negative ignored position; upstream `int` subtraction there is undefined
    // and is not claimed as parity.
    // Complexity review: entity/polymer type lookup is O(E*P); source entity
    // lookup for S sequence rows is O(S*E), matching Gemmi's first-match
    // vector scans. Temporary state and output are O(E + total sequence bytes);
    // each final BIO row is constructed once without cloning metadata.
    let polymer_types = block.find("_entity_poly.", &["entity_id", "type"])?;
    let entity_table = block.find("_entity.", &["id", "?type"])?;
    let mut entities = Vec::with_capacity(entity_table.len());

    for row in entity_table.iter() {
        let source_entity_id = row.decoded(0).unwrap_or_default();
        let mut kind = if row.has(1) {
            entity_kind_from_mmcif_type(&row.decoded(1).unwrap_or_default())
        } else {
            EntityKind::Unknown
        };
        let mut polymer_kind = PolymerKind::Unknown;

        if polymer_types.is_present()
            && let Ok(polymer_row) = polymer_types.find_row(&source_entity_id)
        {
            if kind == EntityKind::Unknown {
                kind = EntityKind::Polymer;
            }
            polymer_kind =
                polymer_kind_from_mmcif_type(&polymer_row.decoded(1).unwrap_or_default());
        }

        entities.push(MmcifEntityBuildState {
            source_entity_id,
            kind,
            polymer_kind,
            full_sequence: Vec::new(),
            dbrefs: Vec::new(),
            subchains: Vec::new(),
            sifts_unp_accessions: Vec::new(),
        });
    }

    let sequence_table = block.find("_entity_poly_seq.", &["entity_id", "num", "mon_id"])?;
    for row in sequence_table.iter() {
        let source_entity_id = row.decoded(0).unwrap_or_default();
        let Some(entity) = entities
            .iter_mut()
            .find(|entity| entity.source_entity_id == source_entity_id)
        else {
            continue;
        };

        let raw_sequence_number = row.get(1).map(|value| value.raw()).unwrap_or_default();
        let sequence_number = if cif_is_null(raw_sequence_number) {
            0_i64
        } else {
            i64::from(cif_as_i32(raw_sequence_number)?)
        };
        let position = sequence_number - 1;
        let sequence_len = entity.full_sequence.len() as i64;
        let monomer = row.decoded(2).unwrap_or_default();

        if position == sequence_len {
            entity.full_sequence.push(monomer);
        } else if position >= 0 && position < sequence_len {
            let sequence = &mut entity.full_sequence[position as usize];
            sequence.push(',');
            sequence.push_str(&monomer);
        }
    }

    let struct_ref = block.find(
        "_struct_ref.",
        &[
            "id",
            "entity_id",
            "db_name",
            "db_code",
            "?pdbx_db_accession",
            "?pdbx_db_isoform",
        ],
    )?;
    let struct_ref_seq = block.find(
        "_struct_ref_seq.",
        &[
            "ref_id",
            "seq_align_beg",
            "seq_align_end",
            "db_align_beg",
            "db_align_end",
            "?pdbx_auth_seq_align_beg",
            "?pdbx_seq_align_beg_ins_code",
            "?pdbx_auth_seq_align_end",
            "?pdbx_seq_align_end_ins_code",
        ],
    )?;
    let mut seen_dbref_keys: Vec<String> = Vec::new();
    for sequence_row in struct_ref_seq.iter() {
        let mut duplicate_key = sequence_row
            .get(0)
            .map(CifValue::raw)
            .unwrap_or_default()
            .to_owned();
        for column in 1..5 {
            duplicate_key.push('\t');
            duplicate_key.push_str(
                sequence_row
                    .get(column)
                    .map(CifValue::raw)
                    .unwrap_or_default(),
            );
        }
        if seen_dbref_keys.iter().any(|seen| seen == &duplicate_key) {
            continue;
        }
        seen_dbref_keys.push(duplicate_key);

        let reference_id = sequence_row.decoded(0).unwrap_or_default();
        let reference_row = struct_ref.find_row(&reference_id)?;
        let entity_id = reference_row.decoded(1).unwrap_or_default();
        let Some(entity) = entities
            .iter_mut()
            .find(|entity| entity.source_entity_id == entity_id)
        else {
            continue;
        };

        let optional_integer =
            |row: CifRow<'_>, column: usize| -> Result<Option<i32>, CifReadError> {
                let raw = row.get(column).map(CifValue::raw).unwrap_or_default();
                if cif_is_null(raw) {
                    Ok(None)
                } else {
                    let number = cif_as_i32(raw)?;
                    Ok((number != i32::MIN).then_some(number))
                }
            };
        let mut dbref = BioEntityDbRef {
            db_name: reference_row.decoded(2).unwrap_or_default(),
            accession_code: String::new(),
            id_code: reference_row.decoded(3).unwrap_or_default(),
            isoform: String::new(),
            seq_begin: PdbSeqId::new(i32::MIN, None),
            seq_end: PdbSeqId::new(i32::MIN, None),
            db_begin: PdbSeqId::new(i32::MIN, None),
            db_end: PdbSeqId::new(i32::MIN, None),
            label_seq_begin: optional_integer(sequence_row, 1)?,
            label_seq_end: optional_integer(sequence_row, 2)?,
        };
        if reference_row.has(4) {
            dbref.accession_code = reference_row.decoded(4).unwrap_or_default();
        }
        if reference_row.has(5) {
            dbref.isoform = reference_row.decoded(5).unwrap_or_default();
        }
        dbref.db_begin =
            PdbSeqId::new(optional_integer(sequence_row, 3)?.unwrap_or(i32::MIN), None);
        dbref.db_end = PdbSeqId::new(optional_integer(sequence_row, 4)?.unwrap_or(i32::MIN), None);
        if sequence_row.has(5) {
            dbref.seq_begin = make_mmcif_seq_id(
                sequence_row.decoded(5).unwrap_or_default(),
                sequence_row.get(6).map(CifValue::raw),
            )?;
        }
        if sequence_row.has(7) {
            dbref.seq_end = make_mmcif_seq_id(
                sequence_row.decoded(7).unwrap_or_default(),
                sequence_row.get(8).map(CifValue::raw),
            )?;
        }
        entity.dbrefs.push(dbref);
    }

    // Gemmi❗✔️:     cif::Table s_asym_table = block.find("_struct_asym.", {"id", "entity_id"});
    // Gemmi❗✔️:     if (s_asym_table.ok()) {
    // Gemmi❗✔️:         for (auto row : s_asym_table)
    // Gemmi❗✔️:             if (Entity* ent = st.get_entity(row.str(1)))
    // Gemmi❗✔️:                 ent->subchains.push_back(row.str(0));
    // Gemmi❗✔️:     } else if (!st.models.empty()) {
    // Gemmi❗✔️:         for (const Chain& chain : st.models[0].chains)
    // Gemmi❗✔️:             for (const ConstResidueSpan& sub : chain.subchains()) {
    // Gemmi❗✔️:                 const Residue& r = sub.front();
    // Gemmi❗✔️:                 if (Entity* ent = st.get_entity(r.entity_id))
    // Gemmi❗✔️:                     if (!in_vector(r.subchain, ent->subchains))
    // Gemmi❗✔️:                         ent->subchains.push_back(r.subchain);
    // Gemmi❗✔️:             }
    // Gemmi❗✔️:     }
    // Gemmi✔️✔️:   bool ok() const { return !positions.empty(); }
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
    // Gemmi✔️✔️:   std::vector<ConstResidueSpan> subchains() const {
    // Gemmi✔️✔️:     return impl::chain_subchains<ConstResidueSpan>(this);
    // Gemmi✔️✔️:   }
    // Gemmi❗✔️: template<typename T, typename Ch> std::vector<T> chain_subchains(Ch* ch) {
    // Gemmi❗✔️:   std::vector<T> v;
    // Gemmi❗✔️:   for (auto start = ch->residues.begin(); start != ch->residues.end(); ) {
    // Gemmi❗✔️:     auto end = start + 1;
    // Gemmi❗✔️:     while (end != ch->residues.end() && end->subchain == start->subchain)
    // Gemmi❗✔️:       ++end;
    // Gemmi❗✔️:     v.push_back(ch->whole().sub(start, end));
    // Gemmi❗✔️:     start = end;
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   return v;
    // Gemmi❗✔️: }
    // Gemmi❗✔️: template <class T>
    // Gemmi❗✔️: bool in_vector(const T& x, const std::vector<T>& v) {
    // Gemmi❗✔️:   return std::find(v.begin(), v.end(), x) != v.end();
    // Gemmi❗✔️: }
    // Behavior review: a present `_struct_asym` table takes precedence even
    // when it has zero rows; direct rows append in source order, retain
    // duplicates, and skip only entity IDs with no first source-order match.
    // Without an available table, fallback reads only the first encountered
    // model and groups each chain's consecutive equal subchain IDs. It resolves
    // entity ID and subchain from the run's first residue, then suppresses
    // repeats only within that destination entity. Missing source IDs map to
    // the empty strings used by Gemmi's default-initialized Residue fields.
    // Complexity review: direct processing uses the same O(A*E) first-match
    // vector lookup as Gemmi. Fallback scans residues and contiguous runs once,
    // then uses the same O(G*E + G*S) entity/duplicate scans. Unlike
    // `chain.subchains()`, it does not allocate a temporary vector of spans;
    // only source subchain values accepted into the destination are copied.
    let struct_asym = block.find("_struct_asym.", &["id", "entity_id"])?;
    if struct_asym.is_present() {
        for row in struct_asym.iter() {
            let entity_id = row.decoded(1).unwrap_or_default();
            if let Some(entity) = entities
                .iter_mut()
                .find(|entity| entity.source_entity_id == entity_id)
            {
                entity
                    .subchains
                    .push(row.decoded(0).unwrap_or_default().to_owned());
            }
        }
    } else if let Some(first_model) = models.and_then(|value| value.models.first()) {
        for chain in &first_model.chains {
            let mut start = 0;
            while start < chain.residues.len() {
                let first_residue = &chain.residues[start];
                let subchain = first_residue.source.subchain_id().unwrap_or_default();
                let entity_id = first_residue.source.label_entity_id().unwrap_or_default();
                let mut end = start + 1;
                while end < chain.residues.len()
                    && chain.residues[end].source.subchain_id().unwrap_or_default() == subchain
                {
                    end += 1;
                }

                if let Some(entity) = entities
                    .iter_mut()
                    .find(|entity| entity.source_entity_id == entity_id)
                    && !entity.subchains.iter().any(|value| value == subchain)
                {
                    entity.subchains.push(subchain.to_owned());
                }
                start = end;
            }
        }
    }

    Ok(entities)
}

fn mmcif_entity_rows(entities: Vec<MmcifEntityBuildState>) -> Vec<BioEntityRow> {
    entities
        .into_iter()
        .map(|entity| {
            BioEntityRow::new(
                entity.kind,
                entity.polymer_kind,
                true,
                entity.full_sequence,
                entity.dbrefs,
                entity.sifts_unp_accessions,
                entity.subchains,
                EntitySourceIds::new(entity.source_entity_id),
            )
        })
        .collect()
}

fn fill_mmcif_residue_entity_type(
    grouping: &mut MmcifAtomSiteGrouping,
    entities: &[MmcifEntityBuildState],
) {
    // Gemmi❗✔️: void fill_residue_entity_type(Structure& st) {
    // Gemmi❗✔️:   for (Model& model : st.models)
    // Gemmi❗✔️:     for (Chain& chain : model.chains)
    // Gemmi❗✔️:       for (ResidueSpan& sub : chain.subchains()) {
    // Gemmi❗✔️:         if (const Entity* ent = st.get_entity_of(sub)) {
    // Gemmi❗✔️:           for (Residue& res : sub)
    // Gemmi❗✔️:             res.entity_type = ent->entity_type;
    // Gemmi❗✔️:         } else {
    // Gemmi❗✔️:           // Don't attempt to distinguish Polymer, Branched and NonPolymer here.
    // Gemmi❗✔️:           // Third-party software may not use the same conventions regarding
    // Gemmi❗✔️:           // label_seq_id and label_asym_id that the PDB uses.
    // Gemmi❗✔️:           for (Residue& res : sub)
    // Gemmi❗✔️:             res.entity_type = res.is_water() ? EntityType::Water : EntityType::Unknown;
    // Gemmi❗✔️:         }
    // Gemmi❗✔️:       }
    // Gemmi❗✔️: }
    // Gemmi❗✔️: inline Entity* find_entity_of_subchain(const std::string& subchain_id,
    // Gemmi❗✔️:                                        std::vector<Entity>& entities) {
    // Gemmi❗✔️:   if (!subchain_id.empty())
    // Gemmi❗✔️:     for (Entity& ent : entities)
    // Gemmi❗✔️:       if (in_vector(subchain_id, ent.subchains))
    // Gemmi❗✔️:         return &ent;
    // Gemmi❗✔️:   return nullptr;
    // Gemmi❗✔️: }
    // Behavior review: each contiguous source subchain run uses the first
    // entity whose subchain list contains that exact nonempty ID; absent an
    // entity, each residue independently uses Gemmi's exact water-name test.
    // Complexity review: scans runs once and retains Gemmi's ordered linear
    // entity/subchain searches. It avoids Gemmi's temporary span vector while
    // retaining only the existing hierarchy grouping.
    for model in &mut grouping.models {
        for chain in &mut model.chains {
            let mut start = 0;
            while start < chain.residues.len() {
                let subchain = chain.residues[start]
                    .source
                    .subchain_id()
                    .unwrap_or_default();
                let mut end = start + 1;
                while end < chain.residues.len()
                    && chain.residues[end].source.subchain_id().unwrap_or_default() == subchain
                {
                    end += 1;
                }

                let entity_kind = if !subchain.is_empty() {
                    entities
                        .iter()
                        .find(|entity| entity.subchains.iter().any(|id| id == subchain))
                        .map(|entity| entity.kind)
                } else {
                    None
                };
                for residue in &mut chain.residues[start..end] {
                    residue.entity_kind = entity_kind.unwrap_or_else(|| {
                        if super::bio_pdb::gemmi_pdb_residue_is_water(residue.address.name()) {
                            EntityKind::Water
                        } else {
                            EntityKind::Unknown
                        }
                    });
                }
                start = end;
            }
        }
    }
}

fn entity_kind_from_mmcif_type(value: &str) -> EntityKind {
    // Gemmi❗✔️: inline EntityType entity_type_from_string(const std::string& t) {
    // Gemmi❗✔️:   if (t == "polymer")     return EntityType::Polymer;
    // Gemmi❗✔️:   if (t == "branched")    return EntityType::Branched;
    // Gemmi❗✔️:   if (t == "non-polymer") return EntityType::NonPolymer;
    // Gemmi❗✔️:   if (t == "water")       return EntityType::Water;
    // Gemmi❗✔️:   return EntityType::Unknown;
    // Gemmi❗✔️: }
    // Behavior review: exact case-sensitive source strings map to the
    // corresponding canonical variants; all other values remain Unknown.
    // Complexity review: the same bounded sequence of at most four string
    // comparisons is used, with no allocation.
    match value {
        "polymer" => EntityKind::Polymer,
        "branched" => EntityKind::Branched,
        "non-polymer" => EntityKind::NonPolymer,
        "water" => EntityKind::Water,
        _ => EntityKind::Unknown,
    }
}

fn polymer_kind_from_mmcif_type(value: &str) -> PolymerKind {
    // Gemmi❗✔️: inline PolymerType polymer_type_from_string(const std::string& t) {
    // Gemmi❗✔️:   if (t == "polypeptide(L)")          return PolymerType::PeptideL;
    // Gemmi❗✔️:   if (t == "polydeoxyribonucleotide") return PolymerType::Dna;
    // Gemmi❗✔️:   if (t == "polyribonucleotide")      return PolymerType::Rna;
    // Gemmi❗✔️:   if (t == "polydeoxyribonucleotide/polyribonucleotide hybrid")
    // Gemmi❗✔️:                                       return PolymerType::DnaRnaHybrid;
    // Gemmi❗✔️:   if (t == "polypeptide(D)")          return PolymerType::PeptideD;
    // Gemmi❗✔️:   if (t == "polysaccharide(D)")       return PolymerType::SaccharideD;
    // Gemmi❗✔️:   if (t == "other")                   return PolymerType::Other;
    // Gemmi❗✔️:   if (t == "peptide nucleic acid")    return PolymerType::Pna;
    // Gemmi❗✔️:   if (t == "cyclic-pseudo-peptide")   return PolymerType::CyclicPseudoPeptide;
    // Gemmi❗✔️:   if (t == "polysaccharide(L)")       return PolymerType::SaccharideL;
    // Gemmi❗✔️:   return PolymerType::Unknown;
    // Gemmi❗✔️: }
    // Behavior review: exact source spelling, including its accepted
    // unquoted decoded hybrid/PNA values, is required; unknown spellings do
    // not gain a guessed polymer kind.
    // Complexity review: the source-ordered bounded comparisons stay O(1)
    // with no temporary strings.
    match value {
        "polypeptide(L)" => PolymerKind::PeptideL,
        "polydeoxyribonucleotide" => PolymerKind::Dna,
        "polyribonucleotide" => PolymerKind::Rna,
        "polydeoxyribonucleotide/polyribonucleotide hybrid" => PolymerKind::DnaRnaHybrid,
        "polypeptide(D)" => PolymerKind::PeptideD,
        "polysaccharide(D)" => PolymerKind::SaccharideD,
        "other" => PolymerKind::Other,
        "peptide nucleic acid" => PolymerKind::Pna,
        "cyclic-pseudo-peptide" => PolymerKind::CyclicPseudoPeptide,
        "polysaccharide(L)" => PolymerKind::SaccharideL,
        _ => PolymerKind::Unknown,
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum AtomSiteIdentityField {
    AsymId,
    CompId,
    AtomId,
    SeqId,
}

impl fmt::Display for AtomSiteIdentityField {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let message = match self {
            Self::AsymId => "Neither _atom_site.label_asym_id nor auth_asym_id found",
            Self::CompId => "Neither _atom_site.label_comp_id nor auth_comp_id found",
            Self::AtomId => "Neither _atom_site.label_atom_id nor auth_atom_id found",
            Self::SeqId => "Neither _atom_site.label_seq_id nor auth_seq_id found",
        };
        f.write_str(message)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct MissingAtomSiteIdentityField(AtomSiteIdentityField);

impl fmt::Display for MissingAtomSiteIdentityField {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt::Display::fmt(&self.0, f)
    }
}

impl std::error::Error for MissingAtomSiteIdentityField {}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct AtomSiteValueColumns {
    primary: usize,
    fallback: Option<usize>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(usize)]
enum AtomSiteColumn {
    Id = 0,
    GroupPdb = 1,
    TypeSymbol = 2,
    LabelAtomId = 3,
    AltId = 4,
    LabelCompId = 5,
    LabelAsymId = 6,
    LabelEntityId = 7,
    LabelSeqId = 8,
    InsertionCode = 9,
    AuthSeqId = 16,
    AuthCompId = 17,
    AuthAsymId = 18,
    AuthAtomId = 19,
    ModelNumber = 20,
    CalcFlag = 21,
    TlsGroupId = 22,
    DeuteriumFraction = 23,
    X = 10,
    Y = 11,
    Z = 12,
    Occupancy = 13,
    BIso = 14,
    FormalCharge = 15,
}

impl AtomSiteColumn {
    const fn index(self) -> usize {
        self as usize
    }
}

impl AtomSiteValueColumns {
    fn new(
        table: &CifTable<'_>,
        auth_column: usize,
        label_column: usize,
        field: AtomSiteIdentityField,
    ) -> Result<Self, MissingAtomSiteIdentityField> {
        // Gemmi❗✔️: RowAccess(const cif::Table& tab, int n1, int n2) {
        // Gemmi❗✔️:   int pos1 = tab.positions.at(n1);
        // Gemmi❗✔️:   int pos2 = tab.positions.at(n2);
        // Gemmi❗✔️:   if (pos1 < 0) {
        // Gemmi❗✔️:     pos1 = pos2;
        // Gemmi❗✔️:     pos2 = -1;
        // Gemmi❗✔️:   }
        // Gemmi❗✔️:   const cif::Loop* loop = const_cast<cif::Table&>(tab).get_loop();
        // Gemmi❗✔️:   if (pos1 >= 0)
        // Gemmi❗✔️:     val = loop ? &loop->values[pos1] : &tab.bloc.items[pos1].pair[1];
        // Gemmi❗✔️:   if (pos2 >= 0)
        // Gemmi❗✔️:     fallback = loop ? &loop->values[pos2] : &tab.bloc.items[pos2].pair[1];
        // Gemmi❗✔️: }
        // Gemmi❗✔️: bool ok() const { return val != nullptr; }
        // Behavior review: resolve column presence once, choosing auth as
        // primary and label as per-row fallback when both are present; if
        // auth is absent, label becomes the sole primary. A missing pair is
        // an explicit source-shaped error rather than a null value.
        // Complexity review: two fixed-index presence checks and fixed-size
        // state construction, with no value allocation or row scan.
        let auth_available = table.has_column(auth_column);
        let label_available = table.has_column(label_column);
        if auth_available {
            Ok(Self {
                primary: auth_column,
                fallback: label_available.then_some(label_column),
            })
        } else if label_available {
            Ok(Self {
                primary: label_column,
                fallback: None,
            })
        } else {
            Err(MissingAtomSiteIdentityField(field))
        }
    }

    fn get<'a>(&self, row: CifRow<'a>) -> &'a str {
        // Gemmi❗✔️: const std::string& get(size_t gap) const {
        // Gemmi❗✔️:   const std::string& r = val[gap];
        // Gemmi❗✔️:   if (!cif::is_null(r) || fallback == nullptr)
        // Gemmi❗✔️:     return r;
        // Gemmi❗✔️:   return fallback[gap];
        // Gemmi❗✔️: }
        // Behavior review: return the raw auth value unless it is CIF-null
        // and a label column exists; when auth is unavailable, `new` binds
        // label as primary so its raw null token is retained. Complexity
        // review: fixed-index row lookup and null check are O(1), with no
        // per-value allocation, matching RowAccess pointer indexing.
        let value = row.get(self.primary).map(CifValue::raw).unwrap_or("");
        if !cif_is_null(value) || self.fallback.is_none() {
            return value;
        }
        row.get(self.fallback.expect("checked fallback column"))
            .map(CifValue::raw)
            .unwrap_or("")
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct AtomSiteIdentityColumns {
    asym_id: AtomSiteValueColumns,
    comp_id: AtomSiteValueColumns,
    atom_id: AtomSiteValueColumns,
    seq_id: AtomSiteValueColumns,
}

fn bind_atom_site_identity_columns(
    atom_table: &CifTable<'_>,
) -> Result<AtomSiteIdentityColumns, MissingAtomSiteIdentityField> {
    // Gemmi❗✔️: enum { kId=0, kGroupPdb, kSymbol, kLabelAtomId, kAltId, kLabelCompId,
    // Gemmi❗✔️:        kLabelAsymId, kLabelEntityId, kLabelSeqId, kInsCode,
    // Gemmi❗✔️:        kX, kY, kZ, kOcc, kBiso, kCharge,
    // Gemmi❗✔️:        kAuthSeqId, kAuthCompId, kAuthAsymId, kAuthAtomId, kModelNum,
    // Gemmi❗✔️:        kCalcFlag, kTlsGroupId, kDeuterium };
    // Gemmi❗✔️: RowAccess asym_id(atom_table, kAuthAsymId, kLabelAsymId);
    // Gemmi❗✔️: RowAccess comp_id(atom_table, kAuthCompId, kLabelCompId);
    // Gemmi❗✔️: RowAccess atom_id(atom_table, kAuthAtomId, kLabelAtomId);
    // Gemmi❗✔️: RowAccess seq_id(atom_table, kAuthSeqId, kLabelSeqId);
    // Gemmi❗✔️: if (!asym_id.ok())
    // Gemmi❗✔️:   fail("Neither _atom_site.label_asym_id nor auth_asym_id found");
    // Gemmi❗✔️: if (!comp_id.ok())
    // Gemmi❗✔️:   fail("Neither _atom_site.label_comp_id nor auth_comp_id found");
    // Gemmi❗✔️: if (!atom_id.ok())
    // Gemmi❗✔️:   fail("Neither _atom_site.label_atom_id nor auth_atom_id found");
    // Gemmi❗✔️: if (!seq_id.ok())
    // Gemmi❗✔️:   fail("Neither _atom_site.label_seq_id nor auth_seq_id found");
    // Behavior review: bind auth first when that column exists, otherwise
    // label; preserve an available label as per-row fallback only when auth
    // exists. Missing pairs return the source-ordered typed error. Complexity
    // review: four constant-time column-presence checks and fixed struct
    // construction; the table owns parsing and ordered lookup.
    let asym_id = AtomSiteValueColumns::new(
        atom_table,
        AtomSiteColumn::AuthAsymId.index(),
        AtomSiteColumn::LabelAsymId.index(),
        AtomSiteIdentityField::AsymId,
    )?;
    let comp_id = AtomSiteValueColumns::new(
        atom_table,
        AtomSiteColumn::AuthCompId.index(),
        AtomSiteColumn::LabelCompId.index(),
        AtomSiteIdentityField::CompId,
    )?;
    let atom_id = AtomSiteValueColumns::new(
        atom_table,
        AtomSiteColumn::AuthAtomId.index(),
        AtomSiteColumn::LabelAtomId.index(),
        AtomSiteIdentityField::AtomId,
    )?;
    let seq_id = AtomSiteValueColumns::new(
        atom_table,
        AtomSiteColumn::AuthSeqId.index(),
        AtomSiteColumn::LabelSeqId.index(),
        AtomSiteIdentityField::SeqId,
    )?;
    Ok(AtomSiteIdentityColumns {
        asym_id,
        comp_id,
        atom_id,
        seq_id,
    })
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum MmcifSequenceIdError {
    InconsistentInsertionCode { sequence_id: String },
    InvalidInsertionCode(CifReadError),
    InvalidSequenceNumber(CifReadError),
}

impl fmt::Display for MmcifSequenceIdError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InconsistentInsertionCode { sequence_id } => {
                write!(f, "Inconsistent insertion code in {sequence_id}")
            }
            Self::InvalidInsertionCode(error) | Self::InvalidSequenceNumber(error) => {
                fmt::Display::fmt(error, f)
            }
        }
    }
}

impl std::error::Error for MmcifSequenceIdError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::InconsistentInsertionCode { .. } => None,
            Self::InvalidInsertionCode(error) | Self::InvalidSequenceNumber(error) => Some(error),
        }
    }
}

fn make_mmcif_seq_id(
    sequence_id: String,
    raw_insertion_code: Option<&str>,
) -> Result<PdbSeqId, MmcifSequenceIdError> {
    // Gemmi❗✔️: SeqId ret;
    // Gemmi❗✔️: if (icode)
    // Gemmi❗✔️:   // the insertion code happens to be always a single letter
    // Gemmi❗✔️:   ret.icode = cif::as_char(*icode, ' ');
    // Gemmi❗✔️: if (!seqid.empty()) {
    // Gemmi❗✔️:   // old mmCIF files have auth_seq_id as number + icode (e.g. 15A)
    // Gemmi❗✔️:   if (seqid.back() >= 'A') {
    // Gemmi❗✔️:     if (ret.icode == ' ')
    // Gemmi❗✔️:       ret.icode = seqid.back();
    // Gemmi❗✔️:     else if (ret.icode != seqid.back())
    // Gemmi❗✔️:       fail("Inconsistent insertion code in " + seqid);
    // Gemmi❗✔️:     seqid.pop_back();
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   // 7pvv has an empty seqnum in tls description, don't throw
    // Gemmi❗✔️:   if (!seqid.empty())
    // Gemmi❗✔️:     ret.num = cif::as_int(seqid, Residue::OptionalNum::None);
    // Gemmi❗✔️: }
    // Gemmi❗✔️: return ret;
    // Behavior review: `sequence_id` is the single CIF `as_string` result, and
    // `raw_insertion_code` is the un-decoded optional row value passed to
    // `as_char`. Gemmi's byte comparison treats `char` as signed in the pinned
    // Linux reference build; only a trailing byte >= 'A' in that signed
    // domain is split. A separate non-space code must equal the embedded byte
    // exactly. `cif_as_char` preserves the source-defined NUL terminator for
    // empty decoded strings. `cif_as_i32` matches defined whitespace/sign/full
    // consumption; only source signed overflow is undefined and not claimed
    // as parity. INT_MIN remains the BIO absent-number sentinel, and space is
    // represented as `None`.
    // Complexity review: one constant-time trailing-byte check, at most one
    // linear integer scan, and one bounded CIF-char decode. The owned string
    // is consumed rather than cloned; no row scan or graph allocation occurs.
    let mut insertion_code = match raw_insertion_code {
        Some(value) => cif_as_char(value, ' ')
            .map(|character| character as u8)
            .map_err(MmcifSequenceIdError::InvalidInsertionCode)?,
        None => b' ',
    };

    let numeric_text = if let Some(last_byte) = sequence_id.as_bytes().last().copied()
        && (b'A'..0x80).contains(&last_byte)
    {
        if insertion_code == b' ' {
            insertion_code = last_byte;
        } else if insertion_code != last_byte {
            return Err(MmcifSequenceIdError::InconsistentInsertionCode { sequence_id });
        }
        &sequence_id[..sequence_id.len() - 1]
    } else {
        &sequence_id
    };

    let sequence_number = if numeric_text.is_empty() {
        i32::MIN
    } else {
        cif_as_i32(numeric_text).map_err(MmcifSequenceIdError::InvalidSequenceNumber)?
    };

    Ok(PdbSeqId::new(
        sequence_number,
        (insertion_code != b' ').then_some(insertion_code),
    ))
}

fn make_mmcif_residue_address(name: ResidueName, sequence_id: PdbSeqId) -> Option<ResidueAddress> {
    // Gemmi❗✔️: inline ResidueId make_resid(const std::string& name,
    // Gemmi❗✔️:                             const std::string& seqid,
    // Gemmi❗✔️:                             const std::string* icode) {
    // Gemmi❗✔️:   return ResidueId{make_seqid(seqid, icode), {}, name};
    // Gemmi❗✔️: }
    // Behavior review: the canonical BIO address receives the exact sequence
    // sentinel/insertion result and Gemmi's empty segment; residue spelling
    // is already represented by the caller's canonical `ResidueName`.
    // Complexity review: fixed scalar copy and the canonical constructor's
    // zero-length segment check; no string clone or traversal.
    ResidueAddress::new(
        Some(sequence_id.seq_num()),
        sequence_id.ins_code(),
        b"",
        name,
    )
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum MmcifAtomSiteGroupingError {
    IdentityField(MissingAtomSiteIdentityField),
    SequenceId(MmcifSequenceIdError),
    AtomNameNotRepresentable {
        value: String,
    },
    AtomScalar(MmcifAtomSiteScalarError),
    ModelNumber(CifReadError),
    LabelSequenceNumber(CifReadError),
    AuthorChainIdNotRepresentable {
        value: String,
    },
    ResidueNameNotRepresentable {
        value: String,
    },
    ResidueAddressNotRepresentable {
        value: String,
    },
    ResidueSourceIdsNotRepresentable,
    UndefinedGroupPdbIndex {
        value: String,
        index: usize,
        length: usize,
    },
    MissingActiveModel {
        row_index: usize,
    },
}

impl fmt::Display for MmcifAtomSiteGroupingError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::IdentityField(error) => fmt::Display::fmt(error, f),
            Self::SequenceId(error) => fmt::Display::fmt(error, f),
            Self::AtomNameNotRepresentable { value } => write!(
                f,
                "atom name is outside the approved BIO representation: {value:?}"
            ),
            Self::AtomScalar(error) => fmt::Display::fmt(error, f),
            Self::ModelNumber(error) => write!(f, "invalid atom-site model number: {error}"),
            Self::LabelSequenceNumber(error) => {
                write!(f, "invalid atom-site label sequence number: {error}")
            }
            Self::AuthorChainIdNotRepresentable { value } => write!(
                f,
                "author chain identifier exceeds the approved BIO representation: {value:?}"
            ),
            Self::ResidueNameNotRepresentable { value } => write!(
                f,
                "residue name exceeds the approved BIO representation: {value:?}"
            ),
            Self::ResidueAddressNotRepresentable { value } => write!(
                f,
                "residue address exceeds the approved BIO representation: {value:?}"
            ),
            Self::ResidueSourceIdsNotRepresentable => {
                f.write_str("mmCIF residue source identifiers are not representable")
            }
            Self::UndefinedGroupPdbIndex {
                value,
                index,
                length,
            } => write!(
                f,
                "Gemmi group_PDB byte access is undefined at index {index} for length {length}: {value:?}"
            ),
            Self::MissingActiveModel { row_index } => {
                write!(f, "atom-site row {row_index} did not select a source model")
            }
        }
    }
}

impl std::error::Error for MmcifAtomSiteGroupingError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::IdentityField(error) => Some(error),
            Self::SequenceId(error) => Some(error),
            Self::AtomScalar(error) => Some(error),
            Self::ModelNumber(error) | Self::LabelSequenceNumber(error) => Some(error),
            Self::AtomNameNotRepresentable { .. }
            | Self::AuthorChainIdNotRepresentable { .. }
            | Self::ResidueNameNotRepresentable { .. }
            | Self::ResidueAddressNotRepresentable { .. }
            | Self::ResidueSourceIdsNotRepresentable
            | Self::UndefinedGroupPdbIndex { .. }
            | Self::MissingActiveModel { .. } => None,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct MmcifResidueGrouping {
    address: ResidueAddress,
    source: ResidueSourceIds,
    het_flag: Option<u8>,
    entity_kind: EntityKind,
    atom_site_rows: Vec<usize>,
    sifts_unp: BioSiftsUnpResidue,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct MmcifChainGrouping {
    source_name: String,
    source: ChainSourceIds,
    residues: Vec<MmcifResidueGrouping>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct MmcifModelGrouping {
    source_model_number: i32,
    chains: Vec<MmcifChainGrouping>,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
struct MmcifAtomSiteGrouping {
    models: Vec<MmcifModelGrouping>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum MmcifLabelAddressError {
    LabelSequenceId(CifReadError),
    ChainNameNotRepresentable { value: String },
}

impl fmt::Display for MmcifLabelAddressError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::LabelSequenceId(error) => fmt::Display::fmt(error, f),
            Self::ChainNameNotRepresentable { value } => write!(
                f,
                "mmCIF label-resolved chain name exceeds the approved BIO representation: {value:?}"
            ),
        }
    }
}

impl std::error::Error for MmcifLabelAddressError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::LabelSequenceId(error) => Some(error),
            Self::ChainNameNotRepresentable { .. } => None,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum MmcifConnectionError {
    Cif(CifReadError),
    SequenceId(MmcifSequenceIdError),
    LabelAddress(MmcifLabelAddressError),
    MissingAddressIdentifiers { partner: usize },
    NoStructuralModels,
    ChainNameNotRepresentable { partner: usize, value: String },
    ResidueNameNotRepresentable { partner: usize, value: String },
    ResidueAddressNotRepresentable { partner: usize, value: String },
    SymmetryOperatorOutsideSourceDefinedRange { value: String },
    SymmetryOperatorNotRepresentable { value: i32 },
}

impl From<CifReadError> for MmcifConnectionError {
    fn from(error: CifReadError) -> Self {
        Self::Cif(error)
    }
}

impl fmt::Display for MmcifConnectionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Cif(error) => fmt::Display::fmt(error, f),
            Self::SequenceId(error) => fmt::Display::fmt(error, f),
            Self::LabelAddress(error) => fmt::Display::fmt(error, f),
            Self::MissingAddressIdentifiers { partner } => write!(
                f,
                "_struct_conn without either _auth_ or _label_ asym_id+seq_id for partner {partner}"
            ),
            Self::NoStructuralModels => f.write_str("no structural models"),
            Self::ChainNameNotRepresentable { partner, value } => write!(
                f,
                "_struct_conn partner {partner} chain identifier exceeds the approved BIO representation: {value:?}"
            ),
            Self::ResidueNameNotRepresentable { partner, value } => write!(
                f,
                "_struct_conn partner {partner} residue name exceeds the approved BIO representation: {value:?}"
            ),
            Self::ResidueAddressNotRepresentable { partner, value } => write!(
                f,
                "_struct_conn partner {partner} residue address is not representable: {value:?}"
            ),
            Self::SymmetryOperatorOutsideSourceDefinedRange { value } => write!(
                f,
                "_struct_conn symmetry operator is outside Gemmi's source-defined integer range: {value:?}"
            ),
            Self::SymmetryOperatorNotRepresentable { value } => write!(
                f,
                "_struct_conn symmetry operator {value} is outside the BIO reported-symmetry representation"
            ),
        }
    }
}

impl std::error::Error for MmcifConnectionError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Cif(error) => Some(error),
            Self::SequenceId(error) => Some(error),
            Self::LabelAddress(error) => Some(error),
            Self::MissingAddressIdentifiers { .. }
            | Self::NoStructuralModels
            | Self::ChainNameNotRepresentable { .. }
            | Self::ResidueNameNotRepresentable { .. }
            | Self::ResidueAddressNotRepresentable { .. }
            | Self::SymmetryOperatorOutsideSourceDefinedRange { .. }
            | Self::SymmetryOperatorNotRepresentable { .. } => None,
        }
    }
}

fn parse_mmcif_connections(
    block: &CifBlock,
    grouping: &MmcifAtomSiteGrouping,
) -> Result<Vec<BioConnection>, MmcifConnectionError> {
    // Gemmi❗✔️: void read_connectivity(cif::Block& block, Structure& st) {
    // Gemmi❗✔️:   enum {
    // Gemmi❗✔️:     kId=0, kConnTypeId=1,
    // Gemmi❗✔️:     kAuthAsymId=2/*-3*/,  kLabelAsymId=4/*-5*/, kLabelCompId=6/*-7*/,
    // Gemmi❗✔️:     kLabelAtomId=8/*-9*/, kLabelAltId=10/*-11*/,
    // Gemmi❗✔️:     kAuthSeqId=12/*-13*/, kLabelSeqId=14/*-15*/, kInsCode=16/*-17*/,
    // Gemmi❗✔️:     kSym1=18, kSym2=19, kDistValue=20, kLinkId=21
    // Gemmi❗✔️:   };
    // Gemmi❗✔️:   // label_ identifiers are not sufficient for HOH:
    // Gemmi❗✔️:   // waters have null label_seq_id so we need auth_seq_id+icode.
    // Gemmi❗✔️:   // And since we need auth_seq_id, we also use auth_asym_id for consistency.
    // Gemmi❗✔️:   // Unless only label_*_id are available.
    // Gemmi❗✔️:   for (const auto row : block.find("_struct_conn.", {
    // Gemmi❗✔️:         "id", "conn_type_id",                                   // 0-1
    // Gemmi❗✔️:         "?ptnr1_auth_asym_id", "?ptnr2_auth_asym_id",           // 2-3
    // Gemmi❗✔️:         "?ptnr1_label_asym_id", "?ptnr2_label_asym_id",         // 4-5
    // Gemmi❗✔️:         "ptnr1_label_comp_id", "ptnr2_label_comp_id",           // 6-7
    // Gemmi❗✔️:         "ptnr1_label_atom_id", "ptnr2_label_atom_id",           // 8-9
    // Gemmi❗✔️:         "?pdbx_ptnr1_label_alt_id", "?pdbx_ptnr2_label_alt_id", // 10-11
    // Gemmi❗✔️:         "?ptnr1_auth_seq_id", "?ptnr2_auth_seq_id",             // 12-13
    // Gemmi❗✔️:         "?ptnr1_label_seq_id", "?ptnr2_label_seq_id",           // 14-15
    // Gemmi❗✔️:         "?pdbx_ptnr1_PDB_ins_code", "?pdbx_ptnr2_PDB_ins_code", // 16-17
    // Gemmi❗✔️:         "?ptnr1_symmetry", "?ptnr2_symmetry",                   // 18-19
    // Gemmi❗✔️:         "?pdbx_dist_value", "?ccp4_link_id"})) {                // 20-21
    // Gemmi❗✔️:     Connection c;
    // Gemmi❗✔️:     c.name = row.str(kId);
    // Gemmi❗✔️:     copy_string(row, kLinkId, c.link_id);
    // Gemmi❗✔️:     c.type = connection_type_from_string(row.str(kConnTypeId));
    // Gemmi❗✔️:     if (row.has2(kSym1) && row.has2(kSym2)) {
    // Gemmi❗✔️:       std::string s1 = row.str(kSym1);
    // Gemmi❗✔️:       std::string s2 = row.str(kSym2);
    // Gemmi❗✔️:       if (s1 == s2) {
    // Gemmi❗✔️:         c.asu = Asu::Same;
    // Gemmi❗✔️:       } else {
    // Gemmi❗✔️:         c.asu = Asu::Different;
    // Gemmi❗✔️:         size_t sep1 = s2.find('_');
    // Gemmi❗✔️:         size_t sep2 = s2.find('_');
    // Gemmi❗✔️:         if (sep1 != std::string::npos && sep1 + 4 == s1.size() &&
    // Gemmi❗✔️:             sep2 != std::string::npos && sep2 + 4 == s2.size()) {
    // Gemmi❗✔️:           if (s1[0] == '1' && s1[1] == '_')  // symop1 is usually 1_555
    // Gemmi❗✔️:             c.reported_sym[0] = no_sign_atoi(s2.c_str());
    // Gemmi❗✔️:           else
    // Gemmi❗✔️:             c.reported_sym[0] = 99;
    // Gemmi❗✔️:           for (size_t i = 1; i <= 3; ++i)
    // Gemmi❗✔️:             c.reported_sym[i] = s2[sep2 + i] - s1[sep1 + i];
    // Gemmi❗✔️:         }
    // Gemmi❗✔️:       }
    // Gemmi❗✔️:     }
    // Gemmi❗✔️:     copy_double(row, kDistValue, c.reported_distance);
    // Gemmi❗✔️:     for (int i = 0; i < 2; ++i) {
    // Gemmi❗✔️:       AtomAddress& a = (i == 0 ? c.partner1 : c.partner2);
    // Gemmi❗✔️:       if (row.has(kAuthAsymId+i) && row.has(kAuthSeqId+i)) {
    // Gemmi❗✔️:         a.chain_name = row.str(kAuthAsymId+i);
    // Gemmi❗✔️:         a.res_id = make_resid(row.str(kLabelCompId+i),
    // Gemmi❗✔️:                               row.str(kAuthSeqId+i), row.ptr_at(kInsCode+i));
    // Gemmi❗✔️:       } else if (row.has(kLabelAsymId+i) && row.has(kLabelSeqId+i)) {
    // Gemmi❗✔️:         set_part_of_address_from_label(a, st.first_model(),
    // Gemmi❗✔️:                                        row.str(kLabelAsymId+i),
    // Gemmi❗✔️:                                        row[kLabelSeqId+i]);
    // Gemmi❗✔️:         a.res_id.name = row.str(kLabelCompId+i);
    // Gemmi❗✔️:       } else {
    // Gemmi❗✔️:         fail("_struct_conn without either _auth_ or _label_ asym_id+seq_id");
    // Gemmi❗✔️:       }
    // Gemmi❗✔️:       a.atom_name = row.str(kLabelAtomId+i);
    // Gemmi❗✔️:       if (row.has2(kLabelAltId+i))
    // Gemmi❗✔️:         a.altloc = cif::as_char(row[kLabelAltId+i], '\0');
    // Gemmi❗✔️:     }
    // Gemmi❗✔️:     st.connections.emplace_back(c);
    // Gemmi❗✔️:   }
    // Gemmi❗✔️: }
    // Gemmi❗✔️: void copy_double(const cif::Table::Row& row, int n, double& dest) {
    // Gemmi❗✔️:   if (row.has2(n))
    // Gemmi❗✔️:     dest = cif::as_number(row[n]);
    // Gemmi❗✔️: }
    // Gemmi✔️✔️: void copy_string(const cif::Table::Row& row, int n, std::string& dest) {
    // Gemmi✔️✔️:   if (row.has2(n))
    // Gemmi✔️✔️:     dest = cif::as_string(row[n]);
    // Gemmi✔️✔️: }
    // Gemmi✔️✔️: inline Connection::Type connection_type_from_string(const std::string& t) {
    // Gemmi✔️✔️:   for (int i = 0; i != Connection::Unknown; ++i)
    // Gemmi✔️✔️:     if (connection_type_to_string(Connection::Type(i)) == t)
    // Gemmi✔️✔️:       return Connection::Type(i);
    // Gemmi✔️✔️:   return Connection::Unknown;
    // Gemmi✔️✔️: }
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
    // Behavior review: preserve row order, defaults, source `has` versus
    // `has2`, exact kind text, both-partner order, source address precedence,
    // and the pinned duplicate lookup of `s2` for `sep1`/`sep2`. Approved
    // chain/residue widths use structured errors. Signed `int` overflow is
    // undefined and conversion to source `short` outside range is
    // implementation-defined; Rust refuses both rather than claiming parity.
    // Empty `as_number`/decoded-empty `as_char` source access is undefined.
    // Complexity review: one source-order table pass and one owned output per
    // row; per-label lookup reuses the existing ordered chain/span scan, with
    // no duplicate hierarchy or temporary connection graph.
    const K_ID: usize = 0;
    const K_CONN_TYPE_ID: usize = 1;
    const K_AUTH_ASYM_ID: usize = 2;
    const K_LABEL_ASYM_ID: usize = 4;
    const K_LABEL_COMP_ID: usize = 6;
    const K_LABEL_ATOM_ID: usize = 8;
    const K_LABEL_ALT_ID: usize = 10;
    const K_AUTH_SEQ_ID: usize = 12;
    const K_LABEL_SEQ_ID: usize = 14;
    const K_INSERTION_CODE: usize = 16;
    const K_SYMMETRY_1: usize = 18;
    const K_SYMMETRY_2: usize = 19;
    const K_DISTANCE: usize = 20;
    const K_LINK_ID: usize = 21;
    const TAGS: [&str; 22] = [
        "id",
        "conn_type_id",
        "?ptnr1_auth_asym_id",
        "?ptnr2_auth_asym_id",
        "?ptnr1_label_asym_id",
        "?ptnr2_label_asym_id",
        "ptnr1_label_comp_id",
        "ptnr2_label_comp_id",
        "ptnr1_label_atom_id",
        "ptnr2_label_atom_id",
        "?pdbx_ptnr1_label_alt_id",
        "?pdbx_ptnr2_label_alt_id",
        "?ptnr1_auth_seq_id",
        "?ptnr2_auth_seq_id",
        "?ptnr1_label_seq_id",
        "?ptnr2_label_seq_id",
        "?pdbx_ptnr1_PDB_ins_code",
        "?pdbx_ptnr2_PDB_ins_code",
        "?ptnr1_symmetry",
        "?ptnr2_symmetry",
        "?pdbx_dist_value",
        "?ccp4_link_id",
    ];
    let table = block.find("_struct_conn.", &TAGS)?;
    let mut connections = Vec::with_capacity(table.len());
    let text_at =
        |row: CifRow<'_>, column| row.get(column).map(CifValue::decoded).unwrap_or_default();

    for row in table.iter() {
        let mut connection = BioConnection::default();
        connection.name = text_at(row, K_ID);
        if row.has_value(K_LINK_ID) {
            connection.link_id = text_at(row, K_LINK_ID);
        }

        connection.kind = match text_at(row, K_CONN_TYPE_ID).as_str() {
            "covale" => BioConnectionKind::Covale,
            "disulf" => BioConnectionKind::Disulf,
            "hydrog" => BioConnectionKind::Hydrog,
            "metalc" => BioConnectionKind::MetalC,
            _ => BioConnectionKind::Unknown,
        };

        if row.has_value(K_SYMMETRY_1) && row.has_value(K_SYMMETRY_2) {
            let symmetry_1 = text_at(row, K_SYMMETRY_1);
            let symmetry_2 = text_at(row, K_SYMMETRY_2);
            if symmetry_1 == symmetry_2 {
                connection.asu = BioAsu::Same;
            } else {
                connection.asu = BioAsu::Different;
                // The pinned source obtains BOTH positions from symmetry_2;
                // preserve the source rather than infer whether this is a typo.
                let separator_1 = symmetry_2.find('_');
                let separator_2 = symmetry_2.find('_');
                if let (Some(separator_1), Some(separator_2)) = (separator_1, separator_2)
                    && separator_1.checked_add(4) == Some(symmetry_1.len())
                    && separator_2.checked_add(4) == Some(symmetry_2.len())
                {
                    let symmetry_bytes_1 = symmetry_1.as_bytes();
                    let symmetry_bytes_2 = symmetry_2.as_bytes();
                    let symmetry_operator = if symmetry_bytes_1.first() == Some(&b'1')
                        && symmetry_bytes_1.get(1) == Some(&b'_')
                    {
                        let (value, _) = super::bio_pdb::gemmi_no_sign_atoi(symmetry_2.as_bytes())
                            .ok_or_else(|| {
                                MmcifConnectionError::SymmetryOperatorOutsideSourceDefinedRange {
                                    value: symmetry_2.clone(),
                                }
                            })?;
                        value
                    } else {
                        99
                    };
                    connection.reported_sym[0] =
                        i16::try_from(symmetry_operator).map_err(|_| {
                            MmcifConnectionError::SymmetryOperatorNotRepresentable {
                                value: symmetry_operator,
                            }
                        })?;
                    for index in 1..=3 {
                        let source_char_1 =
                            i8::from_ne_bytes([symmetry_bytes_1[separator_1 + index]]);
                        let source_char_2 =
                            i8::from_ne_bytes([symmetry_bytes_2[separator_2 + index]]);
                        connection.reported_sym[index] =
                            i16::from(source_char_2) - i16::from(source_char_1);
                    }
                }
            }
        }

        if row.has_value(K_DISTANCE) {
            let raw_distance = row.get(K_DISTANCE).map(CifValue::raw).unwrap_or("");
            connection.reported_distance = cif_as_f64(raw_distance, f64::NAN)?;
        }

        for partner_index in 0..2 {
            let auth_asym_column = K_AUTH_ASYM_ID + partner_index;
            let auth_sequence_column = K_AUTH_SEQ_ID + partner_index;
            let label_asym_column = K_LABEL_ASYM_ID + partner_index;
            let label_sequence_column = K_LABEL_SEQ_ID + partner_index;
            let component_column = K_LABEL_COMP_ID + partner_index;
            let atom_column = K_LABEL_ATOM_ID + partner_index;
            let insertion_column = K_INSERTION_CODE + partner_index;
            let altloc_column = K_LABEL_ALT_ID + partner_index;
            let partner_number = partner_index + 1;

            let (chain_name, sequence_id) =
                if row.has(auth_asym_column) && row.has(auth_sequence_column) {
                    let chain_text = text_at(row, auth_asym_column);
                    let chain_name =
                        PdbChainId::from_ascii(chain_text.as_bytes()).ok_or_else(|| {
                            MmcifConnectionError::ChainNameNotRepresentable {
                                partner: partner_number,
                                value: chain_text.clone(),
                            }
                        })?;
                    let insertion_code = row
                        .has(insertion_column)
                        .then(|| row.get(insertion_column).map(CifValue::raw).unwrap_or(""));
                    let sequence_id =
                        make_mmcif_seq_id(text_at(row, auth_sequence_column), insertion_code)
                            .map_err(MmcifConnectionError::SequenceId)?;
                    (chain_name, sequence_id)
                } else if row.has(label_asym_column) && row.has(label_sequence_column) {
                    let first_model = grouping
                        .models
                        .first()
                        .ok_or(MmcifConnectionError::NoStructuralModels)?;
                    let label_address = resolve_mmcif_label_address(
                        first_model,
                        &text_at(row, label_asym_column),
                        row.get(label_sequence_column)
                            .map(CifValue::raw)
                            .unwrap_or(""),
                    )
                    .map_err(MmcifConnectionError::LabelAddress)?;
                    match label_address {
                        Some((chain_name, Some(sequence_id))) => (chain_name, sequence_id),
                        Some((chain_name, None)) => (chain_name, PdbSeqId::new(i32::MIN, None)),
                        None => (
                            PdbChainId::from_ascii(b"")
                                .expect("empty chain identifier is representable"),
                            PdbSeqId::new(i32::MIN, None),
                        ),
                    }
                } else {
                    return Err(MmcifConnectionError::MissingAddressIdentifiers {
                        partner: partner_number,
                    });
                };

            let component_text = text_at(row, component_column);
            let residue_name =
                ResidueName::from_ascii(component_text.as_bytes()).ok_or_else(|| {
                    MmcifConnectionError::ResidueNameNotRepresentable {
                        partner: partner_number,
                        value: component_text.clone(),
                    }
                })?;
            let residue =
                make_mmcif_residue_address(residue_name, sequence_id).ok_or_else(|| {
                    MmcifConnectionError::ResidueAddressNotRepresentable {
                        partner: partner_number,
                        value: component_text,
                    }
                })?;
            let atom_name = text_at(row, atom_column);
            let altloc = if row.has_value(altloc_column) {
                let raw_altloc = row.get(altloc_column).map(CifValue::raw).unwrap_or("");
                cif_as_char(raw_altloc, '\0')? as u8
            } else {
                0
            };
            let address = AtomAddress::new(chain_name, residue, atom_name, Some(altloc));
            if partner_index == 0 {
                connection.partner1 = address;
            } else {
                connection.partner2 = address;
            }
        }

        connections.push(connection);
    }

    Ok(connections)
}

fn resolve_mmcif_label_address(
    model: &MmcifModelGrouping,
    label_asym_id: &str,
    raw_label_seq_id: &str,
) -> Result<Option<(PdbChainId, Option<PdbSeqId>)>, MmcifLabelAddressError> {
    // Gemmi❗✔️: void set_part_of_address_from_label(AtomAddress& a, const Model& model,
    // Gemmi❗✔️:                                     const std::string& label_asym,
    // Gemmi❗✔️:                                     const std::string& label_seq_id_raw) {
    // Gemmi❗✔️:   int seq = cif::as_int(label_seq_id_raw, SeqId::OptionalNum::None);
    // Gemmi❗✔️:   for (const Chain& chain : model.chains)
    // Gemmi❗✔️:     if (ConstResidueSpan sub = chain.get_subchain(label_asym)) {
    // Gemmi❗✔️:       a.chain_name = chain.name;
    // Gemmi❗✔️:       for (const Residue& res : sub)
    // Gemmi❗✔️:         if (res.label_seq == seq) {
    // Gemmi❗✔️:           a.res_id.seqid = res.seqid;
    // Gemmi❗✔️:           return;
    // Gemmi❗✔️:         }
    // Gemmi❗✔️:     }
    // Gemmi❗✔️: }
    // Gemmi❗✔️: inline int as_int(const std::string& str, int null) {
    // Gemmi❗✔️:   return is_null(str) ? null : as_int(str);
    // Gemmi❗✔️: }
    // Gemmi❗✔️: inline int as_int(const std::string& str) {
    // Gemmi❗✔️:   return string_to_int(str, true);
    // Gemmi❗✔️: }
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
    // Gemmi❗✔️:   for (; (length == 0 || i < length) && is_digit(p[i]); ++i) {
    // Gemmi❗✔️:     n = n * 10 - (p[i] - '0');
    // Gemmi❗✔️:     has_digits = true;
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   if (checked) {
    // Gemmi❗✔️:     while ((length == 0 || i < length) && is_space(p[i]))
    // Gemmi❗✔️:       ++i;
    // Gemmi❗✔️:     if (!has_digits || p[i] != '\0')
    // Gemmi❗✔️:       throw std::invalid_argument("not an integer: " +
    // Gemmi❗✔️:                                   std::string(p, length ? length : i+1));
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   return mult * n;
    // Gemmi❗✔️: }
    // Gemmi❗✔️: ConstResidueSpan get_subchain(const std::string& s) const {
    // Gemmi❗✔️:   return const_cast<Chain*>(this)->get_subchain(s);
    // Gemmi❗✔️: }
    // Gemmi❗✔️: ResidueSpan get_subchain(const std::string& s) {
    // Gemmi❗✔️:   return get_residue_span([&](const Residue& r) { return r.subchain == s; });
    // Gemmi❗✔️: }
    // Gemmi❗✔️: ResidueSpan whole() {
    // Gemmi❗✔️:   Residue* begin = residues.empty() ? nullptr : &residues[0];
    // Gemmi❗✔️:   return ResidueSpan(residues, begin, residues.size());
    // Gemmi❗✔️: }
    // Gemmi❗✔️: template<typename F> ResidueSpan get_residue_span(F&& func) {
    // Gemmi❗✔️:   return whole().subspan(func);
    // Gemmi❗✔️: }
    // Gemmi❗✔️: template<typename F, typename V=Item> Span<V> subspan(F&& func) {
    // Gemmi❗✔️:   iterator group_begin = std::find_if(this->begin(), this->end(), func);
    // Gemmi❗✔️:   iterator group_end = std::find_if_not(group_begin, this->end(), func);
    // Gemmi❗✔️:   return Span<V>(&*group_begin, group_end - group_begin);
    // Gemmi❗✔️: }
    // Behavior review: CIF nulls map to Gemmi's INT_MIN sentinel; checked
    // integer conversion preserves defined sign/whitespace/full-consumption
    // behavior. Chains and each first contiguous subchain are searched in
    // source order. A match returns that chain and the residue's complete
    // author sequence/insertion value; a matching subchain without a sequence
    // match returns its final chain assignment with no sequence replacement;
    // no matching subchain returns no address update. Fixed-width chain-name
    // failures remain structured. Signed C++ integer overflow is undefined
    // and is not claimed as parity.
    // Complexity review: for each chain, the first subchain search, contiguous
    // end search and label-sequence scan are linear, matching Gemmi's find-if,
    // find-if-not and span traversal. No span vector or per-residue allocation
    // is introduced; only the returned bounded chain value is copied.
    let sequence = if cif_is_null(raw_label_seq_id) {
        i32::MIN
    } else {
        cif_as_i32(raw_label_seq_id).map_err(MmcifLabelAddressError::LabelSequenceId)?
    };

    let mut last_matching_chain_name = None;
    for chain in &model.chains {
        let Some(start) = chain
            .residues
            .iter()
            .position(|residue| residue.source.subchain_id().unwrap_or_default() == label_asym_id)
        else {
            continue;
        };
        let after_start = &chain.residues[start..];
        let end = start
            + after_start
                .iter()
                .position(|residue| {
                    residue.source.subchain_id().unwrap_or_default() != label_asym_id
                })
                .unwrap_or(after_start.len());
        let subchain = &chain.residues[start..end];

        if let Some(residue) = subchain
            .iter()
            .find(|residue| residue.source.label_seq_id().unwrap_or(i32::MIN) == sequence)
        {
            let chain_name =
                PdbChainId::from_ascii(chain.source_name.as_bytes()).ok_or_else(|| {
                    MmcifLabelAddressError::ChainNameNotRepresentable {
                        value: chain.source_name.clone(),
                    }
                })?;
            let sequence_id = PdbSeqId::new(
                residue.address.sequence_number().unwrap_or(i32::MIN),
                residue.address.insertion_code(),
            );
            return Ok(Some((chain_name, Some(sequence_id))));
        }

        last_matching_chain_name = Some(chain.source_name.as_str());
    }

    last_matching_chain_name
        .map(|name| {
            PdbChainId::from_ascii(name.as_bytes())
                .map(|chain_name| (chain_name, None))
                .ok_or_else(|| MmcifLabelAddressError::ChainNameNotRepresentable {
                    value: name.to_owned(),
                })
        })
        .transpose()
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum MmcifCisPepError {
    Cif(CifReadError),
    ChainNameNotRepresentable {
        partner: &'static str,
        value: String,
    },
    ResidueNameNotRepresentable {
        partner: &'static str,
        value: String,
    },
    ResidueAddressNotRepresentable {
        partner: &'static str,
    },
    SequenceId {
        partner: &'static str,
        error: MmcifSequenceIdError,
    },
}

impl From<CifReadError> for MmcifCisPepError {
    fn from(error: CifReadError) -> Self {
        Self::Cif(error)
    }
}

impl fmt::Display for MmcifCisPepError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Cif(error) => fmt::Display::fmt(error, f),
            Self::ChainNameNotRepresentable { partner, value } => write!(
                f,
                "{partner} cis-peptide chain name exceeds the approved BIO representation: {value:?}"
            ),
            Self::ResidueNameNotRepresentable { partner, value } => write!(
                f,
                "{partner} cis-peptide residue name exceeds the approved BIO representation: {value:?}"
            ),
            Self::ResidueAddressNotRepresentable { partner } => {
                write!(
                    f,
                    "{partner} cis-peptide residue address is not representable"
                )
            }
            Self::SequenceId { partner, error } => {
                write!(
                    f,
                    "invalid {partner} cis-peptide sequence identifier: {error}"
                )
            }
        }
    }
}

impl std::error::Error for MmcifCisPepError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Cif(error) => Some(error),
            Self::SequenceId { error, .. } => Some(error),
            Self::ChainNameNotRepresentable { .. }
            | Self::ResidueNameNotRepresentable { .. }
            | Self::ResidueAddressNotRepresentable { .. } => None,
        }
    }
}

fn parse_mmcif_cis_peptides(block: &CifBlock) -> Result<Vec<BioCisPep>, MmcifCisPepError> {
    // Gemmi❗✔️: void read_prot_cis(cif::Block& block, Structure& st) {
    // Gemmi❗✔️:   enum {
    // Gemmi❗✔️:     kModelNum=0,
    // Gemmi❗✔️:     kAuthAsymId=1, kAuthSeqId=2, kInsCode=3, kLabelCompId=4, kAuthCompId=5,
    // Gemmi❗✔️:     kAuthAsymId2=6, kAuthSeqId2=7, kInsCode2=8, kLabelCompId2=9, kAuthCompId2=10,
    // Gemmi❗✔️:     kAltId=11, kOmegaAngle=12
    // Gemmi❗✔️:   };
    // Gemmi❗✔️:   // We could use label_seq_id etc and call set_part_of_address_from_label(),
    // Gemmi❗✔️:   // but for now let's assume that auth_seq_id etc are there.
    // Gemmi❗✔️:   for (auto row : block.find("_struct_mon_prot_cis.",
    // Gemmi❗✔️:                              {"pdbx_PDB_model_num",                  // 0
    // Gemmi❗✔️:                               "auth_asym_id",                        // 1
    // Gemmi❗✔️:                               "auth_seq_id", "?pdbx_PDB_ins_code",   // 2-3
    // Gemmi❗✔️:                               "?label_comp_id", "?auth_comp_id",     // 4-5
    // Gemmi❗✔️:                               "?pdbx_auth_asym_id_2",                // 6
    // Gemmi❗✔️:                               "?pdbx_auth_seq_id_2", "?pdbx_PDB_ins_code_2",   // 7-8
    // Gemmi❗✔️:                               "?pdbx_label_comp_id_2", "?pdbx_auth_comp_id_2", // 9-10
    // Gemmi❗✔️:                               "?label_alt_id", "?pdbx_omega_angle"})) {        // 11-12
    // Gemmi❗✔️:     CisPep cispep;
    // Gemmi❗✔️:     cispep.model_num = cif::as_int(row[kModelNum], 0);
    // Gemmi❗✔️:     cispep.partner_c.chain_name = row.str(kAuthAsymId);
    // Gemmi❗✔️:     cispep.partner_c.res_id.seqid = make_seqid(row.str(kAuthSeqId), row.ptr_at(kInsCode));
    // Gemmi❗✔️:     cispep.partner_c.res_id.name = cif::as_string(row.one_of(kAuthCompId, kLabelCompId));
    // Gemmi❗✔️:     if (row.has(kAuthAsymId2))
    // Gemmi❗✔️:       cispep.partner_n.chain_name = row.str(kAuthAsymId2);
    // Gemmi❗✔️:     if (row.has(kAuthSeqId2))
    // Gemmi❗✔️:       cispep.partner_n.res_id.seqid = make_seqid(row.str(kAuthSeqId2), row.ptr_at(kInsCode2));
    // Gemmi❗✔️:     cispep.partner_n.res_id.name = cif::as_string(row.one_of(kAuthCompId2, kLabelCompId2));
    // Gemmi❗✔️:     if (row.has(kAltId))
    // Gemmi❗✔️:       cispep.only_altloc = cif::as_char(row[kAltId], '\0');
    // Gemmi❗✔️:     if (row.has(kOmegaAngle))
    // Gemmi❗✔️:       cispep.reported_angle = cif::as_number(row[kOmegaAngle]);
    // Gemmi❗✔️:     st.cispeps.push_back(cispep);
    // Gemmi❗✔️:   }
    // Gemmi❗✔️: }
    // Behavior review: preserve source row order, required first-partner
    // fields, optional second-partner defaults, one_of auth/label name
    // precedence, and has-based altloc/angle handling. The canonical
    // make_mmcif_seq_id owner carries insertion-code and sentinel behavior;
    // approved bounded BIO names produce structured representation errors.
    // Complexity review: one ordered table pass and one output per row;
    // per-row conversions and bounded name checks do not scan other rows.
    const TAGS: [&str; 13] = [
        "pdbx_PDB_model_num",
        "auth_asym_id",
        "auth_seq_id",
        "?pdbx_PDB_ins_code",
        "?label_comp_id",
        "?auth_comp_id",
        "?pdbx_auth_asym_id_2",
        "?pdbx_auth_seq_id_2",
        "?pdbx_PDB_ins_code_2",
        "?pdbx_label_comp_id_2",
        "?pdbx_auth_comp_id_2",
        "?label_alt_id",
        "?pdbx_omega_angle",
    ];
    let table = block.find("_struct_mon_prot_cis.", &TAGS)?;
    let mut cis_peptides = Vec::with_capacity(table.len());
    for row in table.iter() {
        let mut cis_peptide = BioCisPep::default();
        let raw_model = row.get(0).map(CifValue::raw).unwrap_or("");
        cis_peptide.model_num = if cif_is_null(raw_model) {
            0
        } else {
            cif_as_i32(raw_model)?
        };

        let chain_c_text = row.decoded(1).unwrap_or_default();
        let chain_c = PdbChainId::from_ascii(chain_c_text.as_bytes()).ok_or_else(|| {
            MmcifCisPepError::ChainNameNotRepresentable {
                partner: "C",
                value: chain_c_text,
            }
        })?;
        let sequence_c = make_mmcif_seq_id(
            row.decoded(2).unwrap_or_default(),
            row.has(3)
                .then(|| row.get(3).map(CifValue::raw).unwrap_or("")),
        )
        .map_err(|error| MmcifCisPepError::SequenceId {
            partner: "C",
            error,
        })?;
        let name_c_text = cif_as_string(row.one_of(5, 4)?);
        let name_c = ResidueName::from_ascii(name_c_text.as_bytes()).ok_or_else(|| {
            MmcifCisPepError::ResidueNameNotRepresentable {
                partner: "C",
                value: name_c_text,
            }
        })?;
        let residue_c = make_mmcif_residue_address(name_c, sequence_c)
            .ok_or(MmcifCisPepError::ResidueAddressNotRepresentable { partner: "C" })?;
        cis_peptide.partner_c = AtomAddress::new(chain_c, residue_c, String::new(), None);

        let chain_n_text = if row.has(6) {
            row.decoded(6).unwrap_or_default()
        } else {
            String::new()
        };
        let chain_n = PdbChainId::from_ascii(chain_n_text.as_bytes()).ok_or_else(|| {
            MmcifCisPepError::ChainNameNotRepresentable {
                partner: "N",
                value: chain_n_text,
            }
        })?;
        let sequence_n = if row.has(7) {
            make_mmcif_seq_id(
                row.decoded(7).unwrap_or_default(),
                row.has(8)
                    .then(|| row.get(8).map(CifValue::raw).unwrap_or("")),
            )
            .map_err(|error| MmcifCisPepError::SequenceId {
                partner: "N",
                error,
            })?
        } else {
            PdbSeqId::new(i32::MIN, None)
        };
        let name_n_text = cif_as_string(row.one_of(10, 9)?);
        let name_n = ResidueName::from_ascii(name_n_text.as_bytes()).ok_or_else(|| {
            MmcifCisPepError::ResidueNameNotRepresentable {
                partner: "N",
                value: name_n_text,
            }
        })?;
        let residue_n = make_mmcif_residue_address(name_n, sequence_n)
            .ok_or(MmcifCisPepError::ResidueAddressNotRepresentable { partner: "N" })?;
        cis_peptide.partner_n = AtomAddress::new(chain_n, residue_n, String::new(), None);

        if row.has(11) {
            cis_peptide.only_altloc =
                cif_as_char(row.get(11).map(CifValue::raw).unwrap_or(""), '\0')? as u8;
        }
        if row.has(12) {
            cis_peptide.reported_angle =
                cif_as_f64(row.get(12).map(CifValue::raw).unwrap_or(""), f64::NAN)?;
        }
        cis_peptides.push(cis_peptide);
    }
    Ok(cis_peptides)
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum MmcifModResError {
    Cif(CifReadError),
    ChainNameNotRepresentable { value: String },
    ResidueNameNotRepresentable { value: String },
    ResidueAddressNotRepresentable,
    SequenceId(MmcifSequenceIdError),
}

impl From<CifReadError> for MmcifModResError {
    fn from(error: CifReadError) -> Self {
        Self::Cif(error)
    }
}

impl fmt::Display for MmcifModResError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Cif(error) => fmt::Display::fmt(error, f),
            Self::ChainNameNotRepresentable { value } => {
                write!(
                    f,
                    "modified-residue chain name exceeds the approved BIO representation: {value:?}"
                )
            }
            Self::ResidueNameNotRepresentable { value } => {
                write!(
                    f,
                    "modified-residue name exceeds the approved BIO representation: {value:?}"
                )
            }
            Self::ResidueAddressNotRepresentable => {
                write!(f, "modified-residue address is not representable")
            }
            Self::SequenceId(error) => {
                write!(f, "invalid modified-residue sequence identifier: {error}")
            }
        }
    }
}

impl std::error::Error for MmcifModResError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Cif(error) => Some(error),
            Self::SequenceId(error) => Some(error),
            Self::ChainNameNotRepresentable { .. }
            | Self::ResidueNameNotRepresentable { .. }
            | Self::ResidueAddressNotRepresentable => None,
        }
    }
}

fn parse_mmcif_modified_residues(block: &CifBlock) -> Result<Vec<BioModRes>, MmcifModResError> {
    // Gemmi❗✔️: // MODRES equivalent
    // Gemmi❗✔️: void read_struct_mod_residue(cif::Block& block, Structure& st) {
    // Gemmi❗✔️:   // Here auth_asym_id etc are mandatory and label_asym_id etc optional.
    // Gemmi❗✔️:   for (auto row : block.find("_pdbx_struct_mod_residue.",
    // Gemmi❗✔️:                              {"auth_asym_id",  // 0
    // Gemmi❗✔️:                               "auth_seq_id", "?PDB_ins_code",  // 1-2
    // Gemmi❗✔️:                               "?auth_comp_id", "?label_comp_id",  // 3-4
    // Gemmi❗✔️:                               "?parent_comp_id", "?details",  // 5-6
    // Gemmi❗✔️:                               "?ccp4_mod_id"})) {  // 7
    // Gemmi❗✔️:     ModRes modres;
    // Gemmi❗✔️:     modres.chain_name = row.str(0);
    // Gemmi❗✔️:     modres.res_id.seqid = make_seqid(row.str(1), row.ptr_at(2));
    // Gemmi❗✔️:     modres.res_id.name = row.one_of(3, 4);
    // Gemmi❗✔️:     if (row.has(5))
    // Gemmi❗✔️:       modres.parent_comp_id = row.str(5);
    // Gemmi❗✔️:     if (row.has(6))
    // Gemmi❗✔️:       modres.details = row.str(6);
    // Gemmi❗✔️:     if (row.has(7))
    // Gemmi❗✔️:       modres.mod_id = row.str(7);
    // Gemmi❗✔️:     st.mod_residues.push_back(modres);
    // Gemmi❗✔️:   }
    // Gemmi❗✔️: }
    // Behavior review: required auth chain/sequence gate the table; name
    // selection keeps one_of's raw token rather than decoding it, while the
    // optional parent/details/mod-id columns use source str decoding. Rows
    // append in input order. BIO's approved bounded names fail explicitly.
    // Complexity review: one ordered table pass, one output allocation per
    // row, and fixed-field conversions. No lookup over other residue rows.
    const TAGS: [&str; 8] = [
        "auth_asym_id",
        "auth_seq_id",
        "?PDB_ins_code",
        "?auth_comp_id",
        "?label_comp_id",
        "?parent_comp_id",
        "?details",
        "?ccp4_mod_id",
    ];
    let table = block.find("_pdbx_struct_mod_residue.", &TAGS)?;
    let mut modified_residues = Vec::with_capacity(table.len());
    for row in table.iter() {
        let mut modified = BioModRes::default();
        let chain_text = row.decoded(0).unwrap_or_default();
        modified.chain_name = PdbChainId::from_ascii(chain_text.as_bytes())
            .ok_or(MmcifModResError::ChainNameNotRepresentable { value: chain_text })?;
        let sequence = make_mmcif_seq_id(
            row.decoded(1).unwrap_or_default(),
            row.has(2)
                .then(|| row.get(2).map(CifValue::raw).unwrap_or("")),
        )
        .map_err(MmcifModResError::SequenceId)?;
        let name_text = row.one_of(3, 4)?.to_owned();
        let name = ResidueName::from_ascii(name_text.as_bytes())
            .ok_or(MmcifModResError::ResidueNameNotRepresentable { value: name_text })?;
        modified.res_id = make_mmcif_residue_address(name, sequence)
            .ok_or(MmcifModResError::ResidueAddressNotRepresentable)?;
        if row.has(5) {
            modified.parent_comp_id = row.decoded(5).unwrap_or_default();
        }
        if row.has(6) {
            modified.details = row.decoded(6).unwrap_or_default();
        }
        if row.has(7) {
            modified.mod_id = row.decoded(7).unwrap_or_default();
        }
        modified_residues.push(modified);
    }
    Ok(modified_residues)
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum MmcifOperationExprError {
    SourceUndefinedIntegerOverflow,
    SourceUndefinedRangeIncrement,
}

impl fmt::Display for MmcifOperationExprError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::SourceUndefinedIntegerOverflow => {
                write!(
                    f,
                    "assembly operation range integer overflows Gemmi's signed int"
                )
            }
            Self::SourceUndefinedRangeIncrement => {
                write!(
                    f,
                    "assembly operation range reaches Gemmi's undefined signed increment"
                )
            }
        }
    }
}

impl std::error::Error for MmcifOperationExprError {}

fn parse_mmcif_operation_expr(expr: &str) -> Result<Vec<String>, MmcifOperationExprError> {
    // Gemmi❗✔️: // Operation expression is an item type used for *.oper_expression.
    // Gemmi❗✔️: // Here, to keep it simple, we ignore products such as "(2)(3)".
    // Gemmi❗✔️: // We parse "3", "1,3,5", "one,two", "(3)", "(a)", "(1-60)", "(2,3-8,XY)", etc
    // Gemmi❗✔️: std::vector<std::string> parse_operation_expr(const std::string& expr) {
    // Gemmi❗✔️:   std::vector<std::string> result;
    // Gemmi❗✔️:   std::size_t start = 0;
    // Gemmi❗✔️:   std::size_t close_br = std::string::npos;
    // Gemmi❗✔️:   if (expr[0] == '(') {
    // Gemmi❗✔️:     start = 1;
    // Gemmi❗✔️:     close_br = expr.find(')');
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   for (;;) {
    // Gemmi❗✔️:     std::size_t sep = std::min(expr.find(',', start), close_br);
    // Gemmi❗✔️:     std::size_t minus = expr.find('-', start);
    // Gemmi❗✔️:     if (minus < sep) {
    // Gemmi❗✔️:       int n_min = no_sign_atoi(expr.c_str() + start);
    // Gemmi❗✔️:       int n_max = no_sign_atoi(expr.c_str() + minus + 1);
    // Gemmi❗✔️:       for (int n = n_min; n <= n_max; ++n)
    // Gemmi❗✔️:         result.push_back(std::to_string(n));
    // Gemmi❗✔️:     } else {
    // Gemmi❗✔️:       result.emplace_back(expr, start, sep - start);
    // Gemmi❗✔️:     }
    // Gemmi❗✔️:     if (sep == close_br)
    // Gemmi❗✔️:       break;
    // Gemmi❗✔️:     start = sep + 1;
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   return result;
    // Gemmi❗✔️: }
    // Behavior review: only the first parenthesized group is parsed, even
    // when a product follows; unmatched/empty tokens and no-conversion range
    // endpoints keep source results. The shared no_sign_atoi owner preserves
    // C-locale whitespace/digit-prefix behavior. Signed overflow and a range
    // increment past INT_MAX are source-undefined, represented by typed errors.
    // Complexity review: each separator search scans the remaining expression
    // suffix as the source does; output construction is linear in emitted text.
    // No full-expression tokenization or Cartesian-product allocation occurs.
    let mut result = Vec::new();
    let mut start = 0;
    let mut close_br = usize::MAX;
    if expr.as_bytes().first() == Some(&b'(') {
        start = 1;
        close_br = expr.find(')').unwrap_or(usize::MAX);
    }
    loop {
        let sep = expr[start..]
            .find(',')
            .map(|relative| start + relative)
            .unwrap_or(usize::MAX)
            .min(close_br);
        let minus = expr[start..]
            .find('-')
            .map(|relative| start + relative)
            .unwrap_or(usize::MAX);
        if minus < sep {
            let n_min = super::bio_pdb::gemmi_no_sign_atoi(&expr.as_bytes()[start..])
                .ok_or(MmcifOperationExprError::SourceUndefinedIntegerOverflow)?
                .0;
            let n_max = super::bio_pdb::gemmi_no_sign_atoi(&expr.as_bytes()[minus + 1..])
                .ok_or(MmcifOperationExprError::SourceUndefinedIntegerOverflow)?
                .0;
            if n_min <= n_max && n_max == i32::MAX {
                return Err(MmcifOperationExprError::SourceUndefinedRangeIncrement);
            }
            for number in n_min..=n_max {
                result.push(number.to_string());
            }
        } else {
            result.push(expr[start..sep.min(expr.len())].to_owned());
        }
        if sep == close_br {
            break;
        }
        start = sep + 1;
    }
    Ok(result)
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum MmcifAssemblyError {
    Cif(CifReadError),
    Transform(MmcifCrystalInfoError),
    OperationExpr(MmcifOperationExprError),
}

impl From<CifReadError> for MmcifAssemblyError {
    fn from(error: CifReadError) -> Self {
        Self::Cif(error)
    }
}

impl fmt::Display for MmcifAssemblyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Cif(error) => fmt::Display::fmt(error, f),
            Self::Transform(error) => fmt::Display::fmt(error, f),
            Self::OperationExpr(error) => fmt::Display::fmt(error, f),
        }
    }
}

impl std::error::Error for MmcifAssemblyError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Cif(error) => Some(error),
            Self::Transform(error) => Some(error),
            Self::OperationExpr(error) => Some(error),
        }
    }
}

fn parse_mmcif_assemblies(block: &CifBlock) -> Result<Vec<BioAssembly>, MmcifAssemblyError> {
    // Gemmi❗✔️: std::vector<Assembly> read_assemblies(cif::Block& block) {
    // Gemmi❗✔️:   std::vector<Assembly> assemblies;
    // Gemmi❗✔️:   cif::Table prop_tab = block.find("_pdbx_struct_assembly_prop.",
    // Gemmi❗✔️:                                    {"biol_id", "type", "value"});
    // Gemmi❗✔️:   cif::Table gen_tab = block.find("_pdbx_struct_assembly_gen.",
    // Gemmi❗✔️:                           {"assembly_id", "oper_expression", "asym_id_list"});
    // Gemmi❗✔️:   std::vector<Assembly::Operator> oper_list;
    // Gemmi❗✔️:   std::vector<std::string> oper_list_tags = transform_tags("matrix", "vector");
    // Gemmi❗✔️:   oper_list_tags.emplace_back("id");  // 12
    // Gemmi❗✔️:   oper_list_tags.emplace_back("type");  // 13
    // Gemmi❗✔️:   for (const auto row : block.find("_pdbx_struct_oper_list.", oper_list_tags)) {
    // Gemmi❗✔️:     oper_list.emplace_back();
    // Gemmi❗✔️:     oper_list.back().name = row.str(12);
    // Gemmi❗✔️:     oper_list.back().type = row.str(13);
    // Gemmi❗✔️:     oper_list.back().transform = get_transform_matrix(row);
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   for (const auto row : block.find("_pdbx_struct_assembly.", {
    // Gemmi❗✔️:         "id", "details", "method_details",
    // Gemmi❗✔️:         "oligomeric_details", "oligomeric_count"})) {
    // Gemmi❗✔️:     assemblies.emplace_back(row.str(0));
    // Gemmi❗✔️:     Assembly& a = assemblies.back();
    // Gemmi❗✔️:     std::string detail = row.str(1);
    // Gemmi❗✔️:     if (detail == "author_and_software_defined_assembly")
    // Gemmi❗✔️:       a.author_determined = a.software_determined = true;
    // Gemmi❗✔️:     else if (detail == "author_defined_assembly")
    // Gemmi❗✔️:       a.author_determined = true;
    // Gemmi❗✔️:     else if (detail == "software_defined_assembly")
    // Gemmi❗✔️:       a.software_determined = true;
    // Gemmi❗✔️:     else if (detail == "complete icosahedral assembly")
    // Gemmi❗✔️:       a.special_kind = Assembly::SpecialKind::CompleteIcosahedral;
    // Gemmi❗✔️:     else if (detail == "representative helical assembly")
    // Gemmi❗✔️:       a.special_kind = Assembly::SpecialKind::RepresentativeHelical;
    // Gemmi❗✔️:     else if (detail == "complete point assembly")
    // Gemmi❗✔️:       a.special_kind = Assembly::SpecialKind::CompletePoint;
    // Gemmi❗✔️:     if (!a.author_determined && !a.software_determined &&
    // Gemmi❗✔️:         a.special_kind == Assembly::SpecialKind::NA && !detail.empty()) {
    // Gemmi❗✔️:       assemblies.pop_back();
    // Gemmi❗✔️:       continue;
    // Gemmi❗✔️:     }
    // Gemmi❗✔️:     if (a.software_determined && !cif::is_null(row[2]))
    // Gemmi❗✔️:       a.software_name = row.str(2);  // method_details
    // Gemmi❗✔️:     a.oligomeric_details = row.str(3);
    // Gemmi❗✔️:     a.oligomeric_count = cif::as_int(row[4], 0);
    // Gemmi❗✔️:     for (const auto row_p : prop_tab)
    // Gemmi❗✔️:       if (row_p.str(0) == a.name) {
    // Gemmi❗✔️:         std::string type = row_p.str(1);
    // Gemmi❗✔️:         double value = cif::as_number(row_p[2]);
    // Gemmi❗✔️:         if (type == "ABSA (A^2)")
    // Gemmi❗✔️:           a.absa = value;
    // Gemmi❗✔️:         else if (type == "SSA (A^2)")
    // Gemmi❗✔️:           a.ssa = value;
    // Gemmi❗✔️:         else if (type == "MORE")
    // Gemmi❗✔️:           a.more = value;
    // Gemmi❗✔️:       }
    // Gemmi❗✔️:     for (const auto row_g : gen_tab)
    // Gemmi❗✔️:       if (row_g.str(0) == a.name) {
    // Gemmi❗✔️:         a.generators.emplace_back();
    // Gemmi❗✔️:         Assembly::Gen& gen = a.generators.back();
    // Gemmi❗✔️:         split_str_into(row_g.str(2), ',', gen.subchains);
    // Gemmi❗✔️:         for (const std::string& name : parse_operation_expr(row_g.str(1)))
    // Gemmi❗✔️:           if (const Assembly::Operator* oper = impl::find_or_null(oper_list, name))
    // Gemmi❗✔️:             gen.operators.push_back(*oper);
    // Gemmi❗✔️:       }
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   return assemblies;
    // Gemmi❗✔️: }
    // Behavior review: operators are loaded before assemblies; unknown detail
    // rows are skipped before scalar/prop/gen conversions. Properties and
    // generators retain source row order and duplicate updates; missing
    // operator references are skipped, and duplicate IDs resolve to the first
    // operator. Expression products are deliberately not composed by source.
    // Split subchain fields retain empty tokens. Numeric boundaries reuse the
    // canonical CIF/transform owners and typed source-undefined expression
    // errors. No BioStructureData validation is implied by this detached parse.
    // Complexity review: fixed 14-column operator conversion followed by
    // source-shaped assembly×property and assembly×generator scans, plus a
    // linear first-match operator scan per expanded name. Operators are cloned
    // when referenced, matching source vector copies; no index changes ties.
    const PROP_TAGS: [&str; 3] = ["biol_id", "type", "value"];
    const GEN_TAGS: [&str; 3] = ["assembly_id", "oper_expression", "asym_id_list"];
    const ASSEMBLY_TAGS: [&str; 5] = [
        "id",
        "details",
        "method_details",
        "oligomeric_details",
        "oligomeric_count",
    ];
    let properties = block.find("_pdbx_struct_assembly_prop.", &PROP_TAGS)?;
    let generators = block.find("_pdbx_struct_assembly_gen.", &GEN_TAGS)?;
    let mut operator_tags = mmcif_transform_tags("matrix", "vector");
    operator_tags.push("id".to_owned());
    operator_tags.push("type".to_owned());
    let operator_refs: Vec<&str> = operator_tags.iter().map(String::as_str).collect();
    let operator_table = block.find("_pdbx_struct_oper_list.", &operator_refs)?;
    let mut operators = Vec::with_capacity(operator_table.len());
    for row in operator_table.iter() {
        operators.push(BioAssemblyOperator::new(
            Some(row.decoded(12).unwrap_or_default()),
            Some(row.decoded(13).unwrap_or_default()),
            mmcif_transform_from_row(row).map_err(MmcifAssemblyError::Transform)?,
        ));
    }

    let assembly_table = block.find("_pdbx_struct_assembly.", &ASSEMBLY_TAGS)?;
    let mut assemblies = Vec::with_capacity(assembly_table.len());
    for row in assembly_table.iter() {
        let name = row.decoded(0).unwrap_or_default();
        let detail = row.decoded(1).unwrap_or_default();
        let mut author_determined = false;
        let mut software_determined = false;
        let special_kind = match detail.as_str() {
            "author_and_software_defined_assembly" => {
                author_determined = true;
                software_determined = true;
                BioAssemblySpecialKind::NotApplicable
            }
            "author_defined_assembly" => {
                author_determined = true;
                BioAssemblySpecialKind::NotApplicable
            }
            "software_defined_assembly" => {
                software_determined = true;
                BioAssemblySpecialKind::NotApplicable
            }
            "complete icosahedral assembly" => BioAssemblySpecialKind::CompleteIcosahedral,
            "representative helical assembly" => BioAssemblySpecialKind::RepresentativeHelical,
            "complete point assembly" => BioAssemblySpecialKind::CompletePoint,
            "" => BioAssemblySpecialKind::NotApplicable,
            _ => continue,
        };
        let software_name =
            if software_determined && !cif_is_null(row.get(2).map(CifValue::raw).unwrap_or("")) {
                row.decoded(2).unwrap_or_default()
            } else {
                String::new()
            };
        let oligomeric_details = row.decoded(3).unwrap_or_default();
        let raw_count = row.get(4).map(CifValue::raw).unwrap_or("");
        let oligomeric_count = if cif_is_null(raw_count) {
            0
        } else {
            cif_as_i32(raw_count)?
        };
        let mut absa = f64::NAN;
        let mut ssa = f64::NAN;
        let mut more = f64::NAN;
        for property in properties.iter() {
            if property.decoded(0).as_deref() != Some(name.as_str()) {
                continue;
            }
            let kind = property.decoded(1).unwrap_or_default();
            let value = cif_as_f64(property.get(2).map(CifValue::raw).unwrap_or(""), f64::NAN)?;
            match kind.as_str() {
                "ABSA (A^2)" => absa = value,
                "SSA (A^2)" => ssa = value,
                "MORE" => more = value,
                _ => {}
            }
        }
        let mut assembly_generators = Vec::new();
        for generator in generators.iter() {
            if generator.decoded(0).as_deref() != Some(name.as_str()) {
                continue;
            }
            let subchains = generator
                .decoded(2)
                .unwrap_or_default()
                .split(',')
                .map(str::to_owned)
                .collect();
            let mut selected_operators = Vec::new();
            for requested in parse_mmcif_operation_expr(&generator.decoded(1).unwrap_or_default())
                .map_err(MmcifAssemblyError::OperationExpr)?
            {
                if let Some(operator) = operators
                    .iter()
                    .find(|operator| operator.name.as_deref() == Some(requested.as_str()))
                {
                    selected_operators.push(operator.clone());
                }
            }
            assembly_generators.push(BioAssemblyGenerator::new(
                Vec::new(),
                subchains,
                selected_operators,
            ));
        }
        assemblies.push(BioAssembly::new(
            name,
            author_determined,
            software_determined,
            special_kind,
            oligomeric_count,
            oligomeric_details,
            software_name,
            absa,
            ssa,
            more,
            assembly_generators,
        ));
    }
    Ok(assemblies)
}

fn parse_mmcif_sifts_unp(
    block: &CifBlock,
    entities: &mut [MmcifEntityBuildState],
    models: &mut [MmcifModelGrouping],
) -> Result<(), MmcifSiftsUnpError> {
    // Gemmi❗✔️: void read_sifts_unp(cif::Block& block, Structure& st) {
    // Gemmi❗✔️:   enum { kEntityId=0, kAsymId, kSeqIdOrdinal, kSeqId, kObserved,
    // Gemmi❗✔️:          kUnpRes, kUnpNum, kUnpAcc };
    // Gemmi❗✔️:   cif::Table table = block.find("_pdbx_sifts_xref_db.", {
    // Gemmi❗✔️:       "entity_id", "asym_id", "seq_id_ordinal", "seq_id", "observed",
    // Gemmi❗✔️:       "unp_res", "unp_num", "unp_acc"});
    // Gemmi❗✔️:   if (!table.ok())
    // Gemmi❗✔️:     return;
    // Gemmi❗✔️:   for (Model& model : st.models) {
    // Gemmi❗✔️:     Entity* ent = nullptr;
    // Gemmi❗✔️:     ResidueSpan polymer;
    // Gemmi❗✔️:     Residue* res = nullptr;
    // Gemmi❗✔️:     std::string unp_acc;
    // Gemmi❗✔️:     SiftsUnpResidue unp;
    // Gemmi❗✔️:     for (const auto row : table) {
    // Gemmi❗✔️:       if (row[kSeqIdOrdinal] != "1" || row[kObserved][0] != 'y')
    // Gemmi❗✔️:         continue;
    // Gemmi❗✔️:       if (cif::is_null(row[kUnpAcc]) || cif::is_null(row[kUnpNum]))
    // Gemmi❗✔️:         continue;
    // Gemmi❗✔️:       bool update_acc_index = false;
    // Gemmi❗✔️:       if (!ent || row[kEntityId] != ent->name) {
    // Gemmi❗✔️:         ent = st.get_entity(row[kEntityId]);
    // Gemmi❗✔️:         if (!ent)
    // Gemmi❗✔️:           fail("_pdbx_sifts_xref_db: entity_id not found: " + row[kEntityId]);
    // Gemmi❗✔️:         update_acc_index = true;
    // Gemmi❗✔️:       }
    // Gemmi❗✔️:       if (row[kUnpAcc] != unp_acc) {
    // Gemmi❗✔️:         unp_acc = row.str(kUnpAcc);
    // Gemmi❗✔️:         update_acc_index = true;
    // Gemmi❗✔️:       }
    // Gemmi❗✔️:       if (update_acc_index) {
    // Gemmi❗✔️:         auto& vec = ent->sifts_unp_acc;
    // Gemmi❗✔️:         auto it = std::find(vec.begin(), vec.end(), unp_acc);
    // Gemmi❗✔️:         unp.acc_index = std::uint8_t(it - vec.begin());
    // Gemmi❗✔️:         if (it == vec.end())
    // Gemmi❗✔️:           vec.push_back(unp_acc);
    // Gemmi❗✔️:       }
    // Gemmi❗✔️:       if (!polymer || row[kAsymId] != polymer.front().subchain) {
    // Gemmi❗✔️:         polymer = model.get_subchain(row[kAsymId]);
    // Gemmi❗✔️:         if (!polymer)
    // Gemmi❗✔️:           fail("_pdbx_sifts_xref_db: asym_id not found: " + row[kAsymId]);
    // Gemmi❗✔️:         res = polymer.begin();
    // Gemmi❗✔️:       } else if (res == polymer.end()) {
    // Gemmi❗✔️:         res = polymer.begin();
    // Gemmi❗✔️:       }
    // Gemmi❗✔️:       int label_seq = cif::as_int(row[kSeqId]);
    // Gemmi❗✔️:       if (res->label_seq != label_seq) {
    // Gemmi❗✔️:         res = polymer.begin();
    // Gemmi❗✔️:         while (res->label_seq != label_seq) {
    // Gemmi❗✔️:           ++res;
    // Gemmi❗✔️:           if (res == polymer.end())
    // Gemmi❗✔️:             fail("_pdbx_sifts_xref_db: seq_id not found: " + row[kSeqId]);
    // Gemmi❗✔️:         }
    // Gemmi❗✔️:       }
    // Gemmi❗✔️:       unp.res = cif::as_char(row[kUnpRes], '\0');
    // Gemmi❗✔️:       int num = cif::as_int(row[kUnpNum]);
    // Gemmi❗✔️:       unp.num = (std::uint16_t) num;
    // Gemmi❗✔️:       if (num != (int)unp.num)
    // Gemmi❗✔️:         fail("_pdbx_sifts_xref_db.unp_num: " + row[kUnpNum]);
    // Gemmi❗✔️:       while (res->label_seq == label_seq && res != polymer.end()) {
    // Gemmi❗✔️:         res->sifts_unp = unp;
    // Gemmi❗✔️:         ++res;
    // Gemmi❗✔️:       }
    // Gemmi❗✔️:     }
    // Gemmi❗✔️:   }
    // Gemmi❗✔️: }
    // Gemmi✔️✔️: bool ok() const { return !positions.empty(); }
    // Gemmi❗✔️: inline std::string& Table::Row::operator[](size_t n) {
    // Gemmi❗✔️:   int pos = tab.positions[n];
    // Gemmi❗✔️:   if (Loop* loop = tab.get_loop()) {
    // Gemmi❗✔️:     if (row_index == -1) // tags
    // Gemmi❗✔️:       return loop->tags[pos];
    // Gemmi❗✔️:     return loop->values[loop->width() * row_index + pos];
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   return tab.bloc.items[pos].pair[row_index == -1 ? 0 : 1];
    // Gemmi❗✔️: }
    // Gemmi✔️✔️: inline bool is_null(const std::string& value) {
    // Gemmi✔️✔️:   return value.size() == 1 && (value[0] == '?' || value[0] == '.');
    // Gemmi✔️✔️: }
    // Gemmi✔️✔️: inline std::string as_string(const std::string& value) {
    // Gemmi✔️✔️:   if (value.empty() || is_null(value))
    // Gemmi✔️✔️:     return "";
    // Gemmi✔️✔️:   if (value[0] == '"' || value[0] == '\'')
    // Gemmi✔️✔️:     return std::string(value.begin() + 1, value.end() - 1);
    // Gemmi✔️✔️:   if (value[0] == ';' && value.size() > 2 && *(value.end() - 2) == '\n') {
    // Gemmi✔️✔️:     bool crlf = *(value.end() - 3) == '\r';
    // Gemmi✔️✔️:     return std::string(value.begin() + 1, value.end() - (crlf ? 3 : 2));
    // Gemmi✔️✔️:   }
    // Gemmi✔️✔️:   return value;
    // Gemmi✔️✔️: }
    // Gemmi❗✔️: inline char as_char(const std::string& value, char null) {
    // Gemmi❗✔️:   if (is_null(value))
    // Gemmi❗✔️:     return null;
    // Gemmi❗✔️:   if (value.size() < 2)
    // Gemmi❗✔️:     return value[0];
    // Gemmi❗✔️:   const std::string s = as_string(value);
    // Gemmi❗✔️:   if (s.size() < 2)
    // Gemmi❗✔️:     return s[0];
    // Gemmi❗✔️:   fail("Not a single character: " + value);
    // Gemmi❗✔️: }
    // Gemmi❗✔️: inline int as_int(const std::string& str) {
    // Gemmi❗✔️:   return string_to_int(str, true);
    // Gemmi❗✔️: }
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
    // Gemmi✔️✔️: inline bool is_digit(char c) {
    // Gemmi✔️✔️:   return c >= '0' && c <= '9';
    // Gemmi✔️✔️: }
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
    // Gemmi❗✔️:   for (; (length == 0 || i < length) && is_digit(p[i]); ++i) {
    // Gemmi❗✔️:     n = n * 10 - (p[i] - '0');
    // Gemmi❗✔️:     has_digits = true;
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   if (checked) {
    // Gemmi❗✔️:     while ((length == 0 || i < length) && is_space(p[i]))
    // Gemmi❗✔️:       ++i;
    // Gemmi❗✔️:     if (!has_digits || p[i] != '\0')
    // Gemmi❗✔️:       throw std::invalid_argument("not an integer: " +
    // Gemmi❗✔️:                                   std::string(p, length ? length : i+1));
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   return mult * n;
    // Gemmi❗✔️: }
    // Gemmi❗✔️: inline int string_to_int(const std::string& str, bool checked) {
    // Gemmi❗✔️:   return string_to_int(str.c_str(), checked);
    // Gemmi❗✔️: }
    // Gemmi✔️✔️: inline int as_int(const std::string& str, int null) {
    // Gemmi✔️✔️:   return is_null(str) ? null : as_int(str);
    // Gemmi✔️✔️: }
    // Gemmi✔️✔️: template<typename T>
    // Gemmi✔️✔️: auto get_id(const T& m) -> decltype(m.name) {
    // Gemmi✔️✔️:   return m.name;
    // Gemmi✔️✔️: }
    // Gemmi❗✔️: template<typename Vec, typename S>
    // Gemmi❗✔️: auto find_iter_(Vec& vec, const S& name) {
    // Gemmi❗✔️:   return std::find_if(vec.begin(), vec.end(), [&name](const auto& m) { return get_id(m) == name; });
    // Gemmi❗✔️: }
    // Gemmi❗✔️: template<typename T, typename S>
    // Gemmi❗✔️: T* find_or_null(std::vector<T>& vec, const S& name) {
    // Gemmi❗✔️:   auto it = find_iter_(vec, name);
    // Gemmi❗✔️:   return it != vec.end() ? &*it : nullptr;
    // Gemmi❗✔️: }
    // Gemmi❗✔️: Entity* Structure::get_entity(const std::string& ent_id) {
    // Gemmi❗✔️:   return impl::find_or_null(entities, ent_id);
    // Gemmi❗✔️: }
    // Gemmi❗✔️: ResidueSpan Model::get_subchain(const std::string& sub_name) {
    // Gemmi❗✔️:   for (Chain& chain : chains)
    // Gemmi❗✔️:     if (ResidueSpan sub = chain.get_subchain(sub_name))
    // Gemmi❗✔️:       return sub;
    // Gemmi❗✔️:   return ResidueSpan();
    // Gemmi❗✔️: }
    // Gemmi❗✔️: ResidueSpan Chain::get_subchain(const std::string& s) {
    // Gemmi❗✔️:   return get_residue_span([&](const Residue& r) { return r.subchain == s; });
    // Gemmi❗✔️: }
    // Gemmi❗✔️: template<typename F> ResidueSpan Chain::get_residue_span(F&& func) {
    // Gemmi❗✔️:   return whole().subspan(func);
    // Gemmi❗✔️: }
    // Gemmi❗✔️: template<typename F, typename V=Item> Span<V> Span::subspan(F&& func) {
    // Gemmi❗✔️:   iterator group_begin = std::find_if(this->begin(), this->end(), func);
    // Gemmi❗✔️:   iterator group_end = std::find_if_not(group_begin, this->end(), func);
    // Gemmi❗✔️:   return Span<V>(&*group_begin, group_end - group_begin);
    // Gemmi❗✔️: }
    // Behavior review: table and row cells remain in Gemmi's raw form for
    // ordinal/observed tests, null screening, entity/asym comparisons, and the
    // raw-accession cache comparison; only `row.str(unp_acc)` decodes the
    // accession. First-match entity/subchain lookup, source-order overwrites,
    // per-model caches, contiguous label-sequence assignment, source `u8`
    // accession narrowing, and the checked `u16` roundtrip are retained.
    // Rust checks span/cursor bounds before indexing. Gemmi's source expression
    // dereferences an end iterator before its end test, `Span::subspan` forms
    // `&*end` for no match, and `as_char`/checked `string_to_int` include
    // undefined empty/overflow cases; those are not copied as unsafe behavior.
    // Complexity review: table rows are rescanned once per model; entity,
    // subchain, and accession lookups are the same source-order linear scans;
    // each sequence seek restarts at the selected span's beginning. This uses
    // no span-vector allocation or hash index and stores assignments directly
    // on the existing residue staging rows.
    let table = block.find(
        "_pdbx_sifts_xref_db.",
        &[
            "entity_id",
            "asym_id",
            "seq_id_ordinal",
            "seq_id",
            "observed",
            "unp_res",
            "unp_num",
            "unp_acc",
        ],
    )?;
    if !table.is_present() {
        return Ok(());
    }

    for model in models {
        let mut entity_index: Option<usize> = None;
        let mut polymer_span: Option<(usize, usize, usize)> = None;
        let mut residue_cursor = None;
        let mut unp_accession = String::new();
        let mut accession_index = 0_u8;

        for row in table.iter() {
            let raw = |column| row.get(column).map(CifValue::raw).unwrap_or_default();
            let raw_ordinal = raw(2);
            let raw_observed = raw(4);
            if raw_ordinal != "1" || raw_observed.as_bytes().first() != Some(&b'y') {
                continue;
            }

            let raw_accession = raw(7);
            let raw_number = raw(6);
            if cif_is_null(raw_accession) || cif_is_null(raw_number) {
                continue;
            }

            let raw_entity_id = raw(0);
            let previous_entity_index = entity_index;
            let current_entity_index = match entity_index {
                Some(index)
                    if entities
                        .get(index)
                        .is_some_and(|entity| entity.source_entity_id == raw_entity_id) =>
                {
                    index
                }
                _ => entities
                    .iter()
                    .position(|entity| entity.source_entity_id == raw_entity_id)
                    .ok_or_else(|| MmcifSiftsUnpError::EntityNotFound {
                        entity_id: raw_entity_id.to_owned(),
                    })?,
            };
            entity_index = Some(current_entity_index);
            let mut update_accession_index = previous_entity_index != Some(current_entity_index);

            if raw_accession != unp_accession {
                unp_accession = cif_as_string(raw_accession);
                update_accession_index = true;
            }

            if update_accession_index {
                let entity = entities.get_mut(current_entity_index).ok_or(
                    MmcifSiftsUnpError::StagingInvariant {
                        stage: "resolved entity index is outside the entity rows",
                    },
                )?;
                let found = entity
                    .sifts_unp_accessions
                    .iter()
                    .position(|accession| accession == &unp_accession);
                accession_index = if let Some(index) = found {
                    index as u8
                } else {
                    let index = entity.sifts_unp_accessions.len() as u8;
                    entity.sifts_unp_accessions.push(unp_accession.clone());
                    index
                };
            }

            let raw_asym_id = raw(1);
            let span_has_same_asym = polymer_span.is_some_and(|(chain_index, start, _)| {
                model
                    .chains
                    .get(chain_index)
                    .and_then(|chain| chain.residues.get(start))
                    .is_some_and(|residue| {
                        residue.source.subchain_id().unwrap_or_default() == raw_asym_id
                    })
            });
            let (chain_index, span_start, span_end) = if span_has_same_asym {
                polymer_span.ok_or(MmcifSiftsUnpError::StagingInvariant {
                    stage: "matching polymer span was not retained",
                })?
            } else {
                let mut found_span = None;
                for (chain_index, chain) in model.chains.iter().enumerate() {
                    let Some(start) = chain.residues.iter().position(|residue| {
                        residue.source.subchain_id().unwrap_or_default() == raw_asym_id
                    }) else {
                        continue;
                    };
                    let mut end = start + 1;
                    while end < chain.residues.len()
                        && chain.residues[end].source.subchain_id().unwrap_or_default()
                            == raw_asym_id
                    {
                        end += 1;
                    }
                    found_span = Some((chain_index, start, end));
                    break;
                }
                let span = found_span.ok_or_else(|| MmcifSiftsUnpError::SubchainNotFound {
                    asym_id: raw_asym_id.to_owned(),
                })?;
                polymer_span = Some(span);
                residue_cursor = Some(span.1);
                span
            };
            let label_sequence_raw = raw(3);
            let label_sequence = cif_as_i32(label_sequence_raw)?;
            let mut cursor = residue_cursor.unwrap_or(span_start);
            if cursor == span_end {
                cursor = span_start;
            }
            let current_residue = model
                .chains
                .get(chain_index)
                .and_then(|chain| chain.residues.get(cursor))
                .ok_or(MmcifSiftsUnpError::StagingInvariant {
                    stage: "residue cursor is outside the selected polymer span",
                })?;
            if current_residue.source.label_seq_id().unwrap_or(i32::MIN) != label_sequence {
                cursor = span_start;
                while cursor < span_end
                    && model
                        .chains
                        .get(chain_index)
                        .and_then(|chain| chain.residues.get(cursor))
                        .is_some_and(|residue| {
                            residue.source.label_seq_id().unwrap_or(i32::MIN) != label_sequence
                        })
                {
                    cursor += 1;
                }
                if cursor == span_end {
                    return Err(MmcifSiftsUnpError::SequenceNotFound {
                        seq_id: label_sequence_raw.to_owned(),
                    });
                }
            }

            let residue_code = cif_as_char(raw(5), '\0')?;
            let number = cif_as_i32(raw_number)?;
            let narrowed_number = number as u16;
            if i32::from(narrowed_number) != number {
                return Err(MmcifSiftsUnpError::NumberOutsideU16 {
                    value: raw_number.to_owned(),
                });
            }
            let unp = BioSiftsUnpResidue::new(
                (residue_code != '\0').then_some(residue_code as u8),
                accession_index,
                narrowed_number,
            );

            while cursor < span_end {
                let residue = model
                    .chains
                    .get_mut(chain_index)
                    .and_then(|chain| chain.residues.get_mut(cursor))
                    .ok_or(MmcifSiftsUnpError::StagingInvariant {
                        stage: "assignment cursor is outside the selected polymer span",
                    })?;
                if residue.source.label_seq_id().unwrap_or(i32::MIN) != label_sequence {
                    break;
                }
                residue.sifts_unp = unp;
                cursor += 1;
            }
            residue_cursor = Some(cursor);
        }
    }
    Ok(())
}

fn atom_site_raw_value(row: CifRow<'_>, column: usize) -> &str {
    // Gemmi❗✔️: inline std::string& Table::Row::operator[](size_t n) {
    // Gemmi❗✔️:   int pos = tab.positions[n];
    // Gemmi❗✔️:   if (Loop* loop = tab.get_loop()) {
    // Gemmi❗✔️:     if (row_index == -1) // tags
    // Gemmi❗✔️:       return loop->tags[pos];
    // Gemmi❗✔️:     return loop->values[loop->width() * row_index + pos];
    // Gemmi❗✔️:   }
    // Gemmi❗✔️:   return tab.bloc.items[pos].pair[row_index == -1 ? 0 : 1];
    // Gemmi❗✔️: }
    // Behavior review: an absent Rust pair value is exposed as the source's
    // empty second string; a missing column is still distinguished by row.has.
    // Complexity review: one fixed-index lookup with no allocation.
    row.get(column).map_or("", CifValue::raw)
}

fn atom_site_has_source_value(row: CifRow<'_>, column: usize) -> bool {
    // Gemmi❗✔️: bool has(size_t n) const { return tab.positions.at(n) >= 0; }
    // Gemmi❗✔️: bool has2(size_t n) const { return has(n) && !cif::is_null(operator[](n)); }
    // Behavior review: `None` from the canonical row denotes the source's
    // empty string for a present pair, which is non-null under Gemmi's test.
    // Complexity review: one column-presence check and one fixed value read.
    row.has(column) && !cif_is_null(atom_site_raw_value(row, column))
}

fn has_mmcif_deuterium_fraction_column(atom_table: &CifTable<'_>) -> bool {
    // Gemmi❗✔️:         st.has_d_fraction = atom_table.has_column(kDeuterium);
    // Behavior review: `read_atom_sites` assigns this structure flag only
    // inside its nonempty atom-table branch, and bases it on column presence,
    // not any row's fraction value. The explicit empty check retains the
    // default false state when no atom-site table is selected.
    // Complexity review: both the table-length and fixed-column checks are
    // constant-time; this adds no row scan or value conversion.
    !atom_table.is_empty() && atom_table.has_column(AtomSiteColumn::DeuteriumFraction.index())
}

#[derive(Debug, Clone, PartialEq)]
struct MmcifAtomSiteScalars {
    serial: PdbAtomSerial,
    element: Element,
    isotope_mass_number: Option<u16>,
    altloc: Option<AltLocLabel>,
    formal_charge: i8,
    position: [f64; 3],
    occupancy: f32,
    b_iso: f32,
    calc_flag: BioCalcFlag,
    tls_group_id: i16,
    fraction: f32,
}

#[derive(Debug, Clone, PartialEq)]
struct MmcifParsedAtomSite {
    name: AtomName,
    raw_id: String,
    scalars: MmcifAtomSiteScalars,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum MmcifAtomSiteScalarError {
    SerialOutsideSourceDefinedRange { value: String },
    TlsOutsideSourceDefinedRange { value: String },
    CifCharacter(CifReadError),
    CifInteger(CifReadError),
    CifNumber(CifReadError),
}

impl fmt::Display for MmcifAtomSiteScalarError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::SerialOutsideSourceDefinedRange { value } => write!(
                f,
                "mmCIF atom-site id is outside Gemmi's source-defined integer range: {value:?}"
            ),
            Self::TlsOutsideSourceDefinedRange { value } => write!(
                f,
                "mmCIF atom-site TLS id is outside Gemmi's source-defined integer range: {value:?}"
            ),
            Self::CifCharacter(error) => fmt::Display::fmt(error, f),
            Self::CifInteger(error) => fmt::Display::fmt(error, f),
            Self::CifNumber(error) => fmt::Display::fmt(error, f),
        }
    }
}

impl std::error::Error for MmcifAtomSiteScalarError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::CifCharacter(error) | Self::CifInteger(error) | Self::CifNumber(error) => {
                Some(error)
            }
            Self::SerialOutsideSourceDefinedRange { .. }
            | Self::TlsOutsideSourceDefinedRange { .. } => None,
        }
    }
}

fn project_mmcif_atom_site_scalars(
    row: CifRow<'_>,
    has_d_fraction: bool,
) -> Result<MmcifAtomSiteScalars, MmcifAtomSiteScalarError> {
    // Gemmi❗❌:             atom.altloc = cif::as_char(row[kAltId], '\0');
    // Gemmi❗❌:             atom.charge = row.has2(kCharge) ? cif::as_int(row[kCharge]) : 0;
    // Gemmi❗❌:             atom.element = gemmi::Element(cif::as_string(row[kSymbol]));
    // Gemmi❗❌:             atom.serial = string_to_int(row[kId], false);
    // Gemmi❗❌:             if (st.has_d_fraction)
    // Gemmi❗❌:                 atom.fraction = (float) cif::as_number(row[kDeuterium], 0.);
    // Gemmi❗❌:             if (row.has2(kCalcFlag)) {
    // Gemmi❗❌:                 const std::string& cf = row[kCalcFlag];
    // Gemmi❗❌:                 if (cf[0] == 'c')
    // Gemmi❗❌:                     atom.calc_flag = CalcFlag::Calculated;
    // Gemmi❗❌:                 if (cf[0] == 'd')
    // Gemmi❗❌:                     atom.calc_flag = cf[1] == 'u' ? CalcFlag::Dummy
    // Gemmi❗❌:                                                   : CalcFlag::Determined;
    // Gemmi❗❌:             }
    // Gemmi❗❌:             if (row.has2(kTlsGroupId)) {
    // Gemmi❗❌:                 const char* str = row[kTlsGroupId].c_str();
    // Gemmi❗❌:                 const char* endptr;
    // Gemmi❗❌:                 int tls_id = no_sign_atoi(str, &endptr);
    // Gemmi❗❌:                 if (endptr != str)
    // Gemmi❗❌:                     atom.tls_group_id = (short) tls_id;
    // Gemmi❗❌:             }
    // Gemmi❗❌:             atom.pos.x = cif::as_number(row[kX]);
    // Gemmi❗❌:             atom.pos.y = cif::as_number(row[kY]);
    // Gemmi❗❌:             atom.pos.z = cif::as_number(row[kZ]);
    // Gemmi❗❌:             if (row.has2(kOcc))
    // Gemmi❗❌:                 atom.occ = (float) cif::as_number(row[kOcc]);
    // Gemmi❗❌:             if (row.has2(kBiso))
    // Gemmi❗❌:                 atom.b_iso = (float) cif::as_number(row[kBiso]);
    // Behavior review: CIF null handling is intentionally field-specific:
    // serial is raw prefix text with no-conversion zero, charge uses checked
    // complete integer parsing under source `has2`, coordinates always call
    // `as_number` with NaN fallback, occupancy/B use source defaults only
    // when absent/null, and deuterium fraction uses zero fallback only when
    // the column exists. Source D is represented canonically as H+Some(2).
    // TLS preserves the source end-pointer condition; source int overflow is
    // returned as a typed boundary rather than inventing a value. The source
    // signed-char/short assignments narrow to the corresponding Rust fields.
    let altloc_character = cif_as_char(
        atom_site_raw_value(row, AtomSiteColumn::AltId.index()),
        '\0',
    )
    .map_err(MmcifAtomSiteScalarError::CifCharacter)?;
    let altloc = (altloc_character != '\0').then(|| AltLocLabel::new(altloc_character as u8));

    let charge_column = AtomSiteColumn::FormalCharge.index();
    let formal_charge = if atom_site_has_source_value(row, charge_column) {
        cif_as_i32(atom_site_raw_value(row, charge_column))
            .map_err(MmcifAtomSiteScalarError::CifInteger)? as i8
    } else {
        0
    };

    let symbol = cif_as_string(atom_site_raw_value(row, AtomSiteColumn::TypeSymbol.index()));
    let symbol_bytes = symbol.as_bytes();
    let (element, isotope_mass_number) = super::bio_pdb::gemmi_find_element([
        symbol_bytes.first().copied().unwrap_or_default(),
        symbol_bytes.get(1).copied().unwrap_or_default(),
    ]);

    let id_column = AtomSiteColumn::Id.index();
    let id = atom_site_raw_value(row, id_column);
    let serial = super::bio_pdb::read_int(id.as_bytes()).ok_or_else(|| {
        MmcifAtomSiteScalarError::SerialOutsideSourceDefinedRange {
            value: id.to_owned(),
        }
    })?;

    let fraction = if has_d_fraction {
        number_with_null(row, AtomSiteColumn::DeuteriumFraction.index(), 0.0)? as f32
    } else {
        0.0
    };

    let calc_flag_column = AtomSiteColumn::CalcFlag.index();
    let mut calc_flag = BioCalcFlag::NotSet;
    if atom_site_has_source_value(row, calc_flag_column) {
        let raw = atom_site_raw_value(row, calc_flag_column).as_bytes();
        let first = raw.first().copied().unwrap_or_default();
        if first == b'c' {
            calc_flag = BioCalcFlag::Calculated;
        }
        if first == b'd' {
            calc_flag = if raw.get(1).copied().unwrap_or_default() == b'u' {
                BioCalcFlag::Dummy
            } else {
                BioCalcFlag::Determined
            };
        }
    }

    let tls_column = AtomSiteColumn::TlsGroupId.index();
    let mut tls_group_id = -1_i16;
    if atom_site_has_source_value(row, tls_column) {
        let raw = atom_site_raw_value(row, tls_column);
        let (tls_id, end_offset) =
            super::bio_pdb::gemmi_no_sign_atoi(raw.as_bytes()).ok_or_else(|| {
                MmcifAtomSiteScalarError::TlsOutsideSourceDefinedRange {
                    value: raw.to_owned(),
                }
            })?;
        if end_offset != 0 {
            tls_group_id = tls_id as i16;
        }
    }

    let number = |column| {
        cif_as_f64(atom_site_raw_value(row, column), f64::NAN)
            .map_err(MmcifAtomSiteScalarError::CifNumber)
    };
    let position = [
        number(AtomSiteColumn::X.index())?,
        number(AtomSiteColumn::Y.index())?,
        number(AtomSiteColumn::Z.index())?,
    ];

    let occupancy_column = AtomSiteColumn::Occupancy.index();
    let occupancy = if atom_site_has_source_value(row, occupancy_column) {
        number(occupancy_column)? as f32
    } else {
        1.0
    };
    let b_iso_column = AtomSiteColumn::BIso.index();
    let b_iso = if atom_site_has_source_value(row, b_iso_column) {
        number(b_iso_column)? as f32
    } else {
        20.0
    };

    Ok(MmcifAtomSiteScalars {
        serial: PdbAtomSerial::new(serial),
        element,
        isotope_mass_number,
        altloc,
        formal_charge,
        position,
        occupancy,
        b_iso,
        calc_flag,
        tls_group_id,
        fraction,
    })
}

fn number_with_null(
    row: CifRow<'_>,
    column: usize,
    null: f64,
) -> Result<f64, MmcifAtomSiteScalarError> {
    cif_as_f64(atom_site_raw_value(row, column), null).map_err(MmcifAtomSiteScalarError::CifNumber)
}

fn extract_mmcif_anisotropic_u(
    block: &CifBlock,
) -> Result<HashMap<String, [f32; 6]>, CifReadError> {
    // Gemmi❗❌: template<typename T>
    // Gemmi❗❌: SMat33<T> get_smat33(cif::Table::Row& row, int n) {
    // Gemmi❗❌:   return SMat33<T>{(T) cif::as_number(row[n+0]),
    // Gemmi❗❌:                    (T) cif::as_number(row[n+1]),
    // Gemmi❗❌:                    (T) cif::as_number(row[n+2]),
    // Gemmi❗❌:                    (T) cif::as_number(row[n+3]),
    // Gemmi❗❌:                    (T) cif::as_number(row[n+4]),
    // Gemmi❗❌:                    (T) cif::as_number(row[n+5])};
    // Gemmi❗❌: }
    // Gemmi❗❌: std::unordered_map<std::string, SMat33<float>> get_anisotropic_u(cif::Block& block) {
    // Gemmi❗❌:   cif::Table aniso_tab = block.find("_atom_site_anisotrop.",
    // Gemmi❗❌:                                     {"id", "U[1][1]", "U[2][2]", "U[3][3]",
    // Gemmi❗❌:                                      "U[1][2]", "U[1][3]", "U[2][3]"});
    // Gemmi❗❌:   std::unordered_map<std::string, SMat33<float>> aniso_map;
    // Gemmi❗❌:   for (auto ani : aniso_tab)
    // Gemmi❗❌:     aniso_map.emplace(ani[0], get_smat33<float>(ani, 1));
    // Gemmi❗❌:   return aniso_map;
    // Gemmi❗❌: }
    // Behavior review: keep the raw source ID as the key; `emplace` means the
    // first duplicate wins, and each row's six `as_number` conversions occur
    // before insertion. Components retain SMat33 PDB order and narrow to f32;
    // invalid/null values use the canonical CIF NaN fallback. The later atom
    // consumer must look up the unchanged raw ID, not the parsed serial.
    // Complexity review: one table discovery plus one ordered pass over rows,
    // six O(token-bytes) numeric conversions per row, and expected O(1) map
    // insertion/lookup; keys allocate once and rows are borrowed, without
    // cloning the table. The canonical numeric helper has extra linear scans
    // relative to Gemmi's direct fast_float parse, so the cost axis is marked
    // worse despite matching expected hash-map complexity.
    const ANISOTROPIC_TAGS: [&str; 7] = [
        "id", "U[1][1]", "U[2][2]", "U[3][3]", "U[1][2]", "U[1][3]", "U[2][3]",
    ];
    let table = block.find("_atom_site_anisotrop.", &ANISOTROPIC_TAGS)?;
    let mut anisotropic = HashMap::new();
    for row in table.iter() {
        let raw_id = atom_site_raw_value(row, 0).to_owned();
        let component = |column| {
            cif_as_f64(atom_site_raw_value(row, column), f64::NAN).map(|value| value as f32)
        };
        let values = [
            component(1)?,
            component(2)?,
            component(3)?,
            component(4)?,
            component(5)?,
            component(6)?,
        ];
        anisotropic.entry(raw_id).or_insert(values);
    }
    Ok(anisotropic)
}

fn group_mmcif_atom_sites(
    atom_table: &CifTable<'_>,
) -> Result<MmcifAtomSiteGrouping, MmcifAtomSiteGroupingError> {
    group_mmcif_atom_sites_with(atom_table, |_, _| Ok(())).map(|(grouping, _)| grouping)
}

fn group_mmcif_atom_sites_with<T, F>(
    atom_table: &CifTable<'_>,
    mut process_atom: F,
) -> Result<(MmcifAtomSiteGrouping, Vec<T>), MmcifAtomSiteGroupingError>
where
    F: for<'row> FnMut(
        CifRow<'row>,
        AtomSiteIdentityColumns,
    ) -> Result<T, MmcifAtomSiteGroupingError>,
{
    // Gemmi❗❌:         Model *model = nullptr;
    // Gemmi❗❌:         Chain *chain = nullptr;
    // Gemmi❗❌:         Residue *resi = nullptr;
    // Gemmi❗❌:         std::string model_num;
    // Gemmi❗❌:         if (!atom_table.has_column(kModelNum)) {
    // Gemmi❗❌:             st.models.emplace_back(1);
    // Gemmi❗❌:             model = &st.models[0];
    // Gemmi❗❌:         }
    // Gemmi❗❌:         for (auto row : atom_table) {
    // Gemmi❗❌:             if (row.has(kModelNum) && row[kModelNum] != model_num) {
    // Gemmi❗❌:                 model_num = row[kModelNum];
    // Gemmi❗❌:                 model = &st.find_or_add_model(cif::as_int(model_num, 0));
    // Gemmi❗❌:                 chain = nullptr;
    // Gemmi❗❌:             }
    // Gemmi❗❌:             if (!chain || cif::as_string(asym_id.get(gap)) != chain->name) {
    // Gemmi❗❌:                 model->chains.emplace_back(cif::as_string(asym_id.get(gap)));
    // Gemmi❗❌:                 chain = &model->chains.back();
    // Gemmi❗❌:                 resi = nullptr;
    // Gemmi❗❌:             }
    // Gemmi❗❌:             ResidueId rid = make_resid(cif::as_string(comp_id.get(gap)),
    // Gemmi❗❌:                                        cif::as_string(seq_id.get(gap)),
    // Gemmi❗❌:                                        row.has(kInsCode) ? &row[kInsCode] : nullptr);
    // Gemmi❗❌:             if (!resi || !resi->matches(rid)) {
    // Gemmi❗❌:                 resi = chain->find_or_add_residue(rid);
    // Gemmi❗❌:                 if (resi->atoms.empty()) {
    // Gemmi❗❌:                     if (row.has2(kLabelSeqId))
    // Gemmi❗❌:                         resi->label_seq = cif::as_int(row[kLabelSeqId]);
    // Gemmi❗❌:                     resi->subchain = row.str(kLabelAsymId);
    // Gemmi❗❌:                     if (row.has2(kLabelEntityId))
    // Gemmi❗❌:                         resi->entity_id = row.str(kLabelEntityId);
    // Gemmi❗❌:                     // don't check if group_PDB is consistent, it's not that important
    // Gemmi❗❌:                     if (row.has2(kGroupPdb))
    // Gemmi❗❌:                         for (int i = 0; i < 2; ++i) { // first character could be " or '
    // Gemmi❗❌:                             const char c = alpha_up(row[kGroupPdb][i]);
    // Gemmi❗❌:                             if (c == 'A' || c == 'H' || c == '\0')
    // Gemmi❗❌:                                 resi->het_flag = c;
    // Gemmi❗❌:                         }
    // Gemmi❗❌:                 }
    // Gemmi❗❌:             } else if (resi->seqid != rid.seqid) {
    // Gemmi❗❌:                 fail("Inconsistent sequence ID: " + resi->str() + " / " + rid.str());
    // Gemmi❗❌:             }
    // Gemmi❗✔️: Model& find_or_add_model(int model_num) {
    // Gemmi❗✔️:   return impl::find_or_add(models, model_num);
    // Gemmi❗✔️: }
    // Gemmi❗✔️: T* find_or_null(std::vector<T>& vec, const S& name) {
    // Gemmi❗✔️:   auto it = find_iter_(vec, name);
    // Gemmi❗✔️:   return it != vec.end() ? &*it : nullptr;
    // Gemmi❗✔️: }
    // Gemmi❗✔️: T& find_or_add(std::vector<T>& vec, const S& name) {
    // Gemmi❗✔️:   if (T* ret = find_or_null(vec, name))
    // Gemmi❗✔️:     return *ret;
    // Gemmi❗✔️:   vec.emplace_back(name);
    // Gemmi❗✔️:   return vec.back();
    // Gemmi❗✔️: }
    // Gemmi❗✔️: inline Residue* Chain::find_residue(const ResidueId& rid) {
    // Gemmi❗✔️:   auto it = std::find_if(residues.begin(), residues.end(),
    // Gemmi❗✔️:                          [&](const Residue& r) { return r.matches(rid); });
    // Gemmi❗✔️:   return it != residues.end() ? &*it : nullptr;
    // Gemmi❗✔️: }
    // Gemmi❗✔️: inline Residue* Chain::find_or_add_residue(const ResidueId& rid) {
    // Gemmi❗✔️:   Residue* r = find_residue(rid);
    // Gemmi❗✔️:   if (r)
    // Gemmi❗✔️:     return r;
    // Gemmi❗✔️:   residues.emplace_back(rid);
    // Gemmi❗✔️:   return &residues.back();
    // Gemmi❗✔️: }
    // Behavior review: implementation is not yet accepted until the named
    // C04 regressions run. Preserve Gemmi's exact raw model-token transition test,
    // numeric model reuse, active-chain reset, append-on-chain-transition,
    // and source-order residue reuse. The grouping tree retains auth and label
    // chain identities plus first-created residue label/entity metadata; row
    // indices remain in their matched residue. Canonical BIO width limits
    // produce typed errors and never truncate source identifiers.
    // Complexity review: model and residue identity use the same linear
    // first-match scans as Gemmi; chain transitions compare only the active
    // chain. Staging also retains an intermediate hierarchy and O(A) row
    // indices before later scalar projection, so it uses extra per-row memory
    // compared with Gemmi's direct atom insertion.
    if atom_table.is_empty() {
        return Ok((MmcifAtomSiteGrouping::default(), Vec::new()));
    }

    let identities = bind_atom_site_identity_columns(atom_table)
        .map_err(MmcifAtomSiteGroupingError::IdentityField)?;
    let model_column = AtomSiteColumn::ModelNumber.index();
    let group_pdb_column = AtomSiteColumn::GroupPdb.index();
    let label_entity_column = AtomSiteColumn::LabelEntityId.index();
    let label_sequence_column = AtomSiteColumn::LabelSeqId.index();
    let insertion_code_column = AtomSiteColumn::InsertionCode.index();
    let label_asym_column = AtomSiteColumn::LabelAsymId.index();
    let auth_asym_column = AtomSiteColumn::AuthAsymId.index();
    let has_model_column = atom_table.has_column(model_column);

    let mut grouping = MmcifAtomSiteGrouping::default();
    let mut raw_model_token = String::new();
    let mut current_model_index: Option<usize> = None;
    let mut current_chain_index: Option<usize> = None;
    let mut current_residue_index: Option<usize> = None;

    if !has_model_column {
        grouping.models.push(MmcifModelGrouping {
            source_model_number: 1,
            chains: Vec::new(),
        });
        current_model_index = Some(0);
    }

    let mut parsed_atoms = Vec::with_capacity(atom_table.len());
    for (source_row_index, row) in atom_table.iter().enumerate() {
        if row.has(model_column) {
            let next_raw_model_token = atom_site_raw_value(row, model_column);
            if next_raw_model_token != raw_model_token {
                raw_model_token.clear();
                raw_model_token.push_str(next_raw_model_token);
                let source_model_number = if cif_is_null(next_raw_model_token) {
                    0
                } else {
                    cif_as_i32(next_raw_model_token)
                        .map_err(MmcifAtomSiteGroupingError::ModelNumber)?
                };
                let model_index = grouping
                    .models
                    .iter()
                    .position(|model| model.source_model_number == source_model_number)
                    .unwrap_or_else(|| {
                        grouping.models.push(MmcifModelGrouping {
                            source_model_number,
                            chains: Vec::new(),
                        });
                        grouping.models.len() - 1
                    });
                current_model_index = Some(model_index);
                current_chain_index = None;
            }
        }

        let model_index =
            current_model_index.ok_or(MmcifAtomSiteGroupingError::MissingActiveModel {
                row_index: source_row_index,
            })?;
        let chain_name = cif_as_string(identities.asym_id.get(row));
        let chain_index = match current_chain_index {
            Some(index) if grouping.models[model_index].chains[index].source_name == chain_name => {
                index
            }
            _ => {
                let raw_auth_chain_id = atom_site_raw_value(row, auth_asym_column);
                let auth_chain_id = if atom_site_has_source_value(row, auth_asym_column) {
                    let decoded = cif_as_string(raw_auth_chain_id);
                    Some(PdbChainId::from_ascii(decoded.as_bytes()).ok_or_else(|| {
                        MmcifAtomSiteGroupingError::AuthorChainIdNotRepresentable {
                            value: decoded.clone(),
                        }
                    })?)
                } else {
                    None
                };
                let label_chain_name = cif_as_string(atom_site_raw_value(row, label_asym_column));
                let source = ChainSourceIds::new(auth_chain_id, Some(label_chain_name));
                let chains = &mut grouping.models[model_index].chains;
                chains.push(MmcifChainGrouping {
                    source_name: chain_name,
                    source,
                    residues: Vec::new(),
                });
                current_residue_index = None;
                chains.len() - 1
            }
        };

        let residue_name_text = cif_as_string(identities.comp_id.get(row));
        let residue_name =
            ResidueName::from_ascii(residue_name_text.as_bytes()).ok_or_else(|| {
                MmcifAtomSiteGroupingError::ResidueNameNotRepresentable {
                    value: residue_name_text.clone(),
                }
            })?;
        let sequence_id = make_mmcif_seq_id(
            cif_as_string(identities.seq_id.get(row)),
            row.has(insertion_code_column)
                .then(|| atom_site_raw_value(row, insertion_code_column)),
        )
        .map_err(MmcifAtomSiteGroupingError::SequenceId)?;
        let address = make_mmcif_residue_address(residue_name, sequence_id).ok_or_else(|| {
            MmcifAtomSiteGroupingError::ResidueAddressNotRepresentable {
                value: residue_name_text.clone(),
            }
        })?;

        let residue_index = {
            let residues = &grouping.models[model_index].chains[chain_index].residues;
            current_residue_index
                .filter(|index| residues[*index].address.matches(&address))
                .or_else(|| {
                    residues
                        .iter()
                        .position(|residue| residue.address.matches(&address))
                })
        };
        let residue_index = match residue_index {
            Some(index) => index,
            None => {
                let label_sequence_number =
                    if atom_site_has_source_value(row, label_sequence_column) {
                        let value = cif_as_i32(atom_site_raw_value(row, label_sequence_column))
                            .map_err(MmcifAtomSiteGroupingError::LabelSequenceNumber)?;
                        (value != i32::MIN).then_some(value)
                    } else {
                        None
                    };
                let label_entity_id = atom_site_has_source_value(row, label_entity_column)
                    .then(|| cif_as_string(atom_site_raw_value(row, label_entity_column)));
                let subchain_id = Some(cif_as_string(atom_site_raw_value(row, label_asym_column)));
                let source = ResidueSourceIds::new(
                    Some(sequence_id),
                    label_sequence_number,
                    None,
                    subchain_id,
                    label_entity_id,
                )
                .ok_or(MmcifAtomSiteGroupingError::ResidueSourceIdsNotRepresentable)?;
                let het_flag = if atom_site_has_source_value(row, group_pdb_column) {
                    let raw = atom_site_raw_value(row, group_pdb_column).as_bytes();
                    let mut value = None;
                    for index in 0..2 {
                        let character = match raw.get(index).copied() {
                            Some(character) => character,
                            None if index == raw.len() => 0,
                            None => {
                                return Err(MmcifAtomSiteGroupingError::UndefinedGroupPdbIndex {
                                    value: atom_site_raw_value(row, group_pdb_column).to_owned(),
                                    index,
                                    length: raw.len(),
                                });
                            }
                        }
                        .to_ascii_uppercase();
                        if matches!(character, b'A' | b'H' | 0) {
                            value = Some(character);
                        }
                    }
                    value
                } else {
                    None
                };
                let residues = &mut grouping.models[model_index].chains[chain_index].residues;
                residues.push(MmcifResidueGrouping {
                    address,
                    source,
                    het_flag,
                    entity_kind: EntityKind::Unknown,
                    atom_site_rows: Vec::new(),
                    sifts_unp: BioSiftsUnpResidue::default(),
                });
                residues.len() - 1
            }
        };

        current_chain_index = Some(chain_index);
        current_residue_index = Some(residue_index);
        grouping.models[model_index].chains[chain_index].residues[residue_index]
            .atom_site_rows
            .push(source_row_index);
        parsed_atoms.push(process_atom(row, identities)?);
    }

    Ok((grouping, parsed_atoms))
}

#[derive(Debug)]
struct MmcifAtomSiteRowOutOfBounds {
    row_index: usize,
    parsed_count: usize,
}

impl fmt::Display for MmcifAtomSiteRowOutOfBounds {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "atom-site grouping references row {} but only {} rows were parsed",
            self.row_index, self.parsed_count
        )
    }
}

impl Error for MmcifAtomSiteRowOutOfBounds {}

#[derive(Debug)]
enum MmcifCcdRestoreError {
    Residue(PdbResidueKeyError),
    SequenceUtf8(std::string::FromUtf8Error),
}

impl fmt::Display for MmcifCcdRestoreError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Residue(error) => fmt::Display::fmt(error, f),
            Self::SequenceUtf8(error) => {
                write!(f, "restored entity sequence is not UTF-8: {error}")
            }
        }
    }
}

impl Error for MmcifCcdRestoreError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Residue(error) => Some(error),
            Self::SequenceUtf8(error) => Some(error),
        }
    }
}

fn parse_mmcif_shortened_ccd_codes(block: &CifBlock) -> Result<Vec<PdbCcdAlias>, CifReadError> {
    // Gemmi❗✔️:   cif::Table chem_comp_table = block.find("_chem_comp.", {"id", "three_letter_code"});
    // Gemmi❗✔️:   if (chem_comp_table.ok()) {
    // Gemmi❗✔️:     for (auto row : chem_comp_table) {
    // Gemmi❗✔️:       std::string alias = row.str(0);
    // Gemmi❗✔️:       std::string long_id = row.str(1);
    // Gemmi❗✔️:       if (alias[0] == '~' && long_id[0] != '~' && long_id[0] != '\0')
    // Gemmi❗✔️:         st.shortened_ccd_codes.emplace_back(long_id, alias);
    // Gemmi❗✔️:     }
    // Gemmi❗✔️:     restore_full_ccd_codes(st);
    // Gemmi❗✔️:   }
    // Behavior review: decoded CIF null/empty strings have the source NUL
    // first-byte value; rows are retained only for the exact tilde/nonempty
    // condition and aliases remain in source row order.
    // Complexity review: one table selection and one row pass; owned alias
    // bytes are copied once into the ordered restore list.
    let table = block.find("_chem_comp.", &["id", "three_letter_code"])?;
    let mut aliases = Vec::new();
    for row in table.iter() {
        let alias = row.decoded(0).unwrap_or_default();
        let full_code = row.decoded(1).unwrap_or_default();
        let alias_first = alias.as_bytes().first().copied().unwrap_or_default();
        let full_first = full_code.as_bytes().first().copied().unwrap_or_default();
        if alias_first == b'~' && full_first != b'~' && full_first != 0 {
            aliases.push(PdbCcdAlias {
                full_code: full_code.into_bytes(),
                short_code: alias.into_bytes(),
            });
        }
    }
    Ok(aliases)
}

fn restore_mmcif_ccd_codes(
    aliases: &[PdbCcdAlias],
    grouping: &mut MmcifAtomSiteGrouping,
    entities: &mut [MmcifEntityBuildState],
    connections: &mut [BioConnection],
    cispeps: &mut [BioCisPep],
    helices: &mut [BioHelix],
    sheets: &mut [BioSheet],
    mod_residues: &mut [BioModRes],
) -> Result<(), MmcifCcdRestoreError> {
    // Gemmi❗❌: void restore_full_ccd_codes(Structure& st) {
    // Gemmi❗❌:   for (const auto& item : st.shortened_ccd_codes)
    // Gemmi❗❌:     rename_residues(st, item.second, item.first);
    // Gemmi❗❌:   st.shortened_ccd_codes.clear();
    // Gemmi❗❌: }
    // Gemmi❗❌: inline void rename_residues(Structure& st, const std::string& old_name,
    // Gemmi❗❌:                                            const std::string& new_name) {
    // Gemmi❗❌:   auto update = [&](ResidueId& rid) {
    // Gemmi❗❌:     if (rid.name == old_name)
    // Gemmi❗❌:       rid.name = new_name;
    // Gemmi❗❌:   };
    // Gemmi❗❌:   process_addresses(st, [&](AtomAddress& aa) { update(aa.res_id); });
    // Gemmi❗❌:   for (ModRes& modres : st.mod_residues)
    // Gemmi❗❌:     update(modres.res_id);
    // Gemmi❗❌:   for (Entity& ent : st.entities)
    // Gemmi❗❌:     for (std::string& mon_ids : ent.full_sequence)
    // Gemmi❗❌:       for (size_t start = 0;;) {
    // Gemmi❗❌:         size_t end = mon_ids.find(',', start);
    // Gemmi❗❌:         if (mon_ids.compare(start, end-start, old_name) == 0) {
    // Gemmi❗❌:           mon_ids.replace(start, end-start, new_name);
    // Gemmi❗❌:           if (end != std::string::npos)
    // Gemmi❗❌:             end = start + new_name.size();
    // Gemmi❗❌:         }
    // Gemmi❗❌:         if (end == std::string::npos)
    // Gemmi❗❌:           break;
    // Gemmi❗❌:         start = end + 1;
    // Gemmi❗❌:       }
    // Gemmi❗❌:   for (Model& model : st.models)
    // Gemmi❗❌:     for (Chain& chain : model.chains)
    // Gemmi❗❌:       for (Residue& res : chain.residues)
    // Gemmi❗❌:         update(res);
    // Gemmi❗❌: }
    // Gemmi❗❌: template<typename Func>
    // Gemmi❗❌: void process_addresses(Structure& st, Func func) {
    // Gemmi❗❌:   for (Connection& con : st.connections) {
    // Gemmi❗❌:     func(con.partner1);
    // Gemmi❗❌:     func(con.partner2);
    // Gemmi❗❌:   }
    // Gemmi❗❌:   for (CisPep& cispep : st.cispeps) {
    // Gemmi❗❌:     func(cispep.partner_c);
    // Gemmi❗❌:     func(cispep.partner_n);
    // Gemmi❗❌:   }
    // Gemmi❗❌:   for (Helix& helix : st.helices) {
    // Gemmi❗❌:     func(helix.start);
    // Gemmi❗❌:     func(helix.end);
    // Gemmi❗❌:   }
    // Gemmi❗❌:   for (Sheet& sheet : st.sheets)
    // Gemmi❗❌:     for (Sheet::Strand& strand : sheet.strands) {
    // Gemmi❗❌:       func(strand.start);
    // Gemmi❗❌:       func(strand.end);
    // Gemmi❗❌:       func(strand.hbond_atom2);
    // Gemmi❗❌:       func(strand.hbond_atom1);
    // Gemmi❗❌:     }
    // Gemmi❗❌: }
    // Behavior review: aliases are applied in source order and update every
    // source-address family, MODRES, comma-separated entity sequence token,
    // and hierarchy residue. Existing canonical address-width failures remain
    // typed and the partially transformed detached staging value is discarded.
    // Complexity review: each alias traverses every represented address and
    // row, matching the source's O(alias_count * structure_size) passes; the
    // immutable AtomAddress model requires rebuilding its owned atom-name
    // strings, so this allocation path is slower than Gemmi's in-place writes.
    for alias in aliases {
        let current = std::slice::from_ref(alias);
        for connection in connections.iter_mut() {
            connection.partner1 = mapped_ccd_atom_address(&connection.partner1, current)
                .map_err(MmcifCcdRestoreError::Residue)?;
            connection.partner2 = mapped_ccd_atom_address(&connection.partner2, current)
                .map_err(MmcifCcdRestoreError::Residue)?;
        }
        for cispep in cispeps.iter_mut() {
            cispep.partner_c = mapped_ccd_atom_address(&cispep.partner_c, current)
                .map_err(MmcifCcdRestoreError::Residue)?;
            cispep.partner_n = mapped_ccd_atom_address(&cispep.partner_n, current)
                .map_err(MmcifCcdRestoreError::Residue)?;
        }
        for helix in helices.iter_mut() {
            helix.start = mapped_ccd_atom_address(&helix.start, current)
                .map_err(MmcifCcdRestoreError::Residue)?;
            helix.end = mapped_ccd_atom_address(&helix.end, current)
                .map_err(MmcifCcdRestoreError::Residue)?;
        }
        for sheet in sheets.iter_mut() {
            for strand in &mut sheet.strands {
                strand.start = mapped_ccd_atom_address(&strand.start, current)
                    .map_err(MmcifCcdRestoreError::Residue)?;
                strand.end = mapped_ccd_atom_address(&strand.end, current)
                    .map_err(MmcifCcdRestoreError::Residue)?;
                strand.hbond_atom2 = mapped_ccd_atom_address(&strand.hbond_atom2, current)
                    .map_err(MmcifCcdRestoreError::Residue)?;
                strand.hbond_atom1 = mapped_ccd_atom_address(&strand.hbond_atom1, current)
                    .map_err(MmcifCcdRestoreError::Residue)?;
            }
        }
        for mod_residue in mod_residues.iter_mut() {
            mod_residue.res_id = mapped_ccd_residue_address(mod_residue.res_id, current)
                .map_err(MmcifCcdRestoreError::Residue)?;
        }
        for entity in entities.iter_mut() {
            for sequence in &mut entity.full_sequence {
                let mut bytes = std::mem::take(sequence).into_bytes();
                rename_ccd_sequence_tokens(
                    &mut bytes,
                    alias.short_code.as_slice(),
                    alias.full_code.as_slice(),
                );
                *sequence = String::from_utf8(bytes).map_err(MmcifCcdRestoreError::SequenceUtf8)?;
            }
        }
        for model in &mut grouping.models {
            for chain in &mut model.chains {
                for residue in &mut chain.residues {
                    residue.address = mapped_ccd_residue_address(residue.address, current)
                        .map_err(MmcifCcdRestoreError::Residue)?;
                }
            }
        }
    }
    Ok(())
}

/// Read coordinate mmCIF text into one validated detached BIO structure.
pub fn read_mmcif_bio_structure(
    text: &str,
    source_name: &str,
) -> Result<BioStructureData, BioMmcifReadError> {
    let document = read_cif_document(text, source_name, CifCheckLevel::Default)
        .map_err(|error| BioMmcifReadError::new(BioMmcifReadStage::CifDocument, error))?;
    populate_mmcif_bio_structure_document(&document)
}

pub(crate) fn populate_mmcif_bio_structure_document(
    document: &CifDocument,
) -> Result<BioStructureData, BioMmcifReadError> {
    // Gemmi❗❌: void populate_structure_from_block(const cif::Block& block_, Structure& st) {
    // Gemmi❗❌:   // find() and Table don't have const variants, but we don't change anything.
    // Gemmi❗❌:   cif::Block& block = const_cast<cif::Block&>(block_);
    // Gemmi❗❌:   st.input_format = CoorFormat::Mmcif;
    // Gemmi❗❌:   st.name = block.name;
    // Gemmi❗❌:   impl::set_cell_from_mmcif(block, st.cell);
    // Gemmi❗❌:   st.spacegroup_hm = cif::as_string(impl::find_spacegroup_hm_value(block));
    // Gemmi❗❌:
    // Gemmi❗❌:   read_entry_info(block, st);
    // Gemmi❗❌:   read_audit_author(block, st);
    // Gemmi❗❌:   read_refinement_info(block, st);
    // Gemmi❗❌:   read_tls_info(block, st);
    // Gemmi❗❌:   read_experimental_info(block, st);
    // Gemmi❗❌:   read_reflns_info(block, st);
    // Gemmi❗❌:   read_software_info(block, st);
    // Gemmi❗❌:   read_ncs_info(block, st);
    // Gemmi❗❌:
    // Gemmi❗❌:   // PDBx/mmcif spec defines both _database_PDB_matrix.scale* and
    // Gemmi❗❌:   // _atom_sites.fract_transf_* as equivalent of pdb SCALE, but the former
    // Gemmi❗❌:   // is not used, so we ignore it.
    // Gemmi❗❌:   cif::Table fract_tv = block.find("_atom_sites.fract_transf_",
    // Gemmi❗❌:                                    transform_tags("matrix", "vector"));
    // Gemmi❗❌:   if (fract_tv.length() > 0) {
    // Gemmi❗❌:     Transform fract = get_transform_matrix(fract_tv[0]);
    // Gemmi❗❌:     st.cell.set_matrices_from_fract(fract);
    // Gemmi❗❌:   }
    // Gemmi❗❌:
    // Gemmi❗❌:   // We read/write origx just for completeness, it's not used anywhere.
    // Gemmi❗❌:   cif::Table origx_tv = block.find("_database_PDB_matrix.",
    // Gemmi❗❌:                                    transform_tags("origx", "origx_vector"));
    // Gemmi❗❌:   if (origx_tv.length() > 0) {
    // Gemmi❗❌:     st.has_origx = true;
    // Gemmi❗❌:     st.origx = get_transform_matrix(origx_tv[0]);
    // Gemmi❗❌:   }
    // Gemmi❗❌:
    // Gemmi❗❌:   read_atom_sites(block, st);
    // Gemmi❗❌:   read_entity_and_sequence_info(block, st);
    // Gemmi❗❌:   fill_residue_entity_type(st);
    // Gemmi❗❌:   st.setup_cell_images();
    // Gemmi❗❌:
    // Gemmi❗❌:   st.helices = read_helices(block);
    // Gemmi❗❌:   st.sheets = read_sheets(block);
    // Gemmi❗❌:   read_connectivity(block, st);
    // Gemmi❗❌:   read_prot_cis(block, st);
    // Gemmi❗❌:   read_struct_mod_residue(block, st);
    // Gemmi❗❌:   st.assemblies = read_assemblies(block);
    // Gemmi❗❌:   read_sifts_unp(block, st);
    // Gemmi❗❌:
    // Gemmi❗❌:   cif::Table chem_comp_table = block.find("_chem_comp.", {"id", "three_letter_code"});
    // Gemmi❗❌:   if (chem_comp_table.ok()) {
    // Gemmi❗❌:     for (auto row : chem_comp_table) {
    // Gemmi❗❌:       std::string alias = row.str(0);
    // Gemmi❗❌:       std::string long_id = row.str(1);
    // Gemmi❗❌:       if (alias[0] == '~' && long_id[0] != '~' && long_id[0] != '\0')
    // Gemmi❗❌:         st.shortened_ccd_codes.emplace_back(long_id, alias);
    // Gemmi❗❌:     }
    // Gemmi❗❌:     restore_full_ccd_codes(st);
    // Gemmi❗❌:   }
    // Gemmi❗❌: }
    // Behavior review: preserve the complete source call order, including
    // NCS before fractional transform, ORIGX afterward, SIFTS after assembly
    // creation, and CCD restoration last. The text wrapper uses the default
    // CIF check level; this shared path uses `make_structure` block selection
    // and never enters the separate chemcomp-coordinate reader. Each structured failure
    // aborts this local staging path; no partially built BIO value escapes.
    // Complexity review: each source category is selected/converted once.
    // The detached grouping and parsed-atom vectors retain additional O(A)
    // state versus Gemmi's direct hierarchy writes, and source-address maps
    // retain their existing ordered linear scans; therefore the cost axis is
    // intentionally marked worse until complete-reader profiling/review.
    let block = select_coordinate_block(document.blocks(), document.source())
        .map_err(|error| BioMmcifReadError::new(BioMmcifReadStage::CoordinateBlock, error))?;

    let mut source_state = BioStructureSourceState::default();
    source_state.name = block.name().to_owned();
    let mut metadata = BioMetadata::default();
    let mut crystal = parse_mmcif_crystal_cell_info(block)
        .map_err(|error| BioMmcifReadError::new(BioMmcifReadStage::CrystalCell, error))?;

    parse_mmcif_entry_info(block, &mut source_state);
    append_mmcif_audit_authors(block, &mut metadata);
    parse_mmcif_refinement_info(block, &mut source_state, &mut metadata)
        .map_err(|error| BioMmcifReadError::new(BioMmcifReadStage::Refinement, error))?;
    parse_mmcif_tls_info(block, &mut metadata)
        .map_err(|error| BioMmcifReadError::new(BioMmcifReadStage::Tls, error))?;
    parse_mmcif_experimental_info(block, &mut metadata)
        .map_err(|error| BioMmcifReadError::new(BioMmcifReadStage::Experimental, error))?;
    parse_mmcif_reflns_info(block, &mut metadata)
        .map_err(|error| BioMmcifReadError::new(BioMmcifReadStage::Reflections, error))?;
    parse_mmcif_software_info(block, &mut metadata)
        .map_err(|error| BioMmcifReadError::new(BioMmcifReadStage::Software, error))?;
    let ncs_operators = parse_mmcif_ncs_operators(block, &mut source_state)
        .map_err(|error| BioMmcifReadError::new(BioMmcifReadStage::Ncs, error))?;
    apply_mmcif_fractional_transform(block, &mut crystal)
        .map_err(|error| BioMmcifReadError::new(BioMmcifReadStage::FractionalTransform, error))?;
    parse_mmcif_origx(block, &mut source_state)
        .map_err(|error| BioMmcifReadError::new(BioMmcifReadStage::Origx, error))?;

    let anisotropic = extract_mmcif_anisotropic_u(block)
        .map_err(|error| BioMmcifReadError::new(BioMmcifReadStage::AnisotropicU, error))?;
    const ATOM_SITE_TAGS: [&str; 24] = [
        "id",
        "?group_PDB",
        "type_symbol",
        "?label_atom_id",
        "label_alt_id",
        "?label_comp_id",
        "label_asym_id",
        "?label_entity_id",
        "?label_seq_id",
        "?pdbx_PDB_ins_code",
        "Cartn_x",
        "Cartn_y",
        "Cartn_z",
        "?occupancy",
        "?B_iso_or_equiv",
        "?pdbx_formal_charge",
        "?auth_seq_id",
        "?auth_comp_id",
        "?auth_asym_id",
        "?auth_atom_id",
        "?pdbx_PDB_model_num",
        "?calc_flag",
        "?pdbx_tls_group_id",
        "?ccp4_deuterium_fraction",
    ];
    let atom_table = block
        .find("_atom_site.", &ATOM_SITE_TAGS)
        .map_err(|error| BioMmcifReadError::new(BioMmcifReadStage::AtomSites, error))?;
    source_state.has_d_fraction = has_mmcif_deuterium_fraction_column(&atom_table);
    let (mut grouping, parsed_atoms) =
        group_mmcif_atom_sites_with(&atom_table, |row, identities| {
            let name_text = cif_as_string(identities.atom_id.get(row));
            let name = AtomName::from_ascii(name_text.as_bytes()).ok_or_else(|| {
                MmcifAtomSiteGroupingError::AtomNameNotRepresentable {
                    value: name_text.clone(),
                }
            })?;
            let scalars = project_mmcif_atom_site_scalars(row, source_state.has_d_fraction)
                .map_err(MmcifAtomSiteGroupingError::AtomScalar)?;
            Ok(MmcifParsedAtomSite {
                name,
                raw_id: atom_site_raw_value(row, AtomSiteColumn::Id.index()).to_owned(),
                scalars,
            })
        })
        .map_err(|error| BioMmcifReadError::new(BioMmcifReadStage::AtomSites, error))?;

    let mut entities = parse_mmcif_entity_build_state(block, Some(&grouping))
        .map_err(|error| BioMmcifReadError::new(BioMmcifReadStage::EntitySequence, error))?;
    fill_mmcif_residue_entity_type(&mut grouping, &entities);
    setup_cell_images(&mut crystal, &ncs_operators);

    let mut helices = parse_mmcif_helices(block)
        .map_err(|error| BioMmcifReadError::new(BioMmcifReadStage::Helices, error))?;
    let mut sheets = parse_mmcif_sheets(block)
        .map_err(|error| BioMmcifReadError::new(BioMmcifReadStage::Sheets, error))?;
    let mut connections = parse_mmcif_connections(block, &grouping)
        .map_err(|error| BioMmcifReadError::new(BioMmcifReadStage::Connections, error))?;
    let mut cispeps = parse_mmcif_cis_peptides(block)
        .map_err(|error| BioMmcifReadError::new(BioMmcifReadStage::CisPeptides, error))?;
    let mut mod_residues = parse_mmcif_modified_residues(block)
        .map_err(|error| BioMmcifReadError::new(BioMmcifReadStage::ModifiedResidues, error))?;
    let assemblies = parse_mmcif_assemblies(block)
        .map_err(|error| BioMmcifReadError::new(BioMmcifReadStage::Assemblies, error))?;
    parse_mmcif_sifts_unp(block, &mut entities, &mut grouping.models)
        .map_err(|error| BioMmcifReadError::new(BioMmcifReadStage::SiftsUnp, error))?;

    let ccd_aliases = parse_mmcif_shortened_ccd_codes(block)
        .map_err(|error| BioMmcifReadError::new(BioMmcifReadStage::CcdRestoration, error))?;
    restore_mmcif_ccd_codes(
        &ccd_aliases,
        &mut grouping,
        &mut entities,
        &mut connections,
        &mut cispeps,
        &mut helices,
        &mut sheets,
        &mut mod_residues,
    )
    .map_err(|error| BioMmcifReadError::new(BioMmcifReadStage::CcdRestoration, error))?;

    let parts = materialize_mmcif_parts(
        source_state,
        metadata,
        crystal,
        ncs_operators,
        assemblies,
        entities,
        grouping,
        parsed_atoms,
        anisotropic,
        connections,
        cispeps,
        mod_residues,
        helices,
        sheets,
    )?;
    BioStructureData::from_parts(parts)
        .map_err(|error| BioMmcifReadError::new(BioMmcifReadStage::StructureValidation, error))
}

fn materialize_mmcif_parts(
    source_state: BioStructureSourceState,
    metadata: BioMetadata,
    crystal: BioCrystalInfo,
    ncs_operators: Vec<BioNcsOperator>,
    assemblies: Vec<BioAssembly>,
    entities: Vec<MmcifEntityBuildState>,
    grouping: MmcifAtomSiteGrouping,
    parsed_atoms: Vec<MmcifParsedAtomSite>,
    anisotropic: HashMap<String, [f32; 6]>,
    connections: Vec<BioConnection>,
    cispeps: Vec<BioCisPep>,
    mod_residues: Vec<BioModRes>,
    helices: Vec<BioHelix>,
    sheets: Vec<BioSheet>,
) -> Result<BioStructureParts, BioMmcifReadError> {
    let mut models = Vec::with_capacity(grouping.models.len());
    let chain_count = grouping.models.iter().map(|model| model.chains.len()).sum();
    let residue_count = grouping
        .models
        .iter()
        .flat_map(|model| &model.chains)
        .map(|chain| chain.residues.len())
        .sum();
    let atom_count = grouping
        .models
        .iter()
        .flat_map(|model| &model.chains)
        .flat_map(|chain| &chain.residues)
        .map(|residue| residue.atom_site_rows.len())
        .sum();
    let mut chains = Vec::with_capacity(chain_count);
    let mut residues = Vec::with_capacity(residue_count);
    let mut atoms = Vec::with_capacity(atom_count);
    let mut positions = Vec::with_capacity(atom_count);

    for model in &grouping.models {
        let model_id_value = u32::try_from(models.len()).map_err(|_| {
            BioMmcifReadError::new(
                BioMmcifReadStage::Materialization,
                BioStructureError::RowIndexTooLarge {
                    value: models.len(),
                },
            )
        })?;
        let model_id = BioModelId::new(model_id_value);
        let chain_start = chains.len();
        for chain in &model.chains {
            let chain_id_value = u32::try_from(chains.len()).map_err(|_| {
                BioMmcifReadError::new(
                    BioMmcifReadStage::Materialization,
                    BioStructureError::RowIndexTooLarge {
                        value: chains.len(),
                    },
                )
            })?;
            let chain_id = BioChainId::new(chain_id_value);
            let residue_start = residues.len();
            let chain_entity_id = chain
                .source
                .label_asym_id()
                .filter(|subchain| !subchain.is_empty())
                .and_then(|subchain| {
                    entities.iter().position(|entity| {
                        entity
                            .subchains
                            .iter()
                            .any(|candidate| candidate == subchain)
                    })
                })
                .map(|index| {
                    u32::try_from(index).map(BioEntityId::new).map_err(|_| {
                        BioMmcifReadError::new(
                            BioMmcifReadStage::Materialization,
                            BioStructureError::RowIndexTooLarge { value: index },
                        )
                    })
                })
                .transpose()?;

            for residue in &chain.residues {
                let residue_id_value = u32::try_from(residues.len()).map_err(|_| {
                    BioMmcifReadError::new(
                        BioMmcifReadStage::Materialization,
                        BioStructureError::RowIndexTooLarge {
                            value: residues.len(),
                        },
                    )
                })?;
                let residue_id = BioResidueId::new(residue_id_value);
                let atom_start = atoms.len();
                for &source_row_index in &residue.atom_site_rows {
                    let atom = parsed_atoms.get(source_row_index).ok_or_else(|| {
                        BioMmcifReadError::new(
                            BioMmcifReadStage::Materialization,
                            MmcifAtomSiteRowOutOfBounds {
                                row_index: source_row_index,
                                parsed_count: parsed_atoms.len(),
                            },
                        )
                    })?;
                    let anisou = anisotropic
                        .get(&atom.raw_id)
                        .copied()
                        .unwrap_or([0.0; 6])
                        .map(f64::from);
                    atoms.push(BioAtomRow::new(
                        residue_id,
                        atom.name,
                        atom.scalars.element,
                        atom.scalars.isotope_mass_number,
                        atom.scalars.altloc,
                        atom.scalars.formal_charge,
                        atom.scalars.calc_flag,
                        f64::from(atom.scalars.occupancy),
                        f64::from(atom.scalars.b_iso),
                        anisou,
                        atom.scalars.tls_group_id,
                        f64::from(atom.scalars.fraction),
                        AtomSourceIds::new(Some(atom.scalars.serial)),
                    ));
                    positions.push(atom.scalars.position);
                }
                let atom_span = BioRowSpan::from_usize(atom_start, atoms.len() - atom_start)
                    .map_err(|error| {
                        BioMmcifReadError::new(BioMmcifReadStage::Materialization, error)
                    })?;
                let entity_id = residue
                    .source
                    .label_entity_id()
                    .and_then(|source_id| {
                        entities
                            .iter()
                            .position(|entity| entity.source_entity_id == source_id)
                    })
                    .map(|index| {
                        u32::try_from(index).map(BioEntityId::new).map_err(|_| {
                            BioMmcifReadError::new(
                                BioMmcifReadStage::Materialization,
                                BioStructureError::RowIndexTooLarge { value: index },
                            )
                        })
                    })
                    .transpose()?;
                let residue_info_kind = find_residue_info(residue.address.name().as_str()).kind;
                residues.push(BioResidueRow::new(
                    chain_id,
                    atom_span,
                    residue.address.name(),
                    residue_info_kind,
                    residue.entity_kind,
                    entity_id,
                    residue.het_flag,
                    residue.source.clone(),
                    residue.sifts_unp,
                ));
            }
            let residue_span =
                BioRowSpan::from_usize(residue_start, residues.len() - residue_start).map_err(
                    |error| BioMmcifReadError::new(BioMmcifReadStage::Materialization, error),
                )?;
            chains.push(BioChainRow::new(
                model_id,
                chain_entity_id,
                residue_span,
                ChainKind::Unknown,
                chain.source.clone(),
            ));
        }
        let chain_span = BioRowSpan::from_usize(chain_start, chains.len() - chain_start)
            .map_err(|error| BioMmcifReadError::new(BioMmcifReadStage::Materialization, error))?;
        models.push(BioModelRow::new(
            chain_span,
            Some(model.source_model_number),
        ));
    }

    let entities = mmcif_entity_rows(entities);
    Ok(BioStructureParts {
        input_format: BioCoordinateFormat::Mmcif,
        models,
        chains,
        residues,
        atoms,
        entities,
        connections,
        cispeps,
        mod_residues,
        helices,
        sheets,
        metadata,
        source_state,
        coordinates: BioCoordinateBlock::new(positions),
        crystal: Some(crystal),
        ncs_operators,
        assemblies,
    })
}

#[cfg(test)]
mod tests {
    use super::{
        AtomSiteColumn, AtomSiteIdentityField, CoordinateBlockSelectionError,
        MissingAtomSiteIdentityField, MmcifAssemblyError, MmcifAtomSiteGrouping,
        MmcifAtomSiteScalarError, MmcifCisPepError, MmcifConnectionError, MmcifCrystalInfoError,
        MmcifEntityInfoError, MmcifHelixError, MmcifLabelAddressError, MmcifModResError,
        MmcifOperationExprError, MmcifSequenceIdError, MmcifSheetError, MmcifSiftsUnpError,
        MmcifTlsInfoError, append_mmcif_audit_authors, atom_site_raw_value,
        bind_atom_site_identity_columns, extract_mmcif_anisotropic_u, group_mmcif_atom_sites,
        has_mmcif_deuterium_fraction_column, make_mmcif_residue_address, make_mmcif_seq_id,
        parse_mmcif_assemblies, parse_mmcif_cis_peptides, parse_mmcif_connections,
        parse_mmcif_crystal_info, parse_mmcif_entity_polymer_values, parse_mmcif_entry_info,
        parse_mmcif_experimental_info, parse_mmcif_helices, parse_mmcif_modified_residues,
        parse_mmcif_ncs_origx_info, parse_mmcif_operation_expr, parse_mmcif_refinement_info,
        parse_mmcif_reflns_info, parse_mmcif_sheets, parse_mmcif_software_info,
        parse_mmcif_tls_info, populate_mmcif_bio_structure_document,
        project_mmcif_atom_site_scalars, read_mmcif_bio_structure, resolve_mmcif_label_address,
        select_coordinate_block,
    };
    use crate::cif::{CifBlock, CifCheckLevel, CifReadErrorKind, CifTable, read_cif_document};
    use cosmolkit_bio::{
        AltLocLabel, BioAssemblySpecialKind, BioAsu, BioCalcFlag, BioConnectionKind,
        BioCrystalCell, BioDiffractionInfo, BioEntityDbRef, BioEntityRow, BioExperimentInfo,
        BioExperimentalCrystalInfo, BioHelixClass, BioMetadata, BioRefinementInfo,
        BioReflectionsInfo, BioSiftsUnpResidue, BioSoftwareClassification, BioSoftwareItem,
        BioStructureError, BioStructureSourceState, BioTlsGroup, BioTlsSelection, BioTransform,
        EntityKind, PdbChainId, PdbSeqId, PolymerKind, ResidueName,
    };
    use cosmolkit_types::Element;
    use std::collections::BTreeMap;

    const ATOM_SITE_TAGS: [&str; 24] = [
        "id",
        "?group_PDB",
        "type_symbol",
        "?label_atom_id",
        "label_alt_id",
        "?label_comp_id",
        "label_asym_id",
        "?label_entity_id",
        "?label_seq_id",
        "?pdbx_PDB_ins_code",
        "Cartn_x",
        "Cartn_y",
        "Cartn_z",
        "?occupancy",
        "?B_iso_or_equiv",
        "?pdbx_formal_charge",
        "?auth_seq_id",
        "?auth_comp_id",
        "?auth_asym_id",
        "?auth_atom_id",
        "?pdbx_PDB_model_num",
        "?calc_flag",
        "?pdbx_tls_group_id",
        "?ccp4_deuterium_fraction",
    ];

    const STRUCT_CONN_TAGS: [&str; 22] = [
        "id",
        "conn_type_id",
        "?ptnr1_auth_asym_id",
        "?ptnr2_auth_asym_id",
        "?ptnr1_label_asym_id",
        "?ptnr2_label_asym_id",
        "ptnr1_label_comp_id",
        "ptnr2_label_comp_id",
        "ptnr1_label_atom_id",
        "ptnr2_label_atom_id",
        "?pdbx_ptnr1_label_alt_id",
        "?pdbx_ptnr2_label_alt_id",
        "?ptnr1_auth_seq_id",
        "?ptnr2_auth_seq_id",
        "?ptnr1_label_seq_id",
        "?ptnr2_label_seq_id",
        "?pdbx_ptnr1_PDB_ins_code",
        "?pdbx_ptnr2_PDB_ins_code",
        "?ptnr1_symmetry",
        "?ptnr2_symmetry",
        "?pdbx_dist_value",
        "?ccp4_link_id",
    ];

    const CIS_PEPTIDE_TAGS: [&str; 13] = [
        "pdbx_PDB_model_num",
        "auth_asym_id",
        "auth_seq_id",
        "?pdbx_PDB_ins_code",
        "?label_comp_id",
        "?auth_comp_id",
        "?pdbx_auth_asym_id_2",
        "?pdbx_auth_seq_id_2",
        "?pdbx_PDB_ins_code_2",
        "?pdbx_label_comp_id_2",
        "?pdbx_auth_comp_id_2",
        "?label_alt_id",
        "?pdbx_omega_angle",
    ];

    const MODIFIED_RESIDUE_TAGS: [&str; 8] = [
        "auth_asym_id",
        "auth_seq_id",
        "?PDB_ins_code",
        "?auth_comp_id",
        "?label_comp_id",
        "?parent_comp_id",
        "?details",
        "?ccp4_mod_id",
    ];

    fn document(text: &str, source: &str) -> crate::cif::CifDocument {
        read_cif_document(text, source, CifCheckLevel::Syntax).expect("fixed C01 CIF document")
    }

    #[test]
    fn bio_mmcif_c13_cell_uses_exact_keys_and_first_three_length_gate() {
        let legacy_aliases = document(
            concat!(
                "data_c13_legacy\n",
                "_cell_length_a 10\n_cell_length_b 20\n_cell_length_c 30\n",
                "_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n",
            ),
            "c13-legacy-cell.cif",
        );
        let crystal = parse_mmcif_crystal_info(&legacy_aliases.blocks()[0]).unwrap();
        assert_eq!(crystal.cell(), BioCrystalCell::default());
        assert_eq!(crystal.volume(), 1.0);

        let null_length = document(
            concat!(
                "data_c13_null\n",
                "_cell.length_a .\n_cell.length_b 20\n_cell.length_c 30\n",
                "_cell.angle_alpha 180\n_cell.angle_beta 90\n_cell.angle_gamma 90\n",
            ),
            "c13-null-length.cif",
        );
        let crystal = parse_mmcif_crystal_info(&null_length.blocks()[0]).unwrap();
        assert_eq!(crystal.cell(), BioCrystalCell::default());
        assert_eq!(crystal.volume(), 1.0);

        let zero_gamma = document(
            concat!(
                "data_c13_zero_gamma\n",
                "_cell.length_a 10\n_cell.length_b 20\n_cell.length_c 30\n",
                "_cell.angle_alpha 90\n_cell.angle_beta 90\n_cell.angle_gamma 0\n",
            ),
            "c13-zero-gamma.cif",
        );
        let crystal = parse_mmcif_crystal_info(&zero_gamma.blocks()[0]).unwrap();
        assert_eq!(crystal.cell(), BioCrystalCell::default());
        assert_eq!(crystal.volume(), 1.0);
    }

    #[test]
    fn bio_mmcif_c13_cell_pair_loop_and_source_errors_follow_gemmi() {
        let pair = document(
            concat!(
                "data_c13_pair\n",
                "_cell.length_a 2\n_cell.length_b 3\n_cell.length_c 4\n",
                "_cell.angle_alpha 90\n_cell.angle_beta 90\n_cell.angle_gamma 90\n",
            ),
            "c13-cell-pair.cif",
        );
        let crystal = parse_mmcif_crystal_info(&pair.blocks()[0]).unwrap();
        assert_eq!(
            crystal.cell(),
            BioCrystalCell {
                a: 2.0,
                b: 3.0,
                c: 4.0,
                alpha: 90.0,
                beta: 90.0,
                gamma: 90.0,
            }
        );
        assert_eq!(crystal.volume(), 24.0);

        let looped = document(
            concat!(
                "data_c13_loop\nloop_\n",
                "_cell.length_a\n_cell.length_b\n_cell.length_c\n",
                "_cell.angle_alpha\n_cell.angle_beta\n_cell.angle_gamma\n",
                "2 3 4 90 90 90\n",
            ),
            "c13-cell-loop.cif",
        );
        let crystal = parse_mmcif_crystal_info(&looped.blocks()[0]).unwrap();
        assert_eq!(crystal.cell().a, 2.0);
        assert_eq!(crystal.cell().b, 3.0);
        assert_eq!(crystal.cell().c, 4.0);

        let multiple_rows = document(
            concat!(
                "data_c13_multi\nloop_\n",
                "_cell.length_a\n_cell.length_b\n_cell.length_c\n",
                "_cell.angle_alpha\n_cell.angle_beta\n_cell.angle_gamma\n",
                "2 3 4 90 90 90\n5 6 7 90 90 90\n",
            ),
            "c13-cell-multiple-rows.cif",
        );
        let error = parse_mmcif_crystal_info(&multiple_rows.blocks()[0]).unwrap_err();
        let MmcifCrystalInfoError::Cif(error) = error else {
            panic!("Gemmi Table::one failure remains a CIF error");
        };
        assert_eq!(error.kind(), CifReadErrorKind::InvalidValue);

        let impossible_angle = document(
            concat!(
                "data_c13_impossible\n",
                "_cell.length_a 2\n_cell.length_b 3\n_cell.length_c 4\n",
                "_cell.angle_alpha 0\n_cell.angle_beta 90\n_cell.angle_gamma 90\n",
            ),
            "c13-impossible-angle.cif",
        );
        assert!(matches!(
            parse_mmcif_crystal_info(&impossible_angle.blocks()[0]),
            Err(MmcifCrystalInfoError::Bio(
                BioStructureError::ImpossibleCrystalAngle
            ))
        ));

        let one_hundred_eighty = document(
            concat!(
                "data_c13_180_degree_angle\n",
                "_cell.length_a 2\n_cell.length_b 3\n_cell.length_c 4\n",
                "_cell.angle_alpha 180\n_cell.angle_beta 90\n_cell.angle_gamma 90\n",
            ),
            "c13-180-degree-angle.cif",
        );
        let crystal = parse_mmcif_crystal_info(&one_hundred_eighty.blocks()[0]).unwrap();
        assert_eq!(crystal.cell().alpha, 180.0);
        assert_eq!(crystal.volume(), 0.0);
    }

    #[test]
    fn bio_mmcif_c13_space_group_uses_only_exact_source_key_and_precedence() {
        let exact_pair = document(
            concat!(
                "data_c13_hm\n",
                "_symmetry.space_group_name_H-M 'P 1'\n",
                "_space_group.name_H-M 'P 2'\n",
                "_symmetry.space_group_name_H-M_alt 'P 3'\n",
            ),
            "c13-space-group-pair.cif",
        );
        let crystal = parse_mmcif_crystal_info(&exact_pair.blocks()[0]).unwrap();
        assert_eq!(crystal.space_group_hm(), Some("P 1"));

        let aliases_only = document(
            "data_c13_aliases\n_space_group.name_H-M 'P 2'\n_symmetry.space_group_name_H-M_alt 'P 3'\n",
            "c13-space-group-aliases.cif",
        );
        let crystal = parse_mmcif_crystal_info(&aliases_only.blocks()[0]).unwrap();
        assert_eq!(crystal.space_group_hm(), Some(""));

        let multiple_loop_rows = document(
            concat!(
                "data_c13_hm_loop\nloop_\n_symmetry.space_group_name_H-M\n",
                "'P 1'\n'P 2'\n",
            ),
            "c13-space-group-multiple-loop-rows.cif",
        );
        let crystal = parse_mmcif_crystal_info(&multiple_loop_rows.blocks()[0]).unwrap();
        assert_eq!(crystal.space_group_hm(), Some(""));
    }

    #[test]
    fn bio_mmcif_c13_fractional_transform_uses_first_row_and_ignores_pdb_scale() {
        let doc = document(
            concat!(
                "data_c13_fract\n",
                "_cell.length_a 10\n_cell.length_b 20\n_cell.length_c 30\n",
                "_cell.angle_alpha 90\n_cell.angle_beta 90\n_cell.angle_gamma 90\n",
                "_database_PDB_matrix.scale[1][1] 9\n",
                "loop_\n",
                "_atom_sites.fract_transf_matrix[1][1]\n",
                "_atom_sites.fract_transf_matrix[1][2]\n",
                "_atom_sites.fract_transf_matrix[1][3]\n",
                "_atom_sites.fract_transf_vector[1]\n",
                "_atom_sites.fract_transf_matrix[2][1]\n",
                "_atom_sites.fract_transf_matrix[2][2]\n",
                "_atom_sites.fract_transf_matrix[2][3]\n",
                "_atom_sites.fract_transf_vector[2]\n",
                "_atom_sites.fract_transf_matrix[3][1]\n",
                "_atom_sites.fract_transf_matrix[3][2]\n",
                "_atom_sites.fract_transf_matrix[3][3]\n",
                "_atom_sites.fract_transf_vector[3]\n",
                "0.5 0 0 1 0 0.25 0 2 0 0 0.125 3\n",
                "0.75 0 0 4 0 0.5 0 5 0 0 0.25 6\n",
            ),
            "c13-fract.cif",
        );
        let crystal = parse_mmcif_crystal_info(&doc.blocks()[0]).unwrap();
        assert_eq!(
            *crystal.fractional().matrix(),
            [[0.5, 0.0, 0.0], [0.0, 0.25, 0.0], [0.0, 0.0, 0.125]]
        );
        assert_eq!(*crystal.fractional().translation(), [1.0, 2.0, 3.0]);
        assert!(crystal.explicit_matrices());
    }

    #[test]
    fn bio_mmcif_c14_ncs_preserves_source_matrix_order_identity_and_nan_branches() {
        let doc = document(
            concat!(
                "data_c14_ncs\nloop_\n",
                "_struct_ncs_oper.matrix[1][1]\n",
                "_struct_ncs_oper.matrix[1][2]\n",
                "_struct_ncs_oper.matrix[1][3]\n",
                "_struct_ncs_oper.vector[1]\n",
                "_struct_ncs_oper.matrix[2][1]\n",
                "_struct_ncs_oper.matrix[2][2]\n",
                "_struct_ncs_oper.matrix[2][3]\n",
                "_struct_ncs_oper.vector[2]\n",
                "_struct_ncs_oper.matrix[3][1]\n",
                "_struct_ncs_oper.matrix[3][2]\n",
                "_struct_ncs_oper.matrix[3][3]\n",
                "_struct_ncs_oper.vector[3]\n",
                "_struct_ncs_oper.id\n",
                "_struct_ncs_oper.code\n",
                "1 0 0 0 0 1 0 0 0 0 1 0 identity_first given\n",
                "0.5 0.25 -0.5 10 1.5 2.5 3.5 11 4.5 5.5 6.5 12 op_first given\n",
                ". 0 0 0 0 1 0 0 0 0 1 0 nan_row given\n",
                "1 2 3 4 5 6 7 8 9 10 11 12 op_second Given\n",
                "1 0 0 0 0 1 0 0 0 0 1 0 identity_last given\n",
            ),
            "c14-ncs.cif",
        );
        let projected = parse_mmcif_ncs_origx_info(&doc.blocks()[0]).unwrap();

        assert_eq!(projected.ncs_operators.len(), 2);
        let first = &projected.ncs_operators[0];
        assert_eq!(first.id, "op_first");
        assert!(first.given);
        assert_eq!(
            *first.transform.matrix(),
            [[0.5, 0.25, -0.5], [1.5, 2.5, 3.5], [4.5, 5.5, 6.5]]
        );
        assert_eq!(*first.transform.translation(), [10.0, 11.0, 12.0]);
        let second = &projected.ncs_operators[1];
        assert_eq!(second.id, "op_second");
        assert!(!second.given);
        assert_eq!(
            *second.transform.matrix(),
            [[1.0, 2.0, 3.0], [5.0, 6.0, 7.0], [9.0, 10.0, 11.0]]
        );
        assert_eq!(*second.transform.translation(), [4.0, 8.0, 12.0]);
        assert_eq!(
            projected
                .source_state
                .info
                .get("_struct_ncs_oper.id")
                .map(String::as_str),
            Some("identity_last")
        );
        assert!(!projected.source_state.has_origx);
        assert_eq!(projected.source_state.origx, BioTransform::identity());

        let no_code_column = document(
            concat!(
                "data_c14_ncs_without_code\nloop_\n",
                "_struct_ncs_oper.matrix[1][1]\n",
                "_struct_ncs_oper.matrix[1][2]\n",
                "_struct_ncs_oper.matrix[1][3]\n",
                "_struct_ncs_oper.vector[1]\n",
                "_struct_ncs_oper.matrix[2][1]\n",
                "_struct_ncs_oper.matrix[2][2]\n",
                "_struct_ncs_oper.matrix[2][3]\n",
                "_struct_ncs_oper.vector[2]\n",
                "_struct_ncs_oper.matrix[3][1]\n",
                "_struct_ncs_oper.matrix[3][2]\n",
                "_struct_ncs_oper.matrix[3][3]\n",
                "_struct_ncs_oper.vector[3]\n",
                "_struct_ncs_oper.id\n",
                "0.5 0 0 0 0 1 0 0 0 0 1 0 no_code\n",
            ),
            "c14-ncs-no-code.cif",
        );
        let projected = parse_mmcif_ncs_origx_info(&no_code_column.blocks()[0]).unwrap();
        assert_eq!(projected.ncs_operators.len(), 1);
        assert!(!projected.ncs_operators[0].given);
    }

    #[test]
    fn bio_legacy_n02_identity_id_last_wins_and_nonidentity_order_is_stable() {
        const HEADER: &str = concat!(
            "data_ncs\nloop_\n",
            "_struct_ncs_oper.matrix[1][1]\n_struct_ncs_oper.matrix[1][2]\n_struct_ncs_oper.matrix[1][3]\n_struct_ncs_oper.vector[1]\n",
            "_struct_ncs_oper.matrix[2][1]\n_struct_ncs_oper.matrix[2][2]\n_struct_ncs_oper.matrix[2][3]\n_struct_ncs_oper.vector[2]\n",
            "_struct_ncs_oper.matrix[3][1]\n_struct_ncs_oper.matrix[3][2]\n_struct_ncs_oper.matrix[3][3]\n_struct_ncs_oper.vector[3]\n",
            "_struct_ncs_oper.id\n_struct_ncs_oper.code\n",
        );
        const ID_A: &str = "1 0 0 0 0 1 0 0 0 0 1 0 id_a given\n";
        const ID_B: &str = "1 0 0 0 0 1 0 0 0 0 1 0 id_b generated\n";
        const OP_A: &str = "1 0 0 2 0 1 0 0 0 0 1 0 op_a given\n";
        const OP_B: &str = "1 0 0 3 0 1 0 0 0 0 1 0 op_b generated\n";
        const INCOMPLETE: &str = ". 0 0 0 0 1 0 0 0 0 1 0 missing given\n";
        for (rows, identity, operators) in [
            (OP_A.to_owned(), None, &["op_a"][..]),
            (format!("{ID_A}{OP_A}"), Some("id_a"), &["op_a"]),
            (
                format!("{OP_A}{ID_A}{INCOMPLETE}{ID_B}{OP_B}"),
                Some("id_b"),
                &["op_a", "op_b"],
            ),
        ] {
            let content = format!("{HEADER}{rows}");
            let doc = document(&content, "n02.cif");
            let projected = parse_mmcif_ncs_origx_info(&doc.blocks()[0]).unwrap();
            assert_eq!(
                projected
                    .source_state
                    .info
                    .get("_struct_ncs_oper.id")
                    .map(String::as_str),
                identity
            );
            assert_eq!(
                projected
                    .ncs_operators
                    .iter()
                    .map(|op| op.id.as_str())
                    .collect::<Vec<_>>(),
                operators
            );
            assert!(projected.ncs_operators[0].given);
            if projected.ncs_operators.len() > 1 {
                assert!(!projected.ncs_operators[1].given);
            }
        }
    }

    #[test]
    fn bio_mmcif_c14_origx_uses_presence_gate_and_first_transform_row() {
        let absent = document("data_c14_origx_absent\n", "c14-origx-absent.cif");
        let projected = parse_mmcif_ncs_origx_info(&absent.blocks()[0]).unwrap();
        assert!(!projected.source_state.has_origx);
        assert_eq!(projected.source_state.origx, BioTransform::identity());

        let incomplete = document(
            concat!(
                "data_c14_origx_incomplete\n",
                "_database_PDB_matrix.origx[1][1] 0.5\n",
                "_database_PDB_matrix.origx[1][2] 0\n",
                "_database_PDB_matrix.origx[1][3] 0\n",
                "_database_PDB_matrix.origx_vector[1] 10\n",
                "_database_PDB_matrix.origx[2][1] 0\n",
                "_database_PDB_matrix.origx[2][2] 0.25\n",
                "_database_PDB_matrix.origx[2][3] 0\n",
                "_database_PDB_matrix.origx_vector[2] 20\n",
                "_database_PDB_matrix.origx[3][1] 0\n",
                "_database_PDB_matrix.origx[3][2] 0\n",
                "_database_PDB_matrix.origx[3][3] 0.125\n",
            ),
            "c14-origx-incomplete.cif",
        );
        let projected = parse_mmcif_ncs_origx_info(&incomplete.blocks()[0]).unwrap();
        assert!(!projected.source_state.has_origx);
        assert_eq!(projected.source_state.origx, BioTransform::identity());

        let multiple_rows = document(
            concat!(
                "data_c14_origx_rows\nloop_\n",
                "_database_PDB_matrix.origx[1][1]\n",
                "_database_PDB_matrix.origx[1][2]\n",
                "_database_PDB_matrix.origx[1][3]\n",
                "_database_PDB_matrix.origx_vector[1]\n",
                "_database_PDB_matrix.origx[2][1]\n",
                "_database_PDB_matrix.origx[2][2]\n",
                "_database_PDB_matrix.origx[2][3]\n",
                "_database_PDB_matrix.origx_vector[2]\n",
                "_database_PDB_matrix.origx[3][1]\n",
                "_database_PDB_matrix.origx[3][2]\n",
                "_database_PDB_matrix.origx[3][3]\n",
                "_database_PDB_matrix.origx_vector[3]\n",
                "0.5 0 0 10 0 0.25 0 20 0 0 0.125 30\n",
                "0.75 0 0 40 0 0.5 0 50 0 0 0.25 60\n",
            ),
            "c14-origx-rows.cif",
        );
        let projected = parse_mmcif_ncs_origx_info(&multiple_rows.blocks()[0]).unwrap();
        assert!(projected.source_state.has_origx);
        assert_eq!(
            *projected.source_state.origx.matrix(),
            [[0.5, 0.0, 0.0], [0.0, 0.25, 0.0], [0.0, 0.0, 0.125]]
        );
        assert_eq!(
            *projected.source_state.origx.translation(),
            [10.0, 20.0, 30.0]
        );

        let identity = document(
            concat!(
                "data_c14_origx_identity\n",
                "_database_PDB_matrix.origx[1][1] 1\n",
                "_database_PDB_matrix.origx[1][2] 0\n",
                "_database_PDB_matrix.origx[1][3] 0\n",
                "_database_PDB_matrix.origx_vector[1] 0\n",
                "_database_PDB_matrix.origx[2][1] 0\n",
                "_database_PDB_matrix.origx[2][2] 1\n",
                "_database_PDB_matrix.origx[2][3] 0\n",
                "_database_PDB_matrix.origx_vector[2] 0\n",
                "_database_PDB_matrix.origx[3][1] 0\n",
                "_database_PDB_matrix.origx[3][2] 0\n",
                "_database_PDB_matrix.origx[3][3] 1\n",
                "_database_PDB_matrix.origx_vector[3] 0\n",
            ),
            "c14-origx-identity.cif",
        );
        let projected = parse_mmcif_ncs_origx_info(&identity.blocks()[0]).unwrap();
        assert!(projected.source_state.has_origx);
        assert_eq!(projected.source_state.origx, BioTransform::identity());
    }

    #[test]
    fn bio_mmcif_c15_info_decodes_and_aggregates_source_rows_in_key_order() {
        let doc = document(
            concat!(
                "data_c15_info_aggregate\n",
                "_entry.id \"entry 1\"\n",
                "loop_\n_cell.Z_PDB\n1\n2\n",
                "loop_\n_exptl.method\n'X-ray diffraction'\n?\nNeutron\n",
                "loop_\n_struct.title\n'First title'\n.\n'Second title'\n",
                "loop_\n_database_PDB_rev.date_original\n",
                "2020-01-02\n.\n2020-01-03\n",
                "_struct_keywords.pdbx_keywords kinase\n",
                "_struct_keywords.text 'line one; line two'\n",
            ),
            "c15-info-aggregate.cif",
        );
        let mut state = BioStructureSourceState::default();
        parse_mmcif_entry_info(&doc.blocks()[0], &mut state);

        assert_eq!(
            state.info,
            BTreeMap::from([
                ("_cell.Z_PDB".to_owned(), "1; 2".to_owned()),
                (
                    "_database_PDB_rev.date_original".to_owned(),
                    "2020-01-02; 2020-01-03".to_owned()
                ),
                ("_entry.id".to_owned(), "entry 1".to_owned()),
                (
                    "_exptl.method".to_owned(),
                    "X-ray diffraction; Neutron".to_owned()
                ),
                (
                    "_pdbx_database_status.recvd_initial_deposition_date".to_owned(),
                    "2020-01-02; 2020-01-03".to_owned()
                ),
                (
                    "_struct.title".to_owned(),
                    "First title; Second title".to_owned()
                ),
                (
                    "_struct_keywords.pdbx_keywords".to_owned(),
                    "kinase".to_owned()
                ),
                (
                    "_struct_keywords.text".to_owned(),
                    "line one; line two".to_owned()
                ),
            ])
        );
    }

    #[test]
    fn bio_mmcif_c15_info_uses_first_item_and_preserves_missing_or_all_null_keys() {
        let doc = document(
            concat!(
                "data_c15_info_presence\n",
                "_entry.id first_item\n",
                "loop_\n_entry.id\nsecond_item\nthird_item\n",
                "loop_\n_exptl.method\n.\n?\n",
            ),
            "c15-info-presence.cif",
        );
        let mut state = BioStructureSourceState::default();
        state
            .info
            .insert("_exptl.method".to_owned(), "retained".to_owned());
        state
            .info
            .insert("_struct.title".to_owned(), "preexisting".to_owned());
        parse_mmcif_entry_info(&doc.blocks()[0], &mut state);

        assert_eq!(
            state.info.get("_entry.id").map(String::as_str),
            Some("first_item")
        );
        assert_eq!(
            state.info.get("_exptl.method").map(String::as_str),
            Some("retained")
        );
        assert_eq!(
            state.info.get("_struct.title").map(String::as_str),
            Some("preexisting")
        );
        assert!(!state.info.contains_key("_cell.Z_PDB"));
    }

    #[test]
    fn bio_mmcif_c15_info_falls_back_to_old_date_only_when_new_key_is_absent() {
        let null_new_date = document(
            concat!(
                "data_c15_date_null\n",
                "_database_PDB_rev.date_original 2001-02-03\n",
                "loop_\n_pdbx_database_status.recvd_initial_deposition_date\n?\n.\n",
            ),
            "c15-date-null.cif",
        );
        let mut state = BioStructureSourceState::default();
        parse_mmcif_entry_info(&null_new_date.blocks()[0], &mut state);
        assert_eq!(
            state
                .info
                .get("_pdbx_database_status.recvd_initial_deposition_date")
                .map(String::as_str),
            Some("2001-02-03")
        );

        let explicit_new_date = document(
            concat!(
                "data_c15_date_explicit\n",
                "_database_PDB_rev.date_original 2001-02-03\n",
                "loop_\n_pdbx_database_status.recvd_initial_deposition_date\n",
                "2002-03-04\n2003-04-05\n",
            ),
            "c15-date-explicit.cif",
        );
        let mut state = BioStructureSourceState::default();
        parse_mmcif_entry_info(&explicit_new_date.blocks()[0], &mut state);
        assert_eq!(
            state
                .info
                .get("_pdbx_database_status.recvd_initial_deposition_date")
                .map(String::as_str),
            Some("2002-03-04; 2003-04-05")
        );

        let empty_new_date = document(
            "data_c15_date_empty\n_database_PDB_rev.date_original old\n_pdbx_database_status.recvd_initial_deposition_date ''\n",
            "c15-date-empty.cif",
        );
        let mut state = BioStructureSourceState::default();
        parse_mmcif_entry_info(&empty_new_date.blocks()[0], &mut state);
        assert_eq!(
            state
                .info
                .get("_pdbx_database_status.recvd_initial_deposition_date")
                .map(String::as_str),
            Some("")
        );

        let missing_new_date = document(
            "data_c15_date_existing\n_database_PDB_rev.date_original source_old\n",
            "c15-date-existing.cif",
        );
        let mut state = BioStructureSourceState::default();
        state.info.insert(
            "_pdbx_database_status.recvd_initial_deposition_date".to_owned(),
            "existing_new".to_owned(),
        );
        parse_mmcif_entry_info(&missing_new_date.blocks()[0], &mut state);
        assert_eq!(
            state
                .info
                .get("_pdbx_database_status.recvd_initial_deposition_date")
                .map(String::as_str),
            Some("existing_new")
        );
    }

    #[test]
    fn bio_mmcif_c16_authors_preserve_order_nulls_decoding_duplicates_and_existing_values() {
        let doc = document(
            concat!(
                "data_c16_author_order\n",
                "loop_\n_audit_author.id\n_audit_author.name\n",
                "1 'Ada Lovelace'\n",
                "2 .\n",
                "3 ?\n",
                "4 'Ada Lovelace'\n",
                "5 ''\n",
            ),
            "c16-author-order.cif",
        );
        let mut metadata = BioMetadata::default();
        metadata.authors.push("preexisting".to_owned());

        append_mmcif_audit_authors(&doc.blocks()[0], &mut metadata);

        assert_eq!(
            metadata.authors,
            ["preexisting", "Ada Lovelace", "Ada Lovelace", ""]
        );
    }

    #[test]
    fn bio_mmcif_c16_authors_use_first_case_insensitive_matching_item() {
        let doc = document(
            concat!(
                "data_c16_author_first_item\n",
                "_AUDIT_AUTHOR.NAME 'first author'\n",
                "loop_\n_audit_author.name\nsecond_author\nthird_author\n",
            ),
            "c16-author-first-item.cif",
        );
        let mut metadata = BioMetadata::default();

        append_mmcif_audit_authors(&doc.blocks()[0], &mut metadata);

        assert_eq!(metadata.authors, ["first author"]);
    }

    #[test]
    fn bio_mmcif_c16_authors_leave_metadata_unchanged_when_missing_or_all_null() {
        let absent = document(
            "data_c16_author_absent\n_entry.id entry\n",
            "c16-author-absent.cif",
        );
        let all_null = document(
            "data_c16_author_null\nloop_\n_audit_author.name\n.\n?\n",
            "c16-author-null.cif",
        );
        let mut metadata = BioMetadata::default();
        metadata.authors.push("existing".to_owned());

        append_mmcif_audit_authors(&absent.blocks()[0], &mut metadata);
        append_mmcif_audit_authors(&all_null.blocks()[0], &mut metadata);

        assert_eq!(metadata.authors, ["existing"]);
    }

    #[test]
    fn bio_mmcif_c17_software_classification_mapping_and_source_row_order() {
        let doc = document(
            concat!(
                "data_c17_software_order\n",
                "loop_\n",
                "_software.name\n",
                "_software.classification\n",
                "_software.version\n",
                "_software.date\n",
                "_software.description\n",
                "_software.contact_author\n",
                "_software.contact_author_email\n",
                "tool0 'dAtA cOlLeCtIoN' v0 date0 desc0 author0 email0\n",
                "tool1 'DATA EXTRACTION' v1 date1 desc1 author1 email1\n",
                "tool2 'Data Processing' v2 date2 desc2 author2 email2\n",
                "tool3 'data reduction' v3 date3 desc3 author3 email3\n",
                "tool4 'DATA SCALING' v4 date4 desc4 author4 email4\n",
                "tool5 'model building' v5 date5 desc5 author5 email5\n",
                "tool6 phasing v6 date6 desc6 author6 email6\n",
                "tool7 Refinement v7 date7 desc7 author7 email7\n",
                "tool8 'unknown class' v8 date8 desc8 author8 email8\n",
                "tool9 . v9 date9 desc9 author9 email9\n",
                "tool10 ? v10 date10 desc10 author10 email10\n",
            ),
            "c17-software-order.cif",
        );
        let mut metadata = BioMetadata::default();
        let prior = BioSoftwareItem {
            name: "prior".to_owned(),
            version: "prior-version".to_owned(),
            date: "prior-date".to_owned(),
            description: "prior-description".to_owned(),
            contact_author: "prior-author".to_owned(),
            contact_author_email: "prior-email".to_owned(),
            classification: BioSoftwareClassification::DataScaling,
        };
        metadata.software.push(prior.clone());

        parse_mmcif_software_info(&doc.blocks()[0], &mut metadata).unwrap();

        assert_eq!(metadata.software.len(), 12);
        assert_eq!(metadata.software[0], prior);
        let expected_classifications = [
            BioSoftwareClassification::DataCollection,
            BioSoftwareClassification::DataExtraction,
            BioSoftwareClassification::DataProcessing,
            BioSoftwareClassification::DataReduction,
            BioSoftwareClassification::DataScaling,
            BioSoftwareClassification::ModelBuilding,
            BioSoftwareClassification::Phasing,
            BioSoftwareClassification::Refinement,
            BioSoftwareClassification::Unspecified,
            BioSoftwareClassification::Unspecified,
            BioSoftwareClassification::Unspecified,
        ];
        for (index, classification) in expected_classifications.into_iter().enumerate() {
            let item = &metadata.software[index + 1];
            assert_eq!(item.name, format!("tool{index}"));
            assert_eq!(item.version, format!("v{index}"));
            assert_eq!(item.date, format!("date{index}"));
            assert_eq!(item.description, format!("desc{index}"));
            assert_eq!(item.contact_author, format!("author{index}"));
            assert_eq!(item.contact_author_email, format!("email{index}"));
            assert_eq!(item.classification, classification);
        }
    }

    #[test]
    fn bio_mmcif_c17_software_optional_columns_nulls_and_quoted_empty_values() {
        let absent = document(
            "data_c17_software_absent\n_entry.id no-software\n",
            "c17-software-absent.cif",
        );
        let pair_values = document(
            concat!(
                "data_c17_software_pairs\n",
                "_software.name 'Pair tool'\n",
                "_software.version ''\n",
                "_software.date .\n",
                "_software.description ?\n",
                "_software.contact_author 'Ada Lovelace'\n",
            ),
            "c17-software-pairs.cif",
        );
        let mut metadata = BioMetadata::default();
        let prior = BioSoftwareItem {
            name: "existing".to_owned(),
            version: "v-old".to_owned(),
            date: "d-old".to_owned(),
            description: "description-old".to_owned(),
            contact_author: "author-old".to_owned(),
            contact_author_email: "email-old".to_owned(),
            classification: BioSoftwareClassification::Refinement,
        };
        metadata.software.push(prior.clone());

        parse_mmcif_software_info(&absent.blocks()[0], &mut metadata).unwrap();
        assert_eq!(metadata.software, [prior.clone()]);
        parse_mmcif_software_info(&pair_values.blocks()[0], &mut metadata).unwrap();

        assert_eq!(
            metadata.software,
            [
                prior,
                BioSoftwareItem {
                    name: "Pair tool".to_owned(),
                    version: String::new(),
                    date: String::new(),
                    description: String::new(),
                    contact_author: "Ada Lovelace".to_owned(),
                    contact_author_email: String::new(),
                    classification: BioSoftwareClassification::Unspecified,
                },
            ]
        );
    }

    #[test]
    fn bio_mmcif_c18_refinement_rows_preserve_source_order_fields_and_ignore_shell() {
        let doc = document(
            concat!(
                "data_c18_refinement_order\n",
                "loop_\n",
                "_refine.pdbx_refine_id\n",
                "_refine.ls_d_res_high\n",
                "_refine.ls_d_res_low\n",
                "_refine.ls_percent_reflns_obs\n",
                "_refine.ls_number_reflns_obs\n",
                "_refine.ls_number_reflns_R_work\n",
                "_refine.ls_number_reflns_R_free\n",
                "_refine.ls_R_factor_obs\n",
                "_refine.ls_R_factor_R_work\n",
                "_refine.ls_R_factor_R_free\n",
                "independent-A 2.4 3.6 94.5 1200 1100 100 0.21 0.19 0.24\n",
                "unjoined-B 1.8 2.7 98.0 2200 2000 200 0.17 0.16 0.18\n",
                "loop_\n",
                "_refine_ls_shell.pdbx_total_number_of_bins_used\n",
                "_refine_ls_shell.d_res_high\n",
                "3 1.8\n",
                "_em_3d_reconstruction.resolution 0.5\n",
            ),
            "c18-refinement-order.cif",
        );
        let mut source_state = BioStructureSourceState::default();
        let mut metadata = BioMetadata::default();
        let mut prior = BioRefinementInfo::default();
        prior.id = "preexisting".to_owned();
        metadata.refinement.push(prior);

        parse_mmcif_refinement_info(&doc.blocks()[0], &mut source_state, &mut metadata).unwrap();

        assert_eq!(source_state.resolution, 1.8);
        assert_eq!(metadata.refinement.len(), 3);
        assert_eq!(metadata.refinement[0].id, "preexisting");
        assert_eq!(metadata.refinement[1].id, "independent-A");
        assert_eq!(metadata.refinement[2].id, "unjoined-B");

        let first = &metadata.refinement[1];
        assert_eq!(first.basic.resolution_high, 2.4);
        assert_eq!(first.basic.resolution_low, 3.6);
        assert_eq!(first.basic.completeness, 94.5);
        assert_eq!(first.basic.reflection_count, 1200);
        assert_eq!(first.basic.work_set_count, 1100);
        assert_eq!(first.basic.rfree_set_count, 100);
        assert_eq!(first.basic.r_all, 0.21);
        assert_eq!(first.basic.r_work, 0.19);
        assert_eq!(first.basic.r_free, 0.24);

        let second = &metadata.refinement[2];
        assert_eq!(second.basic.resolution_high, 1.8);
        assert_eq!(second.basic.resolution_low, 2.7);
        assert_eq!(second.basic.completeness, 98.0);
        assert_eq!(second.basic.reflection_count, 2200);
        assert_eq!(second.basic.work_set_count, 2000);
        assert_eq!(second.basic.rfree_set_count, 200);
        assert_eq!(second.basic.r_all, 0.17);
        assert_eq!(second.basic.r_work, 0.16);
        assert_eq!(second.basic.r_free, 0.18);

        for refinement in &metadata.refinement {
            assert_eq!(refinement.bin_count, -1);
            assert!(refinement.bins.is_empty());
        }
    }

    #[test]
    fn bio_mmcif_c18_null_missing_and_resolution_fallback_follow_source() {
        let null_values = document(
            concat!(
                "data_c18_null_values\n",
                "loop_\n",
                "_refine.pdbx_refine_id\n",
                "_refine.ls_d_res_high\n",
                "_refine.ls_d_res_low\n",
                "_refine.ls_percent_reflns_obs\n",
                "_refine.ls_number_reflns_obs\n",
                "_refine.ls_number_reflns_R_work\n",
                "_refine.ls_number_reflns_R_free\n",
                "_refine.ls_R_factor_obs\n",
                "_refine.ls_R_factor_R_work\n",
                "_refine.ls_R_factor_R_free\n",
                ". . ? . ? . ? . . .\n",
                "_em_3d_reconstruction.resolution 1.25\n",
            ),
            "c18-null-values.cif",
        );
        let mut source_state = BioStructureSourceState::default();
        let mut metadata = BioMetadata::default();
        parse_mmcif_refinement_info(&null_values.blocks()[0], &mut source_state, &mut metadata)
            .unwrap();

        assert_eq!(metadata.refinement.len(), 1);
        assert_eq!(metadata.refinement[0].id, "");
        assert!(metadata.refinement[0].basic.resolution_high.is_nan());
        assert!(metadata.refinement[0].basic.resolution_low.is_nan());
        assert!(metadata.refinement[0].basic.completeness.is_nan());
        assert_eq!(metadata.refinement[0].basic.reflection_count, -1);
        assert_eq!(metadata.refinement[0].basic.work_set_count, -1);
        assert_eq!(metadata.refinement[0].basic.rfree_set_count, -1);
        assert!(metadata.refinement[0].basic.r_all.is_nan());
        assert!(metadata.refinement[0].basic.r_work.is_nan());
        assert!(metadata.refinement[0].basic.r_free.is_nan());
        assert_eq!(source_state.resolution, 1.25);

        let missing_columns = document(
            concat!(
                "data_c18_missing_columns\n",
                "_refine.pdbx_refine_id 'pair ID'\n",
                "_refine.ls_d_res_high 3.0\n",
                "_em_3d_reconstruction.resolution 0.5\n",
            ),
            "c18-missing-columns.cif",
        );
        let mut source_state = BioStructureSourceState::default();
        source_state.resolution = 2.25;
        let mut metadata = BioMetadata::default();
        parse_mmcif_refinement_info(
            &missing_columns.blocks()[0],
            &mut source_state,
            &mut metadata,
        )
        .unwrap();

        assert_eq!(metadata.refinement.len(), 1);
        assert_eq!(metadata.refinement[0].id, "pair ID");
        assert_eq!(metadata.refinement[0].basic.resolution_high, 3.0);
        assert!(metadata.refinement[0].basic.resolution_low.is_nan());
        assert!(metadata.refinement[0].basic.completeness.is_nan());
        assert_eq!(metadata.refinement[0].basic.reflection_count, -1);
        assert_eq!(metadata.refinement[0].basic.work_set_count, -1);
        assert_eq!(metadata.refinement[0].basic.rfree_set_count, -1);
        assert_eq!(source_state.resolution, 2.25);
    }

    #[test]
    fn bio_mmcif_c18_bad_integer_propagates_after_source_order_updates() {
        let doc = document(
            concat!(
                "data_c18_bad_integer\n",
                "loop_\n",
                "_refine.pdbx_refine_id\n",
                "_refine.ls_d_res_high\n",
                "_refine.ls_d_res_low\n",
                "_refine.ls_percent_reflns_obs\n",
                "_refine.ls_number_reflns_obs\n",
                "_refine.ls_number_reflns_R_work\n",
                "_refine.ls_number_reflns_R_free\n",
                "_refine.ls_R_factor_obs\n",
                "_refine.ls_R_factor_R_work\n",
                "_refine.ls_R_factor_R_free\n",
                "bad-integer 2.0 3.0 80.0 12x 10 2 0.2 0.1 0.3\n",
            ),
            "c18-bad-integer.cif",
        );
        let mut source_state = BioStructureSourceState::default();
        let mut metadata = BioMetadata::default();

        let error = parse_mmcif_refinement_info(&doc.blocks()[0], &mut source_state, &mut metadata)
            .unwrap_err();

        assert_eq!(error.kind(), CifReadErrorKind::InvalidValue);
        assert_eq!(error.source(), "cif-value");
        assert_eq!(metadata.refinement.len(), 1);
        assert_eq!(metadata.refinement[0].id, "bad-integer");
        assert_eq!(metadata.refinement[0].basic.resolution_high, 2.0);
        assert_eq!(metadata.refinement[0].basic.resolution_low, 3.0);
        assert_eq!(metadata.refinement[0].basic.completeness, 80.0);
        assert_eq!(metadata.refinement[0].basic.reflection_count, -1);
        assert_eq!(source_state.resolution, 2.0);
    }

    #[test]
    fn bio_mmcif_c19_tls_groups_route_and_preserve_source_matrix_order() {
        let tensor_values = (1..=24)
            .map(|value| value.to_string())
            .collect::<Vec<_>>()
            .join(" ");
        let zero_values = ["0"; 24].join(" ");
        let mut null_values = ["0"; 24];
        null_values[0] = ".";
        let null_values = null_values.join(" ");
        let doc = document(
            &format!(
                concat!(
                    "data_c19_tls_groups\n",
                    "loop_\n",
                    "_pdbx_refine_tls.id\n",
                    "_pdbx_refine_tls.pdbx_refine_id\n",
                    "_pdbx_refine_tls.T[1][1]\n",
                    "_pdbx_refine_tls.T[2][2]\n",
                    "_pdbx_refine_tls.T[3][3]\n",
                    "_pdbx_refine_tls.T[1][2]\n",
                    "_pdbx_refine_tls.T[1][3]\n",
                    "_pdbx_refine_tls.T[2][3]\n",
                    "_pdbx_refine_tls.L[1][1]\n",
                    "_pdbx_refine_tls.L[2][2]\n",
                    "_pdbx_refine_tls.L[3][3]\n",
                    "_pdbx_refine_tls.L[1][2]\n",
                    "_pdbx_refine_tls.L[1][3]\n",
                    "_pdbx_refine_tls.L[2][3]\n",
                    "_pdbx_refine_tls.S[1][1]\n",
                    "_pdbx_refine_tls.S[1][2]\n",
                    "_pdbx_refine_tls.S[1][3]\n",
                    "_pdbx_refine_tls.S[2][1]\n",
                    "_pdbx_refine_tls.S[2][2]\n",
                    "_pdbx_refine_tls.S[2][3]\n",
                    "_pdbx_refine_tls.S[3][1]\n",
                    "_pdbx_refine_tls.S[3][2]\n",
                    "_pdbx_refine_tls.S[3][3]\n",
                    "_pdbx_refine_tls.origin_x\n",
                    "_pdbx_refine_tls.origin_y\n",
                    "_pdbx_refine_tls.origin_z\n",
                    "12tail R2 {tensor_values}\n",
                    "40000 unmatched {zero_values}\n",
                    "+7 . {null_values}\n",
                ),
                tensor_values = tensor_values,
                zero_values = zero_values,
                null_values = null_values,
            ),
            "c19-tls-groups.cif",
        );
        let mut first = BioRefinementInfo::default();
        first.id = "R1".to_owned();
        let mut second = BioRefinementInfo::default();
        second.id = "R2".to_owned();
        let mut metadata = BioMetadata {
            refinement: vec![first, second],
            ..BioMetadata::default()
        };

        parse_mmcif_tls_info(&doc.blocks()[0], &mut metadata).unwrap();

        assert_eq!(metadata.refinement[0].tls_groups.len(), 2);
        assert_eq!(metadata.refinement[1].tls_groups.len(), 1);
        let matched = &metadata.refinement[1].tls_groups[0];
        assert_eq!(matched.id, "12tail");
        assert_eq!(matched.num_id, 12);
        assert_eq!(matched.t, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        assert_eq!(matched.l, [7.0, 8.0, 9.0, 10.0, 11.0, 12.0]);
        assert_eq!(
            matched.s,
            [[13.0, 14.0, 15.0], [16.0, 17.0, 18.0], [19.0, 20.0, 21.0],]
        );
        assert_eq!(matched.origin, [22.0, 23.0, 24.0]);

        let fallback = &metadata.refinement[0].tls_groups;
        assert_eq!(fallback[0].id, "40000");
        assert_eq!(fallback[0].num_id, -25_536);
        assert_eq!(fallback[0].t, [0.0; 6]);
        assert_eq!(fallback[0].l, [0.0; 6]);
        assert_eq!(fallback[0].s, [[0.0; 3]; 3]);
        assert_eq!(fallback[0].origin, [0.0; 3]);
        assert_eq!(fallback[1].id, "+7");
        assert_eq!(fallback[1].num_id, 0);
        assert!(fallback[1].t[0].is_nan());
        assert_eq!(fallback[1].t[1..], [0.0; 5]);
    }

    #[test]
    fn bio_mmcif_c19_tls_selections_preserve_first_match_order_and_presence() {
        let doc = document(
            concat!(
                "data_c19_tls_selections\n",
                "loop_\n",
                "_pdbx_refine_tls_group.refine_tls_id\n",
                "_pdbx_refine_tls_group.beg_auth_asym_id\n",
                "_pdbx_refine_tls_group.beg_auth_seq_id\n",
                "_pdbx_refine_tls_group.beg_PDB_ins_code\n",
                "_pdbx_refine_tls_group.end_auth_seq_id\n",
                "_pdbx_refine_tls_group.end_PDB_ins_code\n",
                "_pdbx_refine_tls_group.selection_details\n",
                "shared A 15A . 24 . 'first selection'\n",
                "unmatched Z 1 . 2 . ignored\n",
                "shared B 17 . 31 . second\n",
                "shared . . ? . ? .\n",
            ),
            "c19-tls-selections.cif",
        );
        let mut first = BioRefinementInfo::default();
        first.id = "R1".to_owned();
        let mut first_group = BioTlsGroup::default();
        first_group.id = "shared".to_owned();
        first.tls_groups.push(first_group);
        let mut second = BioRefinementInfo::default();
        second.id = "R2".to_owned();
        let mut duplicate_group = BioTlsGroup::default();
        duplicate_group.id = "shared".to_owned();
        second.tls_groups.push(duplicate_group);
        let mut metadata = BioMetadata {
            refinement: vec![first, second],
            ..BioMetadata::default()
        };

        parse_mmcif_tls_info(&doc.blocks()[0], &mut metadata).unwrap();

        let absent_optional_columns = document(
            "data_c19_tls_missing_columns\n_pdbx_refine_tls_group.refine_tls_id shared\n",
            "c19-tls-missing-columns.cif",
        );
        parse_mmcif_tls_info(&absent_optional_columns.blocks()[0], &mut metadata).unwrap();
        assert_eq!(metadata.refinement[0].tls_groups[0].selections.len(), 4);
        let first_matches = &metadata.refinement[0].tls_groups[0].selections;
        assert_eq!(
            first_matches[0],
            BioTlsSelection {
                chain: PdbChainId::from_ascii(b"A").unwrap(),
                res_begin: PdbSeqId::new(15, Some(b'A')),
                res_end: PdbSeqId::new(24, None),
                details: "first selection".to_owned(),
            }
        );
        assert_eq!(
            first_matches[1],
            BioTlsSelection {
                chain: PdbChainId::from_ascii(b"B").unwrap(),
                res_begin: PdbSeqId::new(17, None),
                res_end: PdbSeqId::new(31, None),
                details: "second".to_owned(),
            }
        );
        assert_eq!(
            first_matches[2],
            BioTlsSelection {
                chain: PdbChainId::from_ascii(b"").unwrap(),
                res_begin: PdbSeqId::new(i32::MIN, None),
                res_end: PdbSeqId::new(i32::MIN, None),
                details: String::new(),
            }
        );
        assert!(metadata.refinement[1].tls_groups[0].selections.is_empty());
        assert_eq!(
            metadata.refinement[0].tls_groups[0].selections[3],
            BioTlsSelection {
                chain: PdbChainId::from_ascii(b"").unwrap(),
                res_begin: PdbSeqId::new(i32::MIN, None),
                res_end: PdbSeqId::new(i32::MIN, None),
                details: String::new(),
            }
        );
    }

    #[test]
    fn bio_mmcif_c19_tls_selection_failure_retains_source_ordered_append() {
        let numeric_values = ["0"; 24].join(" ");
        let text = String::from(concat!(
            "data_c19_tls_selection_error\n",
            "loop_\n",
            "_pdbx_refine_tls.id\n",
            "_pdbx_refine_tls.pdbx_refine_id\n",
            "_pdbx_refine_tls.T[1][1]\n",
            "_pdbx_refine_tls.T[2][2]\n",
            "_pdbx_refine_tls.T[3][3]\n",
            "_pdbx_refine_tls.T[1][2]\n",
            "_pdbx_refine_tls.T[1][3]\n",
            "_pdbx_refine_tls.T[2][3]\n",
            "_pdbx_refine_tls.L[1][1]\n",
            "_pdbx_refine_tls.L[2][2]\n",
            "_pdbx_refine_tls.L[3][3]\n",
            "_pdbx_refine_tls.L[1][2]\n",
            "_pdbx_refine_tls.L[1][3]\n",
            "_pdbx_refine_tls.L[2][3]\n",
            "_pdbx_refine_tls.S[1][1]\n",
            "_pdbx_refine_tls.S[1][2]\n",
            "_pdbx_refine_tls.S[1][3]\n",
            "_pdbx_refine_tls.S[2][1]\n",
            "_pdbx_refine_tls.S[2][2]\n",
            "_pdbx_refine_tls.S[2][3]\n",
            "_pdbx_refine_tls.S[3][1]\n",
            "_pdbx_refine_tls.S[3][2]\n",
            "_pdbx_refine_tls.S[3][3]\n",
            "_pdbx_refine_tls.origin_x\n",
            "_pdbx_refine_tls.origin_y\n",
            "_pdbx_refine_tls.origin_z\n",
            "group R1 {values}\n",
            "loop_\n",
            "_pdbx_refine_tls_group.refine_tls_id\n",
            "_pdbx_refine_tls_group.beg_auth_asym_id\n",
            "_pdbx_refine_tls_group.beg_auth_seq_id\n",
            "_pdbx_refine_tls_group.beg_PDB_ins_code\n",
            "_pdbx_refine_tls_group.end_auth_seq_id\n",
            "_pdbx_refine_tls_group.end_PDB_ins_code\n",
            "_pdbx_refine_tls_group.selection_details\n",
            "group A 15A B 20 . detail\n",
        ))
        .replace("{values}", &numeric_values);
        let doc = document(&text, "c19-tls-selection-error.cif");
        let mut refinement = BioRefinementInfo::default();
        refinement.id = "R1".to_owned();
        let mut metadata = BioMetadata {
            refinement: vec![refinement],
            ..BioMetadata::default()
        };

        let error = parse_mmcif_tls_info(&doc.blocks()[0], &mut metadata).unwrap_err();

        assert_eq!(
            error,
            MmcifTlsInfoError::SequenceId(MmcifSequenceIdError::InconsistentInsertionCode {
                sequence_id: "15A".to_owned(),
            })
        );
        let selection = &metadata.refinement[0].tls_groups[0].selections[0];
        assert_eq!(metadata.refinement[0].tls_groups[0].selections.len(), 1);
        assert_eq!(selection.chain, PdbChainId::from_ascii(b"A").unwrap());
        assert_eq!(selection.res_begin, PdbSeqId::new(i32::MIN, None));
        assert_eq!(selection.details, "");
    }

    #[test]
    fn bio_mmcif_c19_tls_empty_refinements_skip_group_value_conversion() {
        let values = ["0"; 24].join(" ");
        let text = format!(
            concat!(
                "data_c19_tls_no_refinement\n",
                "loop_\n",
                "_pdbx_refine_tls.id\n",
                "_pdbx_refine_tls.pdbx_refine_id\n",
                "_pdbx_refine_tls.T[1][1]\n",
                "_pdbx_refine_tls.T[2][2]\n",
                "_pdbx_refine_tls.T[3][3]\n",
                "_pdbx_refine_tls.T[1][2]\n",
                "_pdbx_refine_tls.T[1][3]\n",
                "_pdbx_refine_tls.T[2][3]\n",
                "_pdbx_refine_tls.L[1][1]\n",
                "_pdbx_refine_tls.L[2][2]\n",
                "_pdbx_refine_tls.L[3][3]\n",
                "_pdbx_refine_tls.L[1][2]\n",
                "_pdbx_refine_tls.L[1][3]\n",
                "_pdbx_refine_tls.L[2][3]\n",
                "_pdbx_refine_tls.S[1][1]\n",
                "_pdbx_refine_tls.S[1][2]\n",
                "_pdbx_refine_tls.S[1][3]\n",
                "_pdbx_refine_tls.S[2][1]\n",
                "_pdbx_refine_tls.S[2][2]\n",
                "_pdbx_refine_tls.S[2][3]\n",
                "_pdbx_refine_tls.S[3][1]\n",
                "_pdbx_refine_tls.S[3][2]\n",
                "_pdbx_refine_tls.S[3][3]\n",
                "_pdbx_refine_tls.origin_x\n",
                "_pdbx_refine_tls.origin_y\n",
                "_pdbx_refine_tls.origin_z\n",
                "2147483648 . {values}\n",
                "loop_\n",
                "_pdbx_refine_tls_group.refine_tls_id\n",
                "2147483648\n",
            ),
            values = values,
        );
        let doc = document(&text, "c19-tls-no-refinement.cif");
        let mut metadata = BioMetadata::default();

        parse_mmcif_tls_info(&doc.blocks()[0], &mut metadata).unwrap();

        assert!(metadata.refinement.is_empty());
    }

    #[test]
    fn bio_mmcif_c20_experiments_and_diffractions_follow_source_order_and_first_matches() {
        let doc = document(
            concat!(
                "data_c20_experimental_order\n",
                "loop_\n",
                "_exptl.method\n",
                "_exptl.crystals_number\n",
                "'x-ray diffraction' 2\n",
                "'neutron diffraction' .\n",
                "loop_\n",
                "_exptl_crystal.id\n",
                "_exptl_crystal.description\n",
                "X 'first X'\n",
                "Y 'first Y'\n",
                "X 'duplicate X'\n",
                "loop_\n",
                "_diffrn.id\n",
                "_diffrn.crystal_id\n",
                "_diffrn.ambient_temp\n",
                "shared Y .\n",
                "x-only X 100.5\n",
                "shared X 300.25\n",
                "orphan missing 288\n",
                "loop_\n",
                "_diffrn_detector.diffrn_id\n",
                "_diffrn_detector.pdbx_collection_date\n",
                "_diffrn_detector.detector\n",
                "_diffrn_detector.type\n",
                "_diffrn_detector.details\n",
                "shared 2026-01-02 'pixel detector' CCD 'mirror optics'\n",
                "orphan ignored ignored ignored ignored\n",
                "loop_\n",
                "_diffrn_radiation.diffrn_id\n",
                "_diffrn_radiation.pdbx_scattering_type\n",
                "_diffrn_radiation.pdbx_monochromatic_or_laue_m_l\n",
                "_diffrn_radiation.monochromator\n",
                "shared 'x-ray' M 'Si(111)'\n",
                "loop_\n",
                "_diffrn_source.diffrn_id\n",
                "_diffrn_source.source\n",
                "_diffrn_source.type\n",
                "_diffrn_source.pdbx_synchrotron_site\n",
                "_diffrn_source.pdbx_synchrotron_beamline\n",
                "_diffrn_source.pdbx_wavelength_list\n",
                "shared 'storage ring' synchrotron 'Diamond' I24 '0.98, 1.00'\n",
            ),
            "c20-experimental-order.cif",
        );
        let mut metadata = BioMetadata::default();

        parse_mmcif_experimental_info(&doc.blocks()[0], &mut metadata).unwrap();

        assert_eq!(metadata.experiments.len(), 2);
        assert_eq!(metadata.experiments[0].method, "x-ray diffraction");
        assert_eq!(metadata.experiments[0].number_of_crystals, 2);
        assert_eq!(metadata.experiments[1].method, "neutron diffraction");
        assert_eq!(metadata.experiments[1].number_of_crystals, -1);
        assert_eq!(metadata.crystals.len(), 3);
        assert_eq!(metadata.crystals[0].description, "first X");
        assert_eq!(metadata.crystals[1].description, "first Y");
        assert_eq!(metadata.crystals[2].description, "duplicate X");

        let first_x = &metadata.crystals[0].diffractions;
        assert_eq!(first_x.len(), 2);
        assert_eq!(first_x[0].id, "x-only");
        assert_eq!(first_x[0].temperature, 100.5);
        assert_eq!(first_x[1].id, "shared");
        assert_eq!(first_x[1].temperature, 300.25);
        assert_eq!(first_x[1].collection_date, "2026-01-02");
        assert_eq!(first_x[1].detector, "pixel detector");
        assert_eq!(first_x[1].detector_make, "CCD");
        assert_eq!(first_x[1].optics, "mirror optics");
        assert_eq!(first_x[1].scattering_type, "x-ray");
        assert_eq!(first_x[1].mono_or_laue, b'M');
        assert_eq!(first_x[1].monochromator, "Si(111)");
        assert_eq!(first_x[1].source, "storage ring");
        assert_eq!(first_x[1].source_type, "synchrotron");
        assert_eq!(first_x[1].synchrotron, "Diamond");
        assert_eq!(first_x[1].beamline, "I24");
        assert_eq!(first_x[1].wavelengths, "0.98, 1.00");

        assert_eq!(metadata.crystals[1].diffractions.len(), 1);
        assert_eq!(metadata.crystals[1].diffractions[0].id, "shared");
        assert!(metadata.crystals[1].diffractions[0].temperature.is_nan());
        assert_eq!(metadata.crystals[1].diffractions[0].detector, "");
        assert!(metadata.crystals[2].diffractions.is_empty());
    }

    #[test]
    fn bio_mmcif_c20_missing_null_and_quoted_empty_values_follow_has2() {
        let doc = document(
            concat!(
                "data_c20_optional_values\n",
                "_exptl.method 'single crystal'\n",
                "_diffrn_detector.diffrn_id keep\n",
                "loop_\n",
                "_diffrn_radiation.diffrn_id\n",
                "_diffrn_radiation.pdbx_scattering_type\n",
                "_diffrn_radiation.pdbx_monochromatic_or_laue_m_l\n",
                "_diffrn_radiation.monochromator\n",
                "keep . '' ?\n",
                "loop_\n",
                "_diffrn_source.diffrn_id\n",
                "_diffrn_source.source\n",
                "_diffrn_source.type\n",
                "_diffrn_source.pdbx_synchrotron_site\n",
                "_diffrn_source.pdbx_synchrotron_beamline\n",
                "_diffrn_source.pdbx_wavelength_list\n",
                "keep . '' ? . ''\n",
            ),
            "c20-optional-values.cif",
        );
        let mut diffraction = BioDiffractionInfo {
            id: "keep".to_owned(),
            collection_date: "old date".to_owned(),
            detector: "old detector".to_owned(),
            detector_make: "old make".to_owned(),
            optics: "old optics".to_owned(),
            scattering_type: "old scattering".to_owned(),
            mono_or_laue: b'M',
            monochromator: "old monochromator".to_owned(),
            source: "old source".to_owned(),
            source_type: "old type".to_owned(),
            synchrotron: "old site".to_owned(),
            beamline: "old beamline".to_owned(),
            wavelengths: "old wavelengths".to_owned(),
            ..BioDiffractionInfo::default()
        };
        diffraction.temperature = 123.0;
        let mut crystal = BioExperimentalCrystalInfo {
            id: "existing".to_owned(),
            ..BioExperimentalCrystalInfo::default()
        };
        crystal.diffractions.push(diffraction);
        let mut metadata = BioMetadata {
            crystals: vec![crystal],
            ..BioMetadata::default()
        };

        parse_mmcif_experimental_info(&doc.blocks()[0], &mut metadata).unwrap();

        assert_eq!(metadata.experiments.len(), 1);
        assert_eq!(metadata.experiments[0].method, "single crystal");
        assert_eq!(metadata.experiments[0].number_of_crystals, -1);
        let updated = &metadata.crystals[0].diffractions[0];
        assert_eq!(updated.temperature, 123.0);
        assert_eq!(updated.collection_date, "old date");
        assert_eq!(updated.detector, "old detector");
        assert_eq!(updated.detector_make, "old make");
        assert_eq!(updated.optics, "old optics");
        assert_eq!(updated.scattering_type, "old scattering");
        assert_eq!(updated.mono_or_laue, 0);
        assert_eq!(updated.monochromator, "old monochromator");
        assert_eq!(updated.source, "old source");
        assert_eq!(updated.source_type, "");
        assert_eq!(updated.synchrotron, "old site");
        assert_eq!(updated.beamline, "old beamline");
        assert_eq!(updated.wavelengths, "");
    }

    #[test]
    fn bio_mmcif_c20_bad_crystal_count_retains_source_ordered_experiment_append() {
        let doc = document(
            concat!(
                "data_c20_bad_crystal_count\n",
                "_exptl.method 'invalid count is parsed after method'\n",
                "_exptl.crystals_number not-an-integer\n",
            ),
            "c20-bad-crystal-count.cif",
        );
        let mut metadata = BioMetadata::default();

        let error = parse_mmcif_experimental_info(&doc.blocks()[0], &mut metadata).unwrap_err();

        assert_eq!(error.kind(), CifReadErrorKind::InvalidValue);
        assert_eq!(error.message(), "not an integer: not-an-integer");
        assert_eq!(metadata.experiments.len(), 1);
        assert_eq!(
            metadata.experiments[0].method,
            "invalid count is parsed after method"
        );
        assert_eq!(metadata.experiments[0].number_of_crystals, -1);
        assert!(metadata.crystals.is_empty());
    }

    #[test]
    fn bio_mmcif_c21_reflns_rows_pair_by_order_and_append_split_fields() {
        let doc = document(
            concat!(
                "data_c21_reflns_order\n",
                "loop_\n",
                "_reflns.pdbx_diffrn_id\n",
                "_reflns.number_obs\n",
                "_reflns.d_resolution_high\n",
                "_reflns.d_resolution_low\n",
                "_reflns.percent_possible_obs\n",
                "_reflns.pdbx_redundancy\n",
                "_reflns.pdbx_Rmerge_I_obs\n",
                "_reflns.pdbx_Rsym_value\n",
                "_reflns.pdbx_netI_over_sigmaI\n",
                "',D1,,D2,' 12 2.25 3.75 97.5 1.4 0.02 0.03 8.5\n",
                "D3 34 4.5 5.5 89.0 2.5 0.04 0.05 6.5\n",
            ),
            "c21-reflns-order.cif",
        );
        let mut first = BioExperimentInfo::default();
        first.diffraction_ids.push("preexisting".to_owned());
        let mut metadata = BioMetadata::default();
        metadata.experiments = vec![first, BioExperimentInfo::default()];

        parse_mmcif_reflns_info(&doc.blocks()[0], &mut metadata).unwrap();

        let first = &metadata.experiments[0];
        assert_eq!(
            first.diffraction_ids,
            ["preexisting", "", "D1", "", "D2", ""]
        );
        assert_eq!(first.unique_reflections, 12);
        assert_eq!(first.reflections.resolution_high, 2.25);
        assert_eq!(first.reflections.resolution_low, 3.75);
        assert_eq!(first.reflections.completeness, 97.5);
        assert_eq!(first.reflections.redundancy, 1.4);
        assert_eq!(first.reflections.r_merge, 0.02);
        assert_eq!(first.reflections.r_sym, 0.03);
        assert_eq!(first.reflections.mean_i_over_sigma, 8.5);

        let second = &metadata.experiments[1];
        assert_eq!(second.diffraction_ids, ["D3"]);
        assert_eq!(second.unique_reflections, 34);
        assert_eq!(second.reflections.resolution_high, 4.5);
        assert_eq!(second.reflections.resolution_low, 5.5);
        assert_eq!(second.reflections.completeness, 89.0);
        assert_eq!(second.reflections.redundancy, 2.5);
        assert_eq!(second.reflections.r_merge, 0.04);
        assert_eq!(second.reflections.r_sym, 0.05);
        assert_eq!(second.reflections.mean_i_over_sigma, 6.5);
    }

    #[test]
    fn bio_mmcif_c21_missing_null_and_shell_values_preserve_source_state() {
        let missing = document(
            concat!(
                "data_c21_reflns_missing\n",
                "loop_\n",
                "_reflns.pdbx_diffrn_id\n",
                "default-experiment\n",
            ),
            "c21-reflns-missing.cif",
        );
        let mut metadata = BioMetadata::default();
        metadata.experiments.push(BioExperimentInfo::default());
        parse_mmcif_reflns_info(&missing.blocks()[0], &mut metadata).unwrap();

        let defaults = &metadata.experiments[0];
        assert_eq!(defaults.diffraction_ids, ["default-experiment"]);
        assert_eq!(defaults.unique_reflections, -1);
        assert!(defaults.reflections.resolution_high.is_nan());
        assert!(defaults.reflections.resolution_low.is_nan());
        assert!(defaults.reflections.completeness.is_nan());
        assert!(defaults.reflections.redundancy.is_nan());
        assert!(defaults.reflections.r_merge.is_nan());
        assert!(defaults.reflections.r_sym.is_nan());
        assert!(defaults.reflections.mean_i_over_sigma.is_nan());
        assert!(defaults.b_wilson.is_nan());
        assert!(defaults.shells.is_empty());

        let nulls_and_unread_categories = document(
            concat!(
                "data_c21_reflns_nulls\n",
                "loop_\n",
                "_reflns.pdbx_diffrn_id\n",
                "_reflns.number_obs\n",
                "_reflns.d_resolution_high\n",
                "_reflns.d_resolution_low\n",
                "_reflns.percent_possible_obs\n",
                "_reflns.pdbx_redundancy\n",
                "_reflns.pdbx_Rmerge_I_obs\n",
                "_reflns.pdbx_Rsym_value\n",
                "_reflns.pdbx_netI_over_sigmaI\n",
                ". ? . ? . ? . ? .\n",
                "_reflns.B_iso_Wilson_estimate 41.0\n",
                "loop_\n",
                "_reflns_shell.d_res_high\n",
                "0.75\n",
            ),
            "c21-reflns-null-shell.cif",
        );
        let mut experiment = BioExperimentInfo {
            unique_reflections: 73,
            reflections: BioReflectionsInfo {
                resolution_high: 1.25,
                resolution_low: 2.5,
                completeness: 88.0,
                redundancy: 3.0,
                r_merge: 0.11,
                r_sym: 0.12,
                mean_i_over_sigma: 4.0,
            },
            b_wilson: 19.5,
            shells: vec![BioReflectionsInfo {
                resolution_high: 0.8,
                resolution_low: 1.1,
                completeness: 82.0,
                redundancy: 2.0,
                r_merge: 0.08,
                r_sym: 0.09,
                mean_i_over_sigma: 3.0,
            }],
            ..BioExperimentInfo::default()
        };
        experiment.diffraction_ids.push("retained".to_owned());
        let original = experiment.clone();
        let mut metadata = BioMetadata::default();
        metadata.experiments.push(experiment);

        parse_mmcif_reflns_info(&nulls_and_unread_categories.blocks()[0], &mut metadata).unwrap();

        let updated = &metadata.experiments[0];
        assert_eq!(updated.diffraction_ids, ["retained", ""]);
        assert_eq!(updated.unique_reflections, original.unique_reflections);
        assert_eq!(updated.reflections, original.reflections);
        assert_eq!(updated.b_wilson, original.b_wilson);
        assert_eq!(updated.shells, original.shells);
    }

    #[test]
    fn bio_mmcif_c21_excess_rows_break_before_invalid_conversion() {
        let doc = document(
            concat!(
                "data_c21_reflns_excess\n",
                "loop_\n",
                "_reflns.pdbx_diffrn_id\n",
                "_reflns.number_obs\n",
                "_reflns.d_resolution_high\n",
                "_reflns.d_resolution_low\n",
                "_reflns.percent_possible_obs\n",
                "_reflns.pdbx_redundancy\n",
                "_reflns.pdbx_Rmerge_I_obs\n",
                "_reflns.pdbx_Rsym_value\n",
                "_reflns.pdbx_netI_over_sigmaI\n",
                "consumed 8 2.0 3.0 90.0 1.5 0.1 0.2 5.0\n",
                "unconsumed not-an-int not-a-float . . . . . .\n",
            ),
            "c21-reflns-excess.cif",
        );
        let mut metadata = BioMetadata::default();
        metadata.experiments.push(BioExperimentInfo::default());

        parse_mmcif_reflns_info(&doc.blocks()[0], &mut metadata).unwrap();

        assert_eq!(metadata.experiments[0].diffraction_ids, ["consumed"]);
        assert_eq!(metadata.experiments[0].unique_reflections, 8);
        assert_eq!(metadata.experiments[0].reflections.resolution_high, 2.0);
        assert_eq!(metadata.experiments[0].reflections.mean_i_over_sigma, 5.0);
    }

    #[test]
    fn bio_mmcif_c21_bad_integer_fails_after_id_split_append() {
        let doc = document(
            concat!(
                "data_c21_reflns_bad_integer\n",
                "loop_\n",
                "_reflns.pdbx_diffrn_id\n",
                "_reflns.number_obs\n",
                "_reflns.d_resolution_high\n",
                "_reflns.d_resolution_low\n",
                "_reflns.percent_possible_obs\n",
                "_reflns.pdbx_redundancy\n",
                "_reflns.pdbx_Rmerge_I_obs\n",
                "_reflns.pdbx_Rsym_value\n",
                "_reflns.pdbx_netI_over_sigmaI\n",
                "'left,,right' invalid 2.0 3.0 90.0 1.5 0.1 0.2 5.0\n",
            ),
            "c21-reflns-bad-integer.cif",
        );
        let mut metadata = BioMetadata::default();
        metadata.experiments.push(BioExperimentInfo::default());

        let error = parse_mmcif_reflns_info(&doc.blocks()[0], &mut metadata).unwrap_err();

        assert_eq!(error.kind(), CifReadErrorKind::InvalidValue);
        assert_eq!(error.source(), "cif-value");
        assert_eq!(
            metadata.experiments[0].diffraction_ids,
            ["left", "", "right"]
        );
        assert_eq!(metadata.experiments[0].unique_reflections, -1);
        assert!(metadata.experiments[0].reflections.resolution_high.is_nan());
    }

    fn atom_site_table(block: &CifBlock) -> CifTable<'_> {
        block
            .find("_atom_site.", &ATOM_SITE_TAGS)
            .expect("fixed source-shaped atom_site table")
    }

    fn c05_scalar_row() -> [&'static str; 24] {
        [
            "1", "ATOM", "C", "CA", ".", "GLY", "A", "E", "1", ".", "0", "0", "0", ".", ".", ".",
            "1", "GLY", "A", "CA", "1", ".", ".", ".",
        ]
    }

    fn c05_scalar_document(
        rows: &[[&str; 24]],
        omitted_columns: &[usize],
    ) -> crate::cif::CifDocument {
        let mut input = String::from("data_c05_scalars\nloop_\n");
        for (column, tag) in ATOM_SITE_TAGS.iter().enumerate() {
            if omitted_columns.contains(&column) {
                continue;
            }
            input.push_str("_atom_site.");
            input.push_str(tag.strip_prefix('?').unwrap_or(tag));
            input.push('\n');
        }
        for row in rows {
            for (column, value) in row.iter().enumerate() {
                if !omitted_columns.contains(&column) {
                    input.push_str(value);
                    input.push(' ');
                }
            }
            input.push('\n');
        }
        document(&input, "c05-scalars.cif")
    }

    fn c25_connection_row() -> [&'static str; 22] {
        [
            "c1", "covale", "A1", "B2", "x", "y", "LIG", "GLY", "C1", "N", "A", ".", "15A", "16",
            "1", "2", ".", "?", "1_555", "2_565", "1.25(8)", "LINK",
        ]
    }

    fn c25_connection_document(
        rows: &[[&str; 22]],
        omitted_columns: &[usize],
    ) -> crate::cif::CifDocument {
        let mut input = String::from("data_c25_connections\nloop_\n");
        for (column, tag) in STRUCT_CONN_TAGS.iter().enumerate() {
            if omitted_columns.contains(&column) {
                continue;
            }
            input.push_str("_struct_conn.");
            input.push_str(tag.strip_prefix('?').unwrap_or(tag));
            input.push('\n');
        }
        for row in rows {
            for (column, value) in row.iter().enumerate() {
                if !omitted_columns.contains(&column) {
                    input.push_str(value);
                    input.push(' ');
                }
            }
            input.push('\n');
        }
        document(&input, "c25-connections.cif")
    }

    fn c26_cis_peptide_row() -> [&'static str; 13] {
        [
            "2", "A", "15A", ".", "SER", "CYS", "B", "16", "C", "GLY", "ALA", "X", "-2.5",
        ]
    }

    fn c26_cis_peptide_document(
        rows: &[[&str; 13]],
        omitted_columns: &[usize],
    ) -> crate::cif::CifDocument {
        let mut input = String::from("data_c26_cis_peptides\nloop_\n");
        for (column, tag) in CIS_PEPTIDE_TAGS.iter().enumerate() {
            if omitted_columns.contains(&column) {
                continue;
            }
            input.push_str("_struct_mon_prot_cis.");
            input.push_str(tag.strip_prefix('?').unwrap_or(tag));
            input.push('\n');
        }
        for row in rows {
            for (column, value) in row.iter().enumerate() {
                if !omitted_columns.contains(&column) {
                    input.push_str(value);
                    input.push(' ');
                }
            }
            input.push('\n');
        }
        document(&input, "c26-cis-peptides.cif")
    }

    fn c27_modified_residue_row() -> [&'static str; 8] {
        [
            "A",
            "15A",
            ".",
            "MSE",
            "SEC",
            "MET",
            "'oxidized sulfur'",
            "REFMAC1",
        ]
    }

    fn c27_modified_residue_document(
        rows: &[[&str; 8]],
        omitted_columns: &[usize],
    ) -> crate::cif::CifDocument {
        let mut input = String::from("data_c27_modified_residues\nloop_\n");
        for (column, tag) in MODIFIED_RESIDUE_TAGS.iter().enumerate() {
            if omitted_columns.contains(&column) {
                continue;
            }
            input.push_str("_pdbx_struct_mod_residue.");
            input.push_str(tag.strip_prefix('?').unwrap_or(tag));
            input.push('\n');
        }
        for row in rows {
            for (column, value) in row.iter().enumerate() {
                if !omitted_columns.contains(&column) {
                    input.push_str(value);
                    input.push(' ');
                }
            }
            input.push('\n');
        }
        document(&input, "c27-modified-residues.cif")
    }

    fn c29_assembly_document(
        assembly_rows: &[[&str; 5]],
        property_rows: &[[&str; 3]],
        generator_rows: &[[&str; 3]],
        operator_rows: &[(&str, &str, [&str; 3])],
    ) -> crate::cif::CifDocument {
        let mut input = String::from("data_c29_assemblies\n");
        if !property_rows.is_empty() {
            input.push_str("loop_\n");
            for tag in ["biol_id", "type", "value"] {
                input.push_str("_pdbx_struct_assembly_prop.");
                input.push_str(tag);
                input.push('\n');
            }
            for row in property_rows {
                input.push_str(&row.join(" "));
                input.push('\n');
            }
        }
        if !generator_rows.is_empty() {
            input.push_str("loop_\n");
            for tag in ["assembly_id", "oper_expression", "asym_id_list"] {
                input.push_str("_pdbx_struct_assembly_gen.");
                input.push_str(tag);
                input.push('\n');
            }
            for row in generator_rows {
                input.push_str(&row.join(" "));
                input.push('\n');
            }
        }
        if !operator_rows.is_empty() {
            input.push_str("loop_\n");
            for tag in super::mmcif_transform_tags("matrix", "vector") {
                input.push_str("_pdbx_struct_oper_list.");
                input.push_str(&tag);
                input.push('\n');
            }
            input.push_str("_pdbx_struct_oper_list.id\n_pdbx_struct_oper_list.type\n");
            for (name, kind, translation) in operator_rows {
                let row = [
                    "1",
                    "0",
                    "0",
                    translation[0],
                    "0",
                    "1",
                    "0",
                    translation[1],
                    "0",
                    "0",
                    "1",
                    translation[2],
                    *name,
                    *kind,
                ];
                input.push_str(&row.join(" "));
                input.push('\n');
            }
        }
        input.push_str("loop_\n");
        for tag in [
            "id",
            "details",
            "method_details",
            "oligomeric_details",
            "oligomeric_count",
        ] {
            input.push_str("_pdbx_struct_assembly.");
            input.push_str(tag);
            input.push('\n');
        }
        for row in assembly_rows {
            input.push_str(&row.join(" "));
            input.push('\n');
        }
        document(&input, "c29-assemblies.cif")
    }

    fn c11_atom_site_row(
        id: &'static str,
        label_asym_id: &'static str,
        entity_id: &'static str,
        label_seq_id: &'static str,
        auth_seq_id: &'static str,
        auth_asym_id: &'static str,
        model_number: &'static str,
    ) -> [&'static str; 24] {
        let mut row = c05_scalar_row();
        row[0] = id;
        row[6] = label_asym_id;
        row[7] = entity_id;
        row[8] = label_seq_id;
        row[16] = auth_seq_id;
        row[18] = auth_asym_id;
        row[20] = model_number;
        row
    }

    fn c11_atom_site_grouping() -> MmcifAtomSiteGrouping {
        let rows = [
            c11_atom_site_row("1", "x", "e1", "1", "1", "A", "2"),
            c11_atom_site_row("2", "x", "e2", "2", "2", "A", "2"),
            c11_atom_site_row("3", "y", "e2", "3", "3", "A", "2"),
            c11_atom_site_row("4", "x", "e1", "4", "4", "A", "2"),
            c11_atom_site_row("5", "u", "missing", "5", "5", "A", "2"),
            c11_atom_site_row("6", "z", "e1", "6", "6", "A", "2"),
            c11_atom_site_row("7", "x", "e2", "1", "1", "B", "2"),
            c11_atom_site_row("8", "x", "e2", "2", "2", "B", "2"),
            c11_atom_site_row("9", "q", "e2", "1", "1", "C", "1"),
        ];
        let doc = c05_scalar_document(&rows, &[]);
        group_mmcif_atom_sites(&atom_site_table(&doc.blocks()[0])).unwrap()
    }

    fn c12_atom_site_grouping() -> MmcifAtomSiteGrouping {
        let rows = [
            c11_atom_site_row("1", "x", "e1", "1", "1", "A", "2"),
            c11_atom_site_row("2", "x", "e1", "2", "2", "A", "2"),
            c11_atom_site_row("3", "x", "e1", "3", "3", "A", "2"),
            c11_atom_site_row("4", "x", "e1", "1", "1", "A", "1"),
            c11_atom_site_row("5", "x", "e1", "2", "2", "A", "1"),
            c11_atom_site_row("6", "x", "e1", "3", "3", "A", "1"),
        ];
        let doc = c05_scalar_document(&rows, &[]);
        group_mmcif_atom_sites(&atom_site_table(&doc.blocks()[0])).unwrap()
    }

    fn c12_sifts_document(rows: &[&str], source: &str) -> crate::cif::CifDocument {
        let mut input = String::from(
            "data_c12_sifts\n_entity.id e1\n_entity.type polymer\nloop_\n_struct_asym.id\n_struct_asym.entity_id\nx e1\nloop_\n_pdbx_sifts_xref_db.entity_id\n_pdbx_sifts_xref_db.asym_id\n_pdbx_sifts_xref_db.seq_id_ordinal\n_pdbx_sifts_xref_db.seq_id\n_pdbx_sifts_xref_db.observed\n_pdbx_sifts_xref_db.unp_res\n_pdbx_sifts_xref_db.unp_num\n_pdbx_sifts_xref_db.unp_acc\n",
        );
        for row in rows {
            input.push_str(row);
            input.push('\n');
        }
        document(&input, source)
    }

    fn c12_parse_error(row: &str) -> MmcifEntityInfoError {
        let doc = c12_sifts_document(&[row], "c12-sifts-error.cif");
        let mut grouping = c12_atom_site_grouping();
        parse_mmcif_entity_polymer_values(&doc.blocks()[0], Some(&mut grouping)).unwrap_err()
    }

    fn c11_subchain_values(entities: &[BioEntityRow]) -> Vec<Vec<&str>> {
        entities
            .iter()
            .map(|entity| entity.subchains().iter().map(String::as_str).collect())
            .collect()
    }

    fn c06_anisotropic_document(rows: &[[&str; 7]], source: &str) -> crate::cif::CifDocument {
        let mut input = String::from("data_c06_anisotropic\nloop_\n");
        for tag in [
            "id", "U[1][1]", "U[2][2]", "U[3][3]", "U[1][2]", "U[1][3]", "U[2][3]",
        ] {
            input.push_str("_atom_site_anisotrop.");
            input.push_str(tag);
            input.push('\n');
        }
        for row in rows {
            for value in row {
                input.push_str(value);
                input.push(' ');
            }
            input.push('\n');
        }
        document(&input, source)
    }

    fn project_c05_scalar_row(
        values: [&str; 24],
        omitted_columns: &[usize],
        has_d_fraction: bool,
    ) -> Result<super::MmcifAtomSiteScalars, MmcifAtomSiteScalarError> {
        let document = c05_scalar_document(&[values], omitted_columns);
        let table = atom_site_table(&document.blocks()[0]);
        let row = table.row(0).expect("fixed one-row C05 atom_site table");
        project_mmcif_atom_site_scalars(row, has_d_fraction)
    }

    #[test]
    fn bio_read_c01_factored_document_matches_text_entry_and_source_metadata() {
        // Gemmi mmcif.hpp::make_structure selects block zero after checking
        // later blocks. A chemcomp-only block still follows ordinary mmCIF
        // population here, rather than the separate chemcomp coordinate path.
        let text = concat!(
            "data_component_only\n",
            "_entry.id DEMO\n",
            "_chem_comp.id LIG\n",
            "_chem_comp.three_letter_code LIG\n",
            "data_restraints\n",
            "_chem_comp.id OTHER\n",
        );
        let parsed = document(text, "source.cif");
        let direct = populate_mmcif_bio_structure_document(&parsed).unwrap();
        let text_entry = read_mmcif_bio_structure(text, "source.cif").unwrap();
        assert_eq!(direct, text_entry);
        assert_eq!(direct.source_state().name, "component_only");
        assert_eq!(
            direct.input_format(),
            cosmolkit_bio::BioCoordinateFormat::Mmcif
        );
        assert!(direct.atoms().is_empty());
    }

    #[test]
    fn bio_read_c01_empty_and_later_coordinate_errors_preserve_stage_and_source() {
        let empty = read_mmcif_bio_structure("", "empty.cif").unwrap_err();
        assert_eq!(empty.stage(), super::BioMmcifReadStage::CifDocument);
        let text = "data_first\n_entry.id first\ndata_second\n_atom_site.id 1\n";
        let parsed = document(text, "later.cif");
        let direct = populate_mmcif_bio_structure_document(&parsed).unwrap_err();
        let wrapped = read_mmcif_bio_structure(text, "later.cif").unwrap_err();
        assert_eq!(direct.stage(), super::BioMmcifReadStage::CoordinateBlock);
        assert_eq!(direct.stage(), wrapped.stage());
        assert_eq!(direct.to_string(), wrapped.to_string());
        assert!(direct.to_string().contains("block #2: later.cif"));
    }

    #[test]
    fn bio_mmcif_c01_empty_document_retains_out_of_range_category() {
        let error = select_coordinate_block(&[], "empty.cif").unwrap_err();
        assert_eq!(
            error,
            CoordinateBlockSelectionError::EmptyDocument {
                source: "empty.cif".to_owned(),
            }
        );
    }

    #[test]
    fn bio_mmcif_c01_selects_first_block_before_later_restraint_blocks() {
        let doc = document(
            concat!(
                "data_first\n",
                "_atom_site.id 1\n",
                "data_restraints\n",
                "_chem_comp.id COMP\n",
                "data_dictionary\n",
                "_item.name comp_id\n",
            ),
            "restraints.cif",
        );

        let selected = select_coordinate_block(doc.blocks(), doc.source()).unwrap();
        assert_eq!(selected.name(), "first");
    }

    #[test]
    fn bio_mmcif_c01_later_pair_tag_rejects_with_source_block_number() {
        let doc = document(
            concat!(
                "data_first\n",
                "_atom_site.id 1\n",
                "data_restraints\n",
                "_chem_comp.id COMP\n",
                "data_later\n",
                "_ATOM_SITE.ID .\n",
            ),
            "later-pair.cif",
        );

        let error = select_coordinate_block(doc.blocks(), doc.source()).unwrap_err();
        assert_eq!(
            error,
            CoordinateBlockSelectionError::LaterCoordinateBlock {
                block_number: 3,
                source: "later-pair.cif".to_owned(),
            }
        );
        assert_eq!(
            error.to_string(),
            "2+ blocks are ok if only the first one has coordinates;\n_atom_site in block #3: later-pair.cif"
        );
    }

    #[test]
    fn bio_mmcif_c01_later_loop_tag_presence_rejects_null_value() {
        let doc = document(
            concat!(
                "data_first\n",
                "_atom_site.id 1\n",
                "data_later\n",
                "loop_\n",
                "_atom_site.id\n",
                ".\n",
            ),
            "later-loop.cif",
        );

        let error = select_coordinate_block(doc.blocks(), doc.source()).unwrap_err();
        assert_eq!(
            error,
            CoordinateBlockSelectionError::LaterCoordinateBlock {
                block_number: 2,
                source: "later-loop.cif".to_owned(),
            }
        );
        assert!(error.to_string().ends_with("block #2: later-loop.cif"));
    }

    #[test]
    fn bio_mmcif_c02_auth_values_are_primary_with_per_row_label_fallback() {
        let doc = document(
            concat!(
                "data_auth_fallback\n",
                "loop_\n",
                "_atom_site.id\n",
                "_atom_site.type_symbol\n",
                "_atom_site.label_atom_id\n",
                "_atom_site.label_alt_id\n",
                "_atom_site.label_comp_id\n",
                "_atom_site.label_asym_id\n",
                "_atom_site.label_seq_id\n",
                "_atom_site.Cartn_x\n",
                "_atom_site.Cartn_y\n",
                "_atom_site.Cartn_z\n",
                "_atom_site.auth_seq_id\n",
                "_atom_site.auth_comp_id\n",
                "_atom_site.auth_asym_id\n",
                "_atom_site.auth_atom_id\n",
                "1 C L_ATOM0 . L_COMP0 L_ASYM0 L_SEQ0 0 0 0 A_SEQ0 A_COMP0 A_ASYM0 A_ATOM0\n",
                "2 C L_ATOM1 . L_COMP1 L_ASYM1 L_SEQ1 0 0 0 . . ? ?\n",
                "3 C L_ATOM2 . L_COMP2 L_ASYM2 L_SEQ2 0 0 0 ? ? . .\n",
                "4 C ? . . ? . 0 0 0 . ? ? .\n",
            ),
            "auth-fallback.cif",
        );
        let table = atom_site_table(&doc.blocks()[0]);
        let identities = bind_atom_site_identity_columns(&table).unwrap();

        let auth_row = table.row(0).unwrap();
        assert_eq!(identities.asym_id.get(auth_row), "A_ASYM0");
        assert_eq!(identities.comp_id.get(auth_row), "A_COMP0");
        assert_eq!(identities.atom_id.get(auth_row), "A_ATOM0");
        assert_eq!(identities.seq_id.get(auth_row), "A_SEQ0");

        let null_primary_row = table.row(1).unwrap();
        assert_eq!(identities.asym_id.get(null_primary_row), "L_ASYM1");
        assert_eq!(identities.comp_id.get(null_primary_row), "L_COMP1");
        assert_eq!(identities.atom_id.get(null_primary_row), "L_ATOM1");
        assert_eq!(identities.seq_id.get(null_primary_row), "L_SEQ1");

        let alternate_null_row = table.row(2).unwrap();
        assert_eq!(identities.asym_id.get(alternate_null_row), "L_ASYM2");
        assert_eq!(identities.comp_id.get(alternate_null_row), "L_COMP2");
        assert_eq!(identities.atom_id.get(alternate_null_row), "L_ATOM2");
        assert_eq!(identities.seq_id.get(alternate_null_row), "L_SEQ2");

        let null_fallback_row = table.row(3).unwrap();
        assert_eq!(identities.asym_id.get(null_fallback_row), "?");
        assert_eq!(identities.comp_id.get(null_fallback_row), ".");
        assert_eq!(identities.atom_id.get(null_fallback_row), "?");
        assert_eq!(identities.seq_id.get(null_fallback_row), ".");
    }

    #[test]
    fn bio_mmcif_c02_label_only_values_remain_primary_including_nulls() {
        let doc = document(
            concat!(
                "data_label_only\n",
                "loop_\n",
                "_atom_site.id\n",
                "_atom_site.type_symbol\n",
                "_atom_site.label_atom_id\n",
                "_atom_site.label_alt_id\n",
                "_atom_site.label_comp_id\n",
                "_atom_site.label_asym_id\n",
                "_atom_site.label_seq_id\n",
                "_atom_site.Cartn_x\n",
                "_atom_site.Cartn_y\n",
                "_atom_site.Cartn_z\n",
                "1 C ? . . LABEL_ASYM . 0 0 0\n",
            ),
            "label-only.cif",
        );
        let table = atom_site_table(&doc.blocks()[0]);
        let identities = bind_atom_site_identity_columns(&table).unwrap();
        let row = table.row(0).unwrap();

        assert_eq!(identities.asym_id.get(row), "LABEL_ASYM");
        assert_eq!(identities.comp_id.get(row), ".");
        assert_eq!(identities.atom_id.get(row), "?");
        assert_eq!(identities.seq_id.get(row), ".");
    }

    #[test]
    fn bio_mmcif_c02_missing_optional_identity_pairs_fail_in_source_order() {
        let missing_comp = document(
            concat!(
                "data_missing_comp\n",
                "loop_\n",
                "_atom_site.id\n",
                "_atom_site.type_symbol\n",
                "_atom_site.label_alt_id\n",
                "_atom_site.label_asym_id\n",
                "_atom_site.label_atom_id\n",
                "_atom_site.label_seq_id\n",
                "_atom_site.Cartn_x\n",
                "_atom_site.Cartn_y\n",
                "_atom_site.Cartn_z\n",
                "1 C . ASYM ATOM 1 0 0 0\n",
            ),
            "missing-comp.cif",
        );
        let table = atom_site_table(&missing_comp.blocks()[0]);
        let error = bind_atom_site_identity_columns(&table).unwrap_err();
        assert_eq!(
            error,
            MissingAtomSiteIdentityField(AtomSiteIdentityField::CompId)
        );
        assert_eq!(
            error.to_string(),
            "Neither _atom_site.label_comp_id nor auth_comp_id found"
        );

        let missing_atom = document(
            concat!(
                "data_missing_atom\n",
                "loop_\n",
                "_atom_site.id\n",
                "_atom_site.type_symbol\n",
                "_atom_site.label_alt_id\n",
                "_atom_site.label_asym_id\n",
                "_atom_site.label_comp_id\n",
                "_atom_site.label_seq_id\n",
                "_atom_site.Cartn_x\n",
                "_atom_site.Cartn_y\n",
                "_atom_site.Cartn_z\n",
                "1 C . ASYM COMP 1 0 0 0\n",
            ),
            "missing-atom.cif",
        );
        let table = atom_site_table(&missing_atom.blocks()[0]);
        let error = bind_atom_site_identity_columns(&table).unwrap_err();
        assert_eq!(
            error,
            MissingAtomSiteIdentityField(AtomSiteIdentityField::AtomId)
        );
        assert_eq!(
            error.to_string(),
            "Neither _atom_site.label_atom_id nor auth_atom_id found"
        );

        let missing_seq = document(
            concat!(
                "data_missing_seq\n",
                "loop_\n",
                "_atom_site.id\n",
                "_atom_site.type_symbol\n",
                "_atom_site.label_alt_id\n",
                "_atom_site.label_asym_id\n",
                "_atom_site.label_comp_id\n",
                "_atom_site.label_atom_id\n",
                "_atom_site.Cartn_x\n",
                "_atom_site.Cartn_y\n",
                "_atom_site.Cartn_z\n",
                "1 C . ASYM COMP ATOM 0 0 0\n",
            ),
            "missing-seq.cif",
        );
        let table = atom_site_table(&missing_seq.blocks()[0]);
        let error = bind_atom_site_identity_columns(&table).unwrap_err();
        assert_eq!(
            error,
            MissingAtomSiteIdentityField(AtomSiteIdentityField::SeqId)
        );
        assert_eq!(
            error.to_string(),
            "Neither _atom_site.label_seq_id nor auth_seq_id found"
        );
    }

    #[test]
    fn bio_mmcif_c02_missing_both_asym_columns_skips_source_identity_checks() {
        let doc = document(
            concat!(
                "data_missing_asym\n",
                "loop_\n",
                "_atom_site.id\n",
                "_atom_site.type_symbol\n",
                "_atom_site.label_alt_id\n",
                "_atom_site.auth_asym_id\n",
                "_atom_site.Cartn_x\n",
                "_atom_site.Cartn_y\n",
                "_atom_site.Cartn_z\n",
                "1 C . AUTH_ASYM 0 0 0\n",
            ),
            "missing-asym.cif",
        );
        let table = atom_site_table(&doc.blocks()[0]);

        // Gemmi's find() requires label_asym_id, and read_atom_sites skips
        // RowAccess construction when the resulting table has zero length.
        assert_eq!(table.len(), 0);
    }

    #[test]
    fn bio_mmcif_c03_embedded_and_separate_insertion_codes_match_gemmi() {
        let expected = PdbSeqId::new(15, Some(b'A'));
        assert_eq!(make_mmcif_seq_id("15A".to_owned(), None), Ok(expected));
        assert_eq!(
            make_mmcif_seq_id("15".to_owned(), Some("'A'")),
            Ok(expected)
        );
        assert_eq!(make_mmcif_seq_id("15A".to_owned(), Some("A")), Ok(expected));

        assert_eq!(
            make_mmcif_seq_id("A".to_owned(), None),
            Ok(PdbSeqId::new(i32::MIN, Some(b'A')))
        );
        assert_eq!(
            make_mmcif_seq_id("".to_owned(), None),
            Ok(PdbSeqId::new(i32::MIN, None))
        );
        assert_eq!(
            make_mmcif_seq_id("15".to_owned(), Some(".")),
            Ok(PdbSeqId::new(15, None))
        );
        assert_eq!(
            make_mmcif_seq_id("15".to_owned(), Some("?")),
            Ok(PdbSeqId::new(15, None))
        );
        assert_eq!(
            make_mmcif_seq_id("15".to_owned(), Some("' '")),
            Ok(PdbSeqId::new(15, None))
        );

        assert_eq!(
            make_mmcif_seq_id("15A".to_owned(), Some("B")),
            Err(MmcifSequenceIdError::InconsistentInsertionCode {
                sequence_id: "15A".to_owned(),
            })
        );
    }

    #[test]
    fn bio_mmcif_c03_literal_suffix_threshold_and_integer_boundaries_match_gemmi() {
        assert_eq!(
            make_mmcif_seq_id("15a".to_owned(), None),
            Ok(PdbSeqId::new(15, Some(b'a')))
        );
        assert_eq!(
            make_mmcif_seq_id("15_".to_owned(), None),
            Ok(PdbSeqId::new(15, Some(b'_')))
        );
        assert_eq!(
            make_mmcif_seq_id("+15".to_owned(), None),
            Ok(PdbSeqId::new(15, None))
        );
        assert_eq!(
            make_mmcif_seq_id("-15".to_owned(), None),
            Ok(PdbSeqId::new(-15, None))
        );
        assert_eq!(
            make_mmcif_seq_id(" \t15 \r\n".to_owned(), None),
            Ok(PdbSeqId::new(15, None))
        );
        assert_eq!(
            make_mmcif_seq_id("-2147483648".to_owned(), None),
            Ok(PdbSeqId::new(i32::MIN, None))
        );
        assert_eq!(
            make_mmcif_seq_id("2147483647".to_owned(), None),
            Ok(PdbSeqId::new(i32::MAX, None))
        );

        assert!(matches!(
            make_mmcif_seq_id("15A ".to_owned(), None),
            Err(MmcifSequenceIdError::InvalidSequenceNumber(error))
                if error.kind() == CifReadErrorKind::InvalidValue
        ));
        assert!(matches!(
            make_mmcif_seq_id("x?".to_owned(), None),
            Err(MmcifSequenceIdError::InvalidSequenceNumber(error))
                if error.kind() == CifReadErrorKind::InvalidValue
        ));
        assert!(matches!(
            make_mmcif_seq_id("2147483648".to_owned(), None),
            Err(MmcifSequenceIdError::InvalidSequenceNumber(error))
                if error.kind() == CifReadErrorKind::OutOfRange
        ));
        assert!(matches!(
            make_mmcif_seq_id("15".to_owned(), Some("'AB'")),
            Err(MmcifSequenceIdError::InvalidInsertionCode(error))
                if error.kind() == CifReadErrorKind::InvalidValue
        ));
    }

    #[test]
    fn bio_mmcif_c03_empty_decoded_insertion_uses_gemmi_string_terminator() {
        // The pinned C++ source indexes std::string at pos == size() in both
        // branches of cif::as_char. C++ basic_string guarantees the
        // terminating charT() at data() + size(); this is defined NUL, not UB.
        for raw in ["", "''", "\"\"", ";\n;"] {
            assert_eq!(
                make_mmcif_seq_id("15".to_owned(), Some(raw)),
                Ok(PdbSeqId::new(15, Some(0))),
                "raw insertion token {raw:?}"
            );
        }
    }

    #[test]
    fn bio_mmcif_c03_residue_address_maps_gemmi_sentinels_canonically() {
        let name = ResidueName::from_ascii(b"GLY").expect("valid fixed residue name");
        let address = make_mmcif_residue_address(name, PdbSeqId::new(i32::MIN, Some(b'A')))
            .expect("empty segment is representable");
        assert_eq!(address.sequence_number(), None);
        assert_eq!(address.insertion_code(), Some(b'A'));
        assert_eq!(address.segment(), "");
        assert_eq!(address.name(), name);

        let nul_address = make_mmcif_residue_address(name, PdbSeqId::new(15, Some(0)))
            .expect("empty segment is representable");
        assert_eq!(nul_address.sequence_number(), Some(15));
        assert_eq!(nul_address.insertion_code(), Some(0));
    }

    #[test]
    fn bio_mmcif_c04_raw_model_transitions_reuse_numeric_models_in_row_order() {
        let doc = document(
            concat!(
                "data_model_transitions\n",
                "loop_\n",
                "_atom_site.id\n",
                "_atom_site.group_PDB\n",
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
                "_atom_site.auth_comp_id\n",
                "_atom_site.auth_asym_id\n",
                "_atom_site.auth_atom_id\n",
                "_atom_site.pdbx_PDB_model_num\n",
                "1 ATOM C CA . GLY L1 E1 100 . 0 0 0 1 20 0 10 GLY A CA 01\n",
                "2 HETATM C CB . GLY L2 E9 999 . 1 2 3 1 30 0 10 GLY A CB 01\n",
                "3 HETATM C CG . GLY L2 E2 101 . 2 3 4 1 40 0 11 GLY A CG 01\n",
                "4 HETATM C CA . GLY L3 E3 102 . 3 4 5 1 50 0 10 GLY A CA 1\n",
                "5 ATOM C CA . GLY L4 E4 200 . 4 5 6 1 60 0 20 GLY B CA 2\n",
                "6 ATOM C CA . GLY L5 E5 103 . 5 6 7 1 70 0 10 GLY A CA 1\n",
            ),
            "model-transitions.cif",
        );
        let grouping = group_mmcif_atom_sites(&atom_site_table(&doc.blocks()[0])).unwrap();

        assert_eq!(
            grouping
                .models
                .iter()
                .map(|model| model.source_model_number)
                .collect::<Vec<_>>(),
            [1, 2]
        );
        let first_model = &grouping.models[0];
        assert_eq!(first_model.chains.len(), 3);
        assert_eq!(
            first_model
                .chains
                .iter()
                .map(|chain| chain.source_name.as_str())
                .collect::<Vec<_>>(),
            ["A", "A", "A"]
        );
        assert_eq!(
            first_model.chains[0].source.auth_chain_id(),
            Some(PdbChainId::from_ascii(b"A").unwrap())
        );
        assert_eq!(first_model.chains[0].source.label_asym_id(), Some("L1"));
        assert_eq!(first_model.chains[0].residues.len(), 2);

        let first_residue = &first_model.chains[0].residues[0];
        assert_eq!(first_residue.atom_site_rows, [0, 1]);
        assert_eq!(first_residue.address.sequence_number(), Some(10));
        assert_eq!(first_residue.source.label_seq_id(), Some(100));
        assert_eq!(first_residue.source.subchain_id(), Some("L1"));
        assert_eq!(first_residue.source.label_entity_id(), Some("E1"));
        assert_eq!(first_residue.het_flag, Some(b'A'));

        let second_residue = &first_model.chains[0].residues[1];
        assert_eq!(second_residue.atom_site_rows, [2]);
        assert_eq!(second_residue.address.sequence_number(), Some(11));
        assert_eq!(second_residue.source.label_seq_id(), Some(101));
        assert_eq!(second_residue.source.subchain_id(), Some("L2"));
        assert_eq!(second_residue.source.label_entity_id(), Some("E2"));

        assert_eq!(first_model.chains[1].source.label_asym_id(), Some("L3"));
        assert_eq!(first_model.chains[1].residues[0].atom_site_rows, [3]);
        assert_eq!(first_model.chains[1].residues[0].het_flag, Some(b'H'));
        assert_eq!(first_model.chains[2].source.label_asym_id(), Some("L5"));
        assert_eq!(first_model.chains[2].residues[0].atom_site_rows, [5]);

        let second_model = &grouping.models[1];
        assert_eq!(second_model.chains.len(), 1);
        assert_eq!(
            second_model.chains[0]
                .source
                .auth_chain_id()
                .unwrap()
                .as_str(),
            "B"
        );
        assert_eq!(second_model.chains[0].residues[0].atom_site_rows, [4]);
    }

    #[test]
    fn bio_mmcif_c04_absent_model_column_creates_model_one_and_keeps_first_metadata() {
        let doc = document(
            concat!(
                "data_default_model\n",
                "loop_\n",
                "_atom_site.id\n",
                "_atom_site.group_PDB\n",
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
                "_atom_site.auth_comp_id\n",
                "_atom_site.auth_asym_id\n",
                "_atom_site.auth_atom_id\n",
                "1 ATOM C CA . GLY L1 E1 100 . 0 0 0 1 20 0 10 GLY A CA\n",
                "2 HETATM C CB . GLY L2 E2 999 . 1 2 3 1 30 0 10 GLY A CB\n",
            ),
            "default-model.cif",
        );
        let grouping = group_mmcif_atom_sites(&atom_site_table(&doc.blocks()[0])).unwrap();

        assert_eq!(grouping.models.len(), 1);
        assert_eq!(grouping.models[0].source_model_number, 1);
        assert_eq!(grouping.models[0].chains.len(), 1);
        let chain = &grouping.models[0].chains[0];
        assert_eq!(chain.source_name, "A");
        assert_eq!(chain.source.label_asym_id(), Some("L1"));
        assert_eq!(chain.residues.len(), 1);
        let residue = &chain.residues[0];
        assert_eq!(residue.atom_site_rows, [0, 1]);
        assert_eq!(residue.source.label_seq_id(), Some(100));
        assert_eq!(residue.source.subchain_id(), Some("L1"));
        assert_eq!(residue.source.label_entity_id(), Some("E1"));
        assert_eq!(residue.het_flag, Some(b'A'));
    }

    #[test]
    fn bio_mmcif_c05_serial_null_defaults_and_nan_position_match_gemmi() {
        let mut values = c05_scalar_row();
        values[0] = "not-a-number";
        values[10] = "bad";
        values[11] = "-0.0";
        values[12] = "1.234567890123";
        let scalars = project_c05_scalar_row(values, &[], true).unwrap();

        assert_eq!(scalars.serial.value(), 0);
        assert_eq!(scalars.element, Element::C);
        assert_eq!(scalars.isotope_mass_number, None);
        assert_eq!(scalars.altloc, None);
        assert_eq!(scalars.formal_charge, 0);
        assert!(scalars.position[0].is_nan());
        assert_eq!(scalars.position[1].to_bits(), 0x8000_0000_0000_0000);
        assert_eq!(scalars.position[2].to_bits(), 0x3ff3_c0ca_428c_51f2);
        assert_eq!(scalars.occupancy, 1.0);
        assert_eq!(scalars.b_iso, 20.0);
        assert_eq!(scalars.calc_flag, BioCalcFlag::NotSet);
        assert_eq!(scalars.tls_group_id, -1);
        assert_eq!(scalars.fraction, 0.0);

        let missing_optional_columns =
            project_c05_scalar_row(values, &[13, 14, 15, 21, 22, 23], false).unwrap();
        assert_eq!(missing_optional_columns.formal_charge, 0);
        assert_eq!(missing_optional_columns.occupancy, 1.0);
        assert_eq!(missing_optional_columns.b_iso, 20.0);
        assert_eq!(missing_optional_columns.calc_flag, BioCalcFlag::NotSet);
        assert_eq!(missing_optional_columns.tls_group_id, -1);
        assert_eq!(missing_optional_columns.fraction, 0.0);
    }

    #[test]
    fn bio_mmcif_c05_scalar_values_match_pinned_gemmi_assignments() {
        let mut values = c05_scalar_row();
        values[0] = "+17tail";
        values[2] = "D";
        values[4] = "'A'";
        values[10] = "1.234567890123";
        values[11] = "-0.0";
        values[12] = "bad";
        values[13] = "0.123456789";
        values[14] = "0.123456789";
        values[15] = "128";
        values[21] = "dummy";
        values[22] = "12tail";
        values[23] = "0.123456789";
        let scalars = project_c05_scalar_row(values, &[], true).unwrap();

        assert_eq!(scalars.serial.value(), 17);
        assert_eq!(scalars.element, Element::H);
        assert_eq!(scalars.isotope_mass_number, Some(2));
        assert_eq!(scalars.altloc, Some(AltLocLabel::new(b'A')));
        // The pinned C++ compiler converts the source signed-char assignment
        // 128 to -128; this is not a portable C++ guarantee.
        assert_eq!(scalars.formal_charge, -128);
        assert_eq!(scalars.position[0].to_bits(), 0x3ff3_c0ca_428c_51f2);
        assert_eq!(scalars.position[1].to_bits(), 0x8000_0000_0000_0000);
        assert!(scalars.position[2].is_nan());
        assert_eq!(scalars.occupancy.to_bits(), 0x3dfc_d6ea);
        assert_eq!(scalars.b_iso.to_bits(), 0x3dfc_d6ea);
        assert_eq!(scalars.calc_flag, BioCalcFlag::Dummy);
        assert_eq!(scalars.tls_group_id, 12);
        assert_eq!(scalars.fraction.to_bits(), 0x3dfc_d6ea);

        values[15] = "+7";
        assert_eq!(
            project_c05_scalar_row(values, &[], true)
                .unwrap()
                .formal_charge,
            7
        );
        values[15] = "-7";
        assert_eq!(
            project_c05_scalar_row(values, &[], true)
                .unwrap()
                .formal_charge,
            -7
        );

        values[2] = "H";
        let hydrogen = project_c05_scalar_row(values, &[], true).unwrap();
        assert_eq!(hydrogen.element, Element::H);
        assert_eq!(hydrogen.isotope_mass_number, None);

        values[2] = "Qq";
        let unknown = project_c05_scalar_row(values, &[], true).unwrap();
        assert_eq!(unknown.element, Element::DUMMY);
        assert_eq!(unknown.isotope_mass_number, None);
    }

    #[test]
    fn bio_mmcif_c05_calc_flag_tls_and_fraction_keep_source_branches() {
        let mut values = c05_scalar_row();
        values[21] = "calculated";
        values[22] = "12tail";
        let scalars = project_c05_scalar_row(values, &[], false).unwrap();
        assert_eq!(scalars.calc_flag, BioCalcFlag::Calculated);
        assert_eq!(scalars.tls_group_id, 12);
        assert_eq!(scalars.fraction, 0.0);

        values[21] = "d?";
        values[22] = "+12";
        values[23] = "0.75";
        let scalars = project_c05_scalar_row(values, &[], false).unwrap();
        assert_eq!(scalars.calc_flag, BioCalcFlag::Determined);
        assert_eq!(scalars.tls_group_id, -1);
        assert_eq!(scalars.fraction, 0.0);

        values[21] = "DUMMY";
        values[22] = "40000";
        values[23] = ".";
        let scalars = project_c05_scalar_row(values, &[], true).unwrap();
        assert_eq!(scalars.calc_flag, BioCalcFlag::NotSet);
        // The pinned C++ compiler narrows source short 40000 to -25536.
        assert_eq!(scalars.tls_group_id, -25536);
        assert_eq!(scalars.fraction, 0.0);

        let mut absent_fraction = c05_scalar_row();
        absent_fraction[23] = "0.75";
        let scalars = project_c05_scalar_row(absent_fraction, &[23], false).unwrap();
        assert_eq!(scalars.fraction, 0.0);
    }

    #[test]
    fn bio_mmcif_c05_source_conversion_errors_keep_order_and_categories() {
        let mut values = c05_scalar_row();
        values[4] = "'AB'";
        values[15] = "7tail";
        let error = project_c05_scalar_row(values, &[], false).unwrap_err();
        assert!(matches!(error, MmcifAtomSiteScalarError::CifCharacter(_)));

        values[4] = ".";
        let error = project_c05_scalar_row(values, &[], false).unwrap_err();
        assert!(matches!(error, MmcifAtomSiteScalarError::CifInteger(_)));
    }

    #[test]
    fn bio_mmcif_c05_source_undefined_integer_overflow_is_typed() {
        let mut values = c05_scalar_row();
        values[0] = "2147483648";
        let error = project_c05_scalar_row(values, &[], false).unwrap_err();
        assert!(matches!(
            error,
            MmcifAtomSiteScalarError::SerialOutsideSourceDefinedRange { .. }
        ));

        values[0] = "1";
        values[15] = "2147483648";
        let error = project_c05_scalar_row(values, &[], false).unwrap_err();
        assert!(matches!(error, MmcifAtomSiteScalarError::CifInteger(_)));

        values[15] = "0";
        values[22] = "2147483648";
        let error = project_c05_scalar_row(values, &[], false).unwrap_err();
        assert!(matches!(
            error,
            MmcifAtomSiteScalarError::TlsOutsideSourceDefinedRange { .. }
        ));
    }

    #[test]
    fn bio_mmcif_c06_anisotropic_map_uses_raw_ids_and_first_duplicate() {
        let doc = c06_anisotropic_document(
            &[
                ["17", "0.123456789", "2", "3", "4", "5", "6"],
                ["17", "9", "8", "7", "6", "5", "4"],
                ["0017", "-1", "-2", "-3", "-4", "-5", "-6"],
            ],
            "aniso-duplicates.cif",
        );
        let anisotropic = extract_mmcif_anisotropic_u(&doc.blocks()[0]).unwrap();
        assert_eq!(anisotropic.len(), 2);

        let first = anisotropic.get("17").unwrap();
        assert_eq!(
            first.map(|value| value.to_bits()),
            [
                0x3dfc_d6ea,
                0x4000_0000,
                0x4040_0000,
                0x4080_0000,
                0x40a0_0000,
                0x40c0_0000,
            ]
        );

        let mut atom_site = c05_scalar_row();
        atom_site[0] = "0017";
        let atom_doc = c05_scalar_document(&[atom_site], &[]);
        let atom_table = atom_site_table(&atom_doc.blocks()[0]);
        let atom_row = atom_table.row(0).expect("one atom-site row");
        let raw_atom_site_id = atom_site_raw_value(atom_row, AtomSiteColumn::Id.index());
        assert_eq!(raw_atom_site_id, "0017");
        assert_eq!(
            crate::bio_pdb::read_int(raw_atom_site_id.as_bytes()),
            Some(17)
        );
        assert_eq!(
            anisotropic.get(raw_atom_site_id).copied(),
            Some([-1.0, -2.0, -3.0, -4.0, -5.0, -6.0])
        );
        assert_eq!(
            anisotropic.get("17").copied(),
            Some([f32::from_bits(0x3dfc_d6ea), 2.0, 3.0, 4.0, 5.0, 6.0,])
        );
        assert!(anisotropic.get("17tail").is_none());
    }

    #[test]
    fn bio_mmcif_c06_anisotropic_map_handles_pair_tables_and_missing_columns() {
        let pair_doc = document(
            concat!(
                "data_aniso_pairs\n",
                "_atom_site_anisotrop.id 31\n",
                "_atom_site_anisotrop.U[1][1] 1.25\n",
                "_atom_site_anisotrop.U[2][2] 2.5\n",
                "_atom_site_anisotrop.U[3][3] 3.75\n",
                "_atom_site_anisotrop.U[1][2] 4.125\n",
                "_atom_site_anisotrop.U[1][3] 5.25\n",
                "_atom_site_anisotrop.U[2][3] 6.5\n",
            ),
            "aniso-pairs.cif",
        );
        let pair_map = extract_mmcif_anisotropic_u(&pair_doc.blocks()[0]).unwrap();
        assert_eq!(
            pair_map.get("31").copied(),
            Some([1.25, 2.5, 3.75, 4.125, 5.25, 6.5])
        );

        let missing_required = document(
            concat!(
                "data_aniso_missing_required\n",
                "_atom_site_anisotrop.id 31\n",
                "_atom_site_anisotrop.U[1][1] 1\n",
                "_atom_site_anisotrop.U[2][2] 2\n",
                "_atom_site_anisotrop.U[3][3] 3\n",
                "_atom_site_anisotrop.U[1][2] 4\n",
                "_atom_site_anisotrop.U[1][3] 5\n",
            ),
            "aniso-missing-column.cif",
        );
        assert!(
            extract_mmcif_anisotropic_u(&missing_required.blocks()[0])
                .unwrap()
                .is_empty()
        );
        let no_category = document("data_no_aniso\n", "no-aniso.cif");
        assert!(
            extract_mmcif_anisotropic_u(&no_category.blocks()[0])
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn bio_mmcif_c06_anisotropic_values_keep_source_nan_fallback() {
        let doc = c06_anisotropic_document(
            &[["42", ".", "?", "malformed", "4", "5", "6"]],
            "aniso-invalid-components.cif",
        );
        let anisotropic = extract_mmcif_anisotropic_u(&doc.blocks()[0]).unwrap();
        let values = anisotropic.get("42").unwrap();
        assert!(values[0].is_nan());
        assert!(values[1].is_nan());
        assert!(values[2].is_nan());
        assert_eq!(&values[3..], &[4.0, 5.0, 6.0]);
    }

    #[test]
    fn bio_mmcif_c07_fraction_column_presence_is_independent_of_zero_values() {
        let mut deuterium = c05_scalar_row();
        deuterium[2] = "D";
        deuterium[23] = "0";
        let mut hydrogen = c05_scalar_row();
        hydrogen[0] = "2";
        hydrogen[2] = "H";
        hydrogen[23] = "0.0";
        let doc = c05_scalar_document(&[deuterium, hydrogen], &[]);
        let table = atom_site_table(&doc.blocks()[0]);

        let has_d_fraction = has_mmcif_deuterium_fraction_column(&table);
        assert!(has_d_fraction);
        assert_eq!(table.len(), 2);

        let first = project_mmcif_atom_site_scalars(table.row(0).unwrap(), has_d_fraction).unwrap();
        let second =
            project_mmcif_atom_site_scalars(table.row(1).unwrap(), has_d_fraction).unwrap();
        assert_eq!(first.element, Element::H);
        assert_eq!(first.isotope_mass_number, Some(2));
        assert_eq!(first.fraction, 0.0);
        assert_eq!(second.element, Element::H);
        assert_eq!(second.isotope_mass_number, None);
        assert_eq!(second.fraction, 0.0);
    }

    #[test]
    fn bio_mmcif_c07_missing_fraction_column_keeps_default_state() {
        let mut values = c05_scalar_row();
        values[23] = "0.75";
        let doc = c05_scalar_document(&[values], &[23]);
        let table = atom_site_table(&doc.blocks()[0]);

        assert_eq!(table.len(), 1);
        assert!(!has_mmcif_deuterium_fraction_column(&table));
        let row = table.row(0).unwrap();
        assert_eq!(
            project_mmcif_atom_site_scalars(row, false)
                .unwrap()
                .fraction,
            0.0
        );

        let no_category = document("data_no_atom_site\n", "c07-no-atom-site.cif");
        assert!(!has_mmcif_deuterium_fraction_column(&atom_site_table(
            &no_category.blocks()[0]
        )));
    }

    #[test]
    fn bio_mmcif_c07_present_fraction_column_preserves_each_atom_value() {
        let mut deuterium = c05_scalar_row();
        deuterium[2] = "D";
        deuterium[23] = "0.75";
        let mut hydrogen = c05_scalar_row();
        hydrogen[0] = "2";
        hydrogen[2] = "H";
        hydrogen[23] = "0.25";
        let doc = c05_scalar_document(&[deuterium, hydrogen], &[]);
        let table = atom_site_table(&doc.blocks()[0]);
        let has_d_fraction = has_mmcif_deuterium_fraction_column(&table);
        assert!(has_d_fraction);

        let deuterium =
            project_mmcif_atom_site_scalars(table.row(0).unwrap(), has_d_fraction).unwrap();
        let hydrogen =
            project_mmcif_atom_site_scalars(table.row(1).unwrap(), has_d_fraction).unwrap();
        assert_eq!(deuterium.fraction.to_bits(), 0.75_f32.to_bits());
        assert_eq!(deuterium.element, Element::H);
        assert_eq!(deuterium.isotope_mass_number, Some(2));
        assert_eq!(hydrogen.fraction.to_bits(), 0.25_f32.to_bits());
        assert_eq!(hydrogen.element, Element::H);
        assert_eq!(hydrogen.isotope_mass_number, None);
    }

    #[test]
    fn bio_mmcif_c08_entity_and_polymer_types_match_pinned_source_order() {
        let doc = document(
            concat!(
                "data_c08_types\n",
                "loop_\n",
                "_entity.id\n",
                "_entity.type\n",
                "e0 .\n",
                "e1 water\n",
                "e2 non-polymer\n",
                "e3 branched\n",
                "e4 Polymer\n",
                "e5 ?\n",
                "e6 ?\n",
                "e7 ?\n",
                "e8 ?\n",
                "e9 ?\n",
                "e10 ?\n",
                "e11 ?\n",
                "e12 ?\n",
                "e13 ?\n",
                "e14 ?\n",
                "loop_\n",
                "_entity_poly.entity_id\n",
                "_entity_poly.type\n",
                "e0 \"polypeptide(L)\"\n",
                "e1 polyribonucleotide\n",
                "e2 \"polydeoxyribonucleotide/polyribonucleotide hybrid\"\n",
                "e3 \"polypeptide(D)\"\n",
                "e3 \"polypeptide(L)\"\n",
                "e4 Polypeptide\n",
                "e6 polydeoxyribonucleotide\n",
                "e7 polyribonucleotide\n",
                "e8 \"polydeoxyribonucleotide/polyribonucleotide hybrid\"\n",
                "e9 \"polypeptide(D)\"\n",
                "e10 \"polysaccharide(D)\"\n",
                "e11 other\n",
                "e12 \"peptide nucleic acid\"\n",
                "e13 cyclic-pseudo-peptide\n",
                "e14 \"polysaccharide(L)\"\n",
            ),
            "c08-entity-types.cif",
        );

        let entities = parse_mmcif_entity_polymer_values(&doc.blocks()[0], None).unwrap();
        let expected = [
            ("e0", EntityKind::Polymer, PolymerKind::PeptideL),
            ("e1", EntityKind::Water, PolymerKind::Rna),
            ("e2", EntityKind::NonPolymer, PolymerKind::DnaRnaHybrid),
            ("e3", EntityKind::Branched, PolymerKind::PeptideD),
            ("e4", EntityKind::Polymer, PolymerKind::Unknown),
            ("e5", EntityKind::Unknown, PolymerKind::Unknown),
            ("e6", EntityKind::Polymer, PolymerKind::Dna),
            ("e7", EntityKind::Polymer, PolymerKind::Rna),
            ("e8", EntityKind::Polymer, PolymerKind::DnaRnaHybrid),
            ("e9", EntityKind::Polymer, PolymerKind::PeptideD),
            ("e10", EntityKind::Polymer, PolymerKind::SaccharideD),
            ("e11", EntityKind::Polymer, PolymerKind::Other),
            ("e12", EntityKind::Polymer, PolymerKind::Pna),
            ("e13", EntityKind::Polymer, PolymerKind::CyclicPseudoPeptide),
            ("e14", EntityKind::Polymer, PolymerKind::SaccharideL),
        ];
        assert_eq!(entities.len(), expected.len());
        for (entity, (source_id, kind, polymer_kind)) in entities.iter().zip(expected) {
            assert_eq!(entity.source().source_entity_id(), source_id);
            assert_eq!(entity.kind(), kind, "entity {source_id}");
            assert_eq!(entity.polymer_kind(), polymer_kind, "entity {source_id}");
            assert!(entity.reflects_microhetero(), "entity {source_id}");
        }
    }

    #[test]
    fn bio_mmcif_c08_absent_type_and_polymer_row_keep_source_defaults() {
        let no_polymer = document(
            concat!(
                "data_c08_no_polymer\n",
                "loop_\n",
                "_entity.id\n",
                "_entity.type\n",
                "explicit-water water\n",
                "explicit-nonpolymer non-polymer\n",
                "explicit-unknown mystery\n",
            ),
            "c08-no-polymer.cif",
        );
        let no_polymer_entities =
            parse_mmcif_entity_polymer_values(&no_polymer.blocks()[0], None).unwrap();
        assert_eq!(no_polymer_entities.len(), 3);
        assert_eq!(
            no_polymer_entities
                .iter()
                .map(|entity| (entity.kind(), entity.polymer_kind()))
                .collect::<Vec<_>>(),
            [
                (EntityKind::Water, PolymerKind::Unknown),
                (EntityKind::NonPolymer, PolymerKind::Unknown),
                (EntityKind::Unknown, PolymerKind::Unknown),
            ]
        );

        let absent_type = document(
            concat!(
                "data_c08_absent_type\n",
                "loop_\n",
                "_entity.id\n",
                "id-only-a\n",
                "id-only-b\n",
            ),
            "c08-absent-type.cif",
        );
        let absent_type_entities =
            parse_mmcif_entity_polymer_values(&absent_type.blocks()[0], None).unwrap();
        assert_eq!(absent_type_entities.len(), 2);
        assert_eq!(
            absent_type_entities
                .iter()
                .map(|entity| (entity.kind(), entity.polymer_kind()))
                .collect::<Vec<_>>(),
            [
                (EntityKind::Unknown, PolymerKind::Unknown),
                (EntityKind::Unknown, PolymerKind::Unknown),
            ]
        );
    }

    #[test]
    fn bio_mmcif_c08_pair_tables_follow_source_entity_and_polymer_lookup() {
        let doc = document(
            concat!(
                "data_c08_pairs\n",
                "_entity.id pair-entity\n",
                "_entity.type ?\n",
                "_entity_poly.entity_id pair-entity\n",
                "_entity_poly.type \"polypeptide(L)\"\n",
            ),
            "c08-pair-tables.cif",
        );

        let entities = parse_mmcif_entity_polymer_values(&doc.blocks()[0], None).unwrap();
        assert_eq!(entities.len(), 1);
        assert_eq!(entities[0].source().source_entity_id(), "pair-entity");
        assert_eq!(entities[0].kind(), EntityKind::Polymer);
        assert_eq!(entities[0].polymer_kind(), PolymerKind::PeptideL);
        assert!(entities[0].reflects_microhetero());
    }

    #[test]
    fn bio_mmcif_c09_sequences_append_and_join_microheterogeneity_in_source_order() {
        let doc = document(
            concat!(
                "data_c09_ordered_sequence\n",
                "loop_\n",
                "_entity.id\n",
                "_entity.type\n",
                "e1 polymer\n",
                "e2 polymer\n",
                "loop_\n",
                "_entity_poly_seq.entity_id\n",
                "_entity_poly_seq.num\n",
                "_entity_poly_seq.mon_id\n",
                "e1 1 ALA\n",
                "e2 1 GLY\n",
                "e1 2 SER\n",
                "e1 2 THR\n",
                "e2 2 CYS\n",
            ),
            "c09-ordered-sequence.cif",
        );

        let entities = parse_mmcif_entity_polymer_values(&doc.blocks()[0], None).unwrap();
        assert_eq!(entities.len(), 2);
        assert_eq!(entities[0].source().source_entity_id(), "e1");
        assert_eq!(entities[0].full_sequence(), ["ALA", "SER,THR"]);
        assert_eq!(entities[1].source().source_entity_id(), "e2");
        assert_eq!(entities[1].full_sequence(), ["GLY", "CYS"]);
    }

    #[test]
    fn bio_mmcif_c09_sequence_null_negative_holes_and_missing_entity_follow_source() {
        let doc = document(
            concat!(
                "data_c09_sequence_gaps\n",
                "loop_\n",
                "_entity.id\n",
                "_entity.type\n",
                "e1 polymer\n",
                "e2 polymer\n",
                "loop_\n",
                "_entity_poly_seq.entity_id\n",
                "_entity_poly_seq.num\n",
                "_entity_poly_seq.mon_id\n",
                "e1 3 forward-hole\n",
                "missing 1 absent-entity\n",
                "e1 . dot-null\n",
                "e1 ? question-null\n",
                "e1 0 zero-position\n",
                "e1 -1 negative-position\n",
                "e1 1 FIRST\n",
                "e1 2147483647 int-max-hole\n",
                "e1 9 later-hole\n",
                "e1 2 SECOND\n",
            ),
            "c09-sequence-gaps.cif",
        );

        let entities = parse_mmcif_entity_polymer_values(&doc.blocks()[0], None).unwrap();
        assert_eq!(entities.len(), 2);
        assert_eq!(entities[0].full_sequence(), ["FIRST", "SECOND"]);
        assert!(entities[1].full_sequence().is_empty());
    }

    #[test]
    fn bio_mmcif_c09_sequence_invalid_integer_returns_structured_cif_error() {
        let doc = document(
            concat!(
                "data_c09_sequence_invalid_num\n",
                "_entity.id e1\n",
                "_entity.type polymer\n",
                "loop_\n",
                "_entity_poly_seq.entity_id\n",
                "_entity_poly_seq.num\n",
                "_entity_poly_seq.mon_id\n",
                "e1 not-an-integer ALA\n",
            ),
            "c09-invalid-sequence-number.cif",
        );

        let MmcifEntityInfoError::Cif(error) =
            parse_mmcif_entity_polymer_values(&doc.blocks()[0], None).unwrap_err()
        else {
            panic!("invalid C09 sequence integer must retain its CIF error")
        };
        assert_eq!(error.kind(), CifReadErrorKind::InvalidValue);
        assert_eq!(error.message(), "not an integer: not-an-integer");
    }

    #[test]
    fn bio_mmcif_c11_struct_asym_uses_source_order_and_preserves_duplicates() {
        let doc = document(
            concat!(
                "data_c11_direct_struct_asym\n",
                "loop_\n",
                "_entity.id\n",
                "_entity.type\n",
                "e1 polymer\n",
                "e2 polymer\n",
                "loop_\n",
                "_struct_asym.id\n",
                "_struct_asym.entity_id\n",
                "x e2\n",
                "unknown absent\n",
                "y e1\n",
                "x e1\n",
                "x e1\n",
            ),
            "c11-direct-struct-asym.cif",
        );
        let mut grouping = c11_atom_site_grouping();
        assert_eq!(grouping.models[0].source_model_number, 2);

        let entities =
            parse_mmcif_entity_polymer_values(&doc.blocks()[0], Some(&mut grouping)).unwrap();
        assert_eq!(
            entities
                .iter()
                .map(|entity| entity.source().source_entity_id())
                .collect::<Vec<_>>(),
            ["e1", "e2"]
        );
        assert_eq!(
            c11_subchain_values(&entities),
            [vec!["y", "x", "x"], vec!["x"]]
        );
    }

    #[test]
    fn bio_mmcif_c11_first_model_fallback_uses_contiguous_runs_and_entity_dedup() {
        let doc = document(
            concat!(
                "data_c11_fallback\n",
                "loop_\n",
                "_entity.id\n",
                "_entity.type\n",
                "e1 polymer\n",
                "e2 polymer\n",
            ),
            "c11-first-model-fallback.cif",
        );
        let mut grouping = c11_atom_site_grouping();
        assert_eq!(
            grouping
                .models
                .iter()
                .map(|model| model.source_model_number)
                .collect::<Vec<_>>(),
            [2, 1]
        );
        assert_eq!(grouping.models[0].chains.len(), 2);
        assert_eq!(
            grouping.models[0].chains[0]
                .residues
                .iter()
                .map(|residue| residue.source.subchain_id())
                .collect::<Vec<_>>(),
            [
                Some("x"),
                Some("x"),
                Some("y"),
                Some("x"),
                Some("u"),
                Some("z")
            ]
        );

        let entities =
            parse_mmcif_entity_polymer_values(&doc.blocks()[0], Some(&mut grouping)).unwrap();
        assert_eq!(
            c11_subchain_values(&entities),
            [vec!["x", "z"], vec!["y", "x"]]
        );
    }

    #[test]
    fn bio_mmcif_c11_present_empty_struct_asym_suppresses_fallback() {
        let doc = document(
            concat!(
                "data_c11_empty_struct_asym\n",
                "loop_\n",
                "_entity.id\n",
                "_entity.type\n",
                "e1 polymer\n",
                "e2 polymer\n",
                "loop_\n",
                "_struct_asym.id\n",
                "_struct_asym.entity_id\n",
            ),
            "c11-empty-struct-asym.cif",
        );
        let mut grouping = c11_atom_site_grouping();
        let entities =
            parse_mmcif_entity_polymer_values(&doc.blocks()[0], Some(&mut grouping)).unwrap();
        assert_eq!(
            c11_subchain_values(&entities),
            [Vec::<&str>::new(), Vec::new()]
        );
    }

    #[test]
    fn bio_mmcif_c11_absent_or_not_ok_struct_asym_uses_fallback_only_with_model() {
        let absent = document(
            concat!(
                "data_c11_absent_struct_asym\n",
                "loop_\n",
                "_entity.id\n",
                "_entity.type\n",
                "e1 polymer\n",
                "e2 polymer\n",
            ),
            "c11-absent-struct-asym.cif",
        );
        let not_ok = document(
            concat!(
                "data_c11_not_ok_struct_asym\n",
                "loop_\n",
                "_entity.id\n",
                "_entity.type\n",
                "e1 polymer\n",
                "e2 polymer\n",
                "loop_\n",
                "_struct_asym.id\n",
                "x\n",
            ),
            "c11-not-ok-struct-asym.cif",
        );
        let mut grouping = c11_atom_site_grouping();
        for block in [
            absent.blocks().first().unwrap(),
            not_ok.blocks().first().unwrap(),
        ] {
            let entities = parse_mmcif_entity_polymer_values(block, Some(&mut grouping)).unwrap();
            assert_eq!(
                c11_subchain_values(&entities),
                [vec!["x", "z"], vec!["y", "x"]]
            );
        }

        let entities = parse_mmcif_entity_polymer_values(&absent.blocks()[0], None).unwrap();
        assert_eq!(
            c11_subchain_values(&entities),
            [Vec::<&str>::new(), Vec::new()]
        );
    }

    #[test]
    fn bio_mmcif_c10_dbrefs_keep_source_order_raw_duplicate_key_and_first_reference() {
        let doc = document(
            concat!(
                "data_c10_dbrefs_ordered\n",
                "loop_\n",
                "_entity.id\n",
                "_entity.type\n",
                "e1 polymer\n",
                "e2 polymer\n",
                "loop_\n",
                "_struct_ref.id\n",
                "_struct_ref.entity_id\n",
                "_struct_ref.db_name\n",
                "_struct_ref.db_code\n",
                "_struct_ref.pdbx_db_accession\n",
                "_struct_ref.pdbx_db_isoform\n",
                "r1 e1 UniProt P12345 ACC1 isoA\n",
                "r1 e2 Wrong WRONG WACC WISO\n",
                "r2 e2 PDB 2XYZ . .\n",
                "r3 absent EMDB EMD-1 ACC3 .\n",
                "loop_\n",
                "_struct_ref_seq.ref_id\n",
                "_struct_ref_seq.seq_align_beg\n",
                "_struct_ref_seq.seq_align_end\n",
                "_struct_ref_seq.db_align_beg\n",
                "_struct_ref_seq.db_align_end\n",
                "_struct_ref_seq.pdbx_auth_seq_align_beg\n",
                "_struct_ref_seq.pdbx_seq_align_beg_ins_code\n",
                "_struct_ref_seq.pdbx_auth_seq_align_end\n",
                "_struct_ref_seq.pdbx_seq_align_end_ins_code\n",
                "r1 2 10 4 12 5 A 14 B\n",
                "r1 2 10 4 12 9 C 15 D\n",
                "\"r1\" 2 10 4 12 6 F 16 G\n",
                "r1 3 11 8 16 . ? ? ?\n",
                "r2 1 5 2 6 3 . 7 .\n",
                "r3 1 3 1 3 1 ? 2 ?\n",
            ),
            "c10-dbrefs-source-order.cif",
        );

        let entities = parse_mmcif_entity_polymer_values(&doc.blocks()[0], None).unwrap();
        assert_eq!(entities.len(), 2);
        assert_eq!(entities[0].source().source_entity_id(), "e1");
        assert_eq!(
            entities[0].dbrefs(),
            [
                BioEntityDbRef {
                    db_name: "UniProt".to_owned(),
                    accession_code: "ACC1".to_owned(),
                    id_code: "P12345".to_owned(),
                    isoform: "isoA".to_owned(),
                    seq_begin: PdbSeqId::new(5, Some(b'A')),
                    seq_end: PdbSeqId::new(14, Some(b'B')),
                    db_begin: PdbSeqId::new(4, None),
                    db_end: PdbSeqId::new(12, None),
                    label_seq_begin: Some(2),
                    label_seq_end: Some(10),
                },
                BioEntityDbRef {
                    db_name: "UniProt".to_owned(),
                    accession_code: "ACC1".to_owned(),
                    id_code: "P12345".to_owned(),
                    isoform: "isoA".to_owned(),
                    seq_begin: PdbSeqId::new(6, Some(b'F')),
                    seq_end: PdbSeqId::new(16, Some(b'G')),
                    db_begin: PdbSeqId::new(4, None),
                    db_end: PdbSeqId::new(12, None),
                    label_seq_begin: Some(2),
                    label_seq_end: Some(10),
                },
                BioEntityDbRef {
                    db_name: "UniProt".to_owned(),
                    accession_code: "ACC1".to_owned(),
                    id_code: "P12345".to_owned(),
                    isoform: "isoA".to_owned(),
                    seq_begin: PdbSeqId::new(i32::MIN, None),
                    seq_end: PdbSeqId::new(i32::MIN, None),
                    db_begin: PdbSeqId::new(8, None),
                    db_end: PdbSeqId::new(16, None),
                    label_seq_begin: Some(3),
                    label_seq_end: Some(11),
                },
            ]
        );
        assert_eq!(entities[1].source().source_entity_id(), "e2");
        assert_eq!(
            entities[1].dbrefs(),
            [BioEntityDbRef {
                db_name: "PDB".to_owned(),
                accession_code: String::new(),
                id_code: "2XYZ".to_owned(),
                isoform: String::new(),
                seq_begin: PdbSeqId::new(3, None),
                seq_end: PdbSeqId::new(7, None),
                db_begin: PdbSeqId::new(2, None),
                db_end: PdbSeqId::new(6, None),
                label_seq_begin: Some(1),
                label_seq_end: Some(5),
            }]
        );
    }

    #[test]
    fn bio_mmcif_c10_dbrefs_keep_absent_optional_columns_and_null_sentinels() {
        let doc = document(
            concat!(
                "data_c10_dbrefs_optional_absent\n",
                "_entity.id e1\n",
                "_entity.type polymer\n",
                "_struct_ref.id r1\n",
                "_struct_ref.entity_id e1\n",
                "_struct_ref.db_name PDB\n",
                "_struct_ref.db_code 1ABC\n",
                "loop_\n",
                "_struct_ref_seq.ref_id\n",
                "_struct_ref_seq.seq_align_beg\n",
                "_struct_ref_seq.seq_align_end\n",
                "_struct_ref_seq.db_align_beg\n",
                "_struct_ref_seq.db_align_end\n",
                "r1 . ? 12 13\n",
            ),
            "c10-dbrefs-optional-absent.cif",
        );

        let entities = parse_mmcif_entity_polymer_values(&doc.blocks()[0], None).unwrap();
        assert_eq!(entities.len(), 1);
        assert_eq!(
            entities[0].dbrefs(),
            [BioEntityDbRef {
                db_name: "PDB".to_owned(),
                accession_code: String::new(),
                id_code: "1ABC".to_owned(),
                isoform: String::new(),
                seq_begin: PdbSeqId::new(i32::MIN, None),
                seq_end: PdbSeqId::new(i32::MIN, None),
                db_begin: PdbSeqId::new(12, None),
                db_end: PdbSeqId::new(13, None),
                label_seq_begin: None,
                label_seq_end: None,
            }]
        );
    }

    #[test]
    fn bio_mmcif_c10_missing_reference_row_returns_structured_cif_error() {
        let doc = document(
            concat!(
                "data_c10_missing_reference\n",
                "_entity.id e1\n",
                "_entity.type polymer\n",
                "_struct_ref.id present\n",
                "_struct_ref.entity_id e1\n",
                "_struct_ref.db_name PDB\n",
                "_struct_ref.db_code 1ABC\n",
                "loop_\n",
                "_struct_ref_seq.ref_id\n",
                "_struct_ref_seq.seq_align_beg\n",
                "_struct_ref_seq.seq_align_end\n",
                "_struct_ref_seq.db_align_beg\n",
                "_struct_ref_seq.db_align_end\n",
                "missing 1 2 3 4\n",
            ),
            "c10-missing-reference.cif",
        );

        let error = parse_mmcif_entity_polymer_values(&doc.blocks()[0], None).unwrap_err();
        let MmcifEntityInfoError::Cif(error) = error else {
            panic!("missing reference must retain the CIF table lookup error")
        };
        assert_eq!(error.kind(), CifReadErrorKind::InvalidValue);
        assert_eq!(error.message(), "Not found in _struct_ref.id: missing");
    }

    #[test]
    fn bio_mmcif_c10_invalid_db_number_and_insertion_conflict_keep_typed_errors() {
        let invalid_number = document(
            concat!(
                "data_c10_invalid_db_number\n",
                "_entity.id e1\n",
                "_entity.type polymer\n",
                "_struct_ref.id r1\n",
                "_struct_ref.entity_id e1\n",
                "_struct_ref.db_name PDB\n",
                "_struct_ref.db_code 1ABC\n",
                "loop_\n",
                "_struct_ref_seq.ref_id\n",
                "_struct_ref_seq.seq_align_beg\n",
                "_struct_ref_seq.seq_align_end\n",
                "_struct_ref_seq.db_align_beg\n",
                "_struct_ref_seq.db_align_end\n",
                "r1 1 2 not-an-integer 4\n",
            ),
            "c10-invalid-db-number.cif",
        );
        let error =
            parse_mmcif_entity_polymer_values(&invalid_number.blocks()[0], None).unwrap_err();
        assert_eq!(error.kind(), CifReadErrorKind::InvalidValue);
        assert!(error.to_string().contains("not an integer: not-an-integer"));

        let insertion_conflict = document(
            concat!(
                "data_c10_insertion_conflict\n",
                "_entity.id e1\n",
                "_entity.type polymer\n",
                "_struct_ref.id r1\n",
                "_struct_ref.entity_id e1\n",
                "_struct_ref.db_name PDB\n",
                "_struct_ref.db_code 1ABC\n",
                "loop_\n",
                "_struct_ref_seq.ref_id\n",
                "_struct_ref_seq.seq_align_beg\n",
                "_struct_ref_seq.seq_align_end\n",
                "_struct_ref_seq.db_align_beg\n",
                "_struct_ref_seq.db_align_end\n",
                "_struct_ref_seq.pdbx_auth_seq_align_beg\n",
                "_struct_ref_seq.pdbx_seq_align_beg_ins_code\n",
                "r1 1 2 3 4 5A B\n",
            ),
            "c10-insertion-conflict.cif",
        );
        let error =
            parse_mmcif_entity_polymer_values(&insertion_conflict.blocks()[0], None).unwrap_err();
        assert!(matches!(
            error,
            MmcifEntityInfoError::SequenceId(
                MmcifSequenceIdError::InconsistentInsertionCode { sequence_id }
            ) if sequence_id == "5A"
        ));
    }

    #[test]
    fn bio_mmcif_c12_sifts_preserves_row_order_accession_indices_and_model_mapping() {
        let doc = c12_sifts_document(
            &[
                "missing absent 2 1 y Z 1 SKIP-ORDINAL",
                "missing absent 1 1 n Z 1 SKIP-OBSERVED",
                "missing absent 1 1 y Z . .",
                "e1 x 1 1 y A 65535 ACC-A",
                "e1 x 1 2 y B 0 ACC-B",
                "e1 x 1 3 y C 12 ACC-A",
                "e1 x 1 1 y D 9 ACC-A",
            ],
            "c12-sifts-source-order.cif",
        );
        let mut grouping = c12_atom_site_grouping();
        let entities =
            parse_mmcif_entity_polymer_values(&doc.blocks()[0], Some(&mut grouping)).unwrap();

        assert_eq!(entities.len(), 1);
        assert_eq!(entities[0].sifts_unp_accessions(), ["ACC-A", "ACC-B"]);
        assert_eq!(
            grouping
                .models
                .iter()
                .map(|model| model.source_model_number)
                .collect::<Vec<_>>(),
            [2, 1]
        );
        for model in &grouping.models {
            let residues = &model.chains[0].residues;
            assert_eq!(residues.len(), 3);
            assert_eq!(
                residues
                    .iter()
                    .map(|residue| residue.sifts_unp)
                    .collect::<Vec<_>>(),
                [
                    BioSiftsUnpResidue::new(Some(b'D'), 0, 9),
                    BioSiftsUnpResidue::new(Some(b'B'), 1, 0),
                    BioSiftsUnpResidue::new(Some(b'C'), 0, 12),
                ]
            );
        }
    }

    #[test]
    fn bio_mmcif_c12_sifts_absent_empty_and_incomplete_tables_are_noops() {
        let absent = document(
            concat!(
                "data_c12_sifts_absent\n",
                "_entity.id e1\n",
                "_entity.type polymer\n",
                "loop_\n",
                "_struct_asym.id\n",
                "_struct_asym.entity_id\n",
                "x e1\n",
            ),
            "c12-sifts-absent.cif",
        );
        let empty = c12_sifts_document(&[], "c12-sifts-empty.cif");
        let incomplete = document(
            concat!(
                "data_c12_sifts_incomplete\n",
                "_entity.id e1\n",
                "_entity.type polymer\n",
                "loop_\n",
                "_pdbx_sifts_xref_db.entity_id\n",
                "e1\n",
            ),
            "c12-sifts-incomplete.cif",
        );

        for doc in [&absent, &empty, &incomplete] {
            let mut grouping = c12_atom_site_grouping();
            let entities =
                parse_mmcif_entity_polymer_values(&doc.blocks()[0], Some(&mut grouping)).unwrap();
            assert!(entities[0].sifts_unp_accessions().is_empty());
            assert!(grouping.models.iter().all(|model| {
                model.chains.iter().all(|chain| {
                    chain
                        .residues
                        .iter()
                        .all(|residue| residue.sifts_unp == BioSiftsUnpResidue::default())
                })
            }));
        }
    }

    #[test]
    fn bio_mmcif_c12_sifts_reports_typed_missing_rows_and_numeric_errors() {
        let missing_entity = c12_parse_error("missing x 1 1 y A 1 ACC");
        assert!(matches!(
            &missing_entity,
            MmcifEntityInfoError::Sifts(MmcifSiftsUnpError::EntityNotFound { entity_id })
                if entity_id == "missing"
        ));
        assert_eq!(
            missing_entity.to_string(),
            "_pdbx_sifts_xref_db: entity_id not found: missing"
        );

        let missing_subchain = c12_parse_error("e1 absent 1 1 y A 1 ACC");
        assert!(matches!(
            &missing_subchain,
            MmcifEntityInfoError::Sifts(MmcifSiftsUnpError::SubchainNotFound { asym_id })
                if asym_id == "absent"
        ));
        assert_eq!(
            missing_subchain.to_string(),
            "_pdbx_sifts_xref_db: asym_id not found: absent"
        );

        let missing_sequence = c12_parse_error("e1 x 1 9 y A 1 ACC");
        assert!(matches!(
            &missing_sequence,
            MmcifEntityInfoError::Sifts(MmcifSiftsUnpError::SequenceNotFound { seq_id })
                if seq_id == "9"
        ));
        assert_eq!(
            missing_sequence.to_string(),
            "_pdbx_sifts_xref_db: seq_id not found: 9"
        );

        for number in ["65536", "-1"] {
            let row = format!("e1 x 1 1 y A {number} ACC");
            let error = c12_parse_error(&row);
            assert!(matches!(
                &error,
                MmcifEntityInfoError::Sifts(MmcifSiftsUnpError::NumberOutsideU16 { value })
                    if value == number
            ));
            assert_eq!(
                error.to_string(),
                format!("_pdbx_sifts_xref_db.unp_num: {number}")
            );
        }

        let invalid_number = c12_parse_error("e1 x 1 1 y A invalid ACC");
        assert!(matches!(
            &invalid_number,
            MmcifEntityInfoError::Sifts(MmcifSiftsUnpError::Cif(error))
                if error.kind() == CifReadErrorKind::InvalidValue
                    && error.message() == "not an integer: invalid"
        ));
        let invalid_residue = c12_parse_error("e1 x 1 1 y AB 1 ACC");
        assert!(matches!(
            &invalid_residue,
            MmcifEntityInfoError::Sifts(MmcifSiftsUnpError::Cif(error))
                if error.kind() == CifReadErrorKind::InvalidValue
                    && error.message() == "Not a single character: AB"
        ));
    }

    #[test]
    fn bio_mmcif_c22_helices_filter_by_first_byte_and_preserve_source_order() {
        let doc = document(
            concat!(
                "data_c22_helix_order\n",
                "loop_\n",
                "_struct_conf.conf_type_id\n",
                "_struct_conf.beg_auth_asym_id\n",
                "_struct_conf.beg_label_comp_id\n",
                "_struct_conf.beg_auth_seq_id\n",
                "_struct_conf.pdbx_beg_PDB_ins_code\n",
                "_struct_conf.end_auth_asym_id\n",
                "_struct_conf.end_label_comp_id\n",
                "_struct_conf.end_auth_seq_id\n",
                "_struct_conf.pdbx_end_PDB_ins_code\n",
                "_struct_conf.pdbx_PDB_helix_class\n",
                "_struct_conf.pdbx_PDB_helix_length\n",
                "HELX_P A ALA 15A . B GLY 24 B 1 7\n",
                "h B SER -5 ? A THR 8 ? 99 .\n",
                "TURN A XXX invalid . B YYY invalid . invalid invalid\n",
            ),
            "c22-helix-order.cif",
        );

        let helices = parse_mmcif_helices(&doc.blocks()[0]).unwrap();
        assert_eq!(helices.len(), 2);
        assert_eq!(helices[0].start.chain_name().as_str(), "A");
        assert_eq!(helices[0].start.residue().name().as_str(), "ALA");
        assert_eq!(helices[0].start.residue().sequence_number(), Some(15));
        assert_eq!(helices[0].start.residue().insertion_code(), Some(b'A'));
        assert_eq!(helices[0].start.logical_atom_name(), "");
        assert_eq!(helices[0].start.altloc(), 0);
        assert_eq!(helices[0].end.chain_name().as_str(), "B");
        assert_eq!(helices[0].end.residue().name().as_str(), "GLY");
        assert_eq!(helices[0].end.residue().sequence_number(), Some(24));
        assert_eq!(helices[0].end.residue().insertion_code(), Some(b'B'));
        assert_eq!(helices[0].pdb_helix_class, BioHelixClass::RAlpha);
        assert_eq!(helices[0].length, 7);

        assert_eq!(helices[1].start.chain_name().as_str(), "B");
        assert_eq!(helices[1].start.residue().sequence_number(), Some(-5));
        assert_eq!(helices[1].start.residue().insertion_code(), None);
        assert_eq!(helices[1].end.chain_name().as_str(), "A");
        assert_eq!(helices[1].end.residue().sequence_number(), Some(8));
        assert_eq!(helices[1].end.residue().insertion_code(), None);
        assert_eq!(helices[1].pdb_helix_class, BioHelixClass::UnknownHelix);
        assert_eq!(helices[1].length, -1);
    }

    #[test]
    fn bio_mmcif_c22_helices_distinguish_missing_null_and_empty_optional_fields() {
        let missing_columns = document(
            concat!(
                "data_c22_helix_missing_columns\n",
                "loop_\n",
                "_struct_conf.conf_type_id\n",
                "_struct_conf.beg_auth_asym_id\n",
                "_struct_conf.beg_label_comp_id\n",
                "_struct_conf.beg_auth_seq_id\n",
                "_struct_conf.end_auth_asym_id\n",
                "_struct_conf.end_label_comp_id\n",
                "_struct_conf.end_auth_seq_id\n",
                "hElX_P A ALA 12A B GLY 13\n",
            ),
            "c22-helix-missing-columns.cif",
        );
        let absent = parse_mmcif_helices(&missing_columns.blocks()[0]).unwrap();
        assert_eq!(absent.len(), 1);
        assert_eq!(absent[0].start.residue().sequence_number(), Some(12));
        assert_eq!(absent[0].start.residue().insertion_code(), Some(b'A'));
        assert_eq!(absent[0].end.residue().sequence_number(), Some(13));
        assert_eq!(absent[0].end.residue().insertion_code(), None);
        assert_eq!(absent[0].pdb_helix_class, BioHelixClass::UnknownHelix);
        assert_eq!(absent[0].length, -1);

        let present_columns = document(
            concat!(
                "data_c22_helix_present_columns\n",
                "loop_\n",
                "_struct_conf.conf_type_id\n",
                "_struct_conf.beg_auth_asym_id\n",
                "_struct_conf.beg_label_comp_id\n",
                "_struct_conf.beg_auth_seq_id\n",
                "_struct_conf.pdbx_beg_PDB_ins_code\n",
                "_struct_conf.end_auth_asym_id\n",
                "_struct_conf.end_label_comp_id\n",
                "_struct_conf.end_auth_seq_id\n",
                "_struct_conf.pdbx_end_PDB_ins_code\n",
                "_struct_conf.pdbx_PDB_helix_class\n",
                "_struct_conf.pdbx_PDB_helix_length\n",
                "HELX_P A ALA 12 . B GLY 13 ? 1 .\n",
                "HELX_P A ALA 14 '' B GLY 15B B . 8\n",
            ),
            "c22-helix-present-columns.cif",
        );
        let present = parse_mmcif_helices(&present_columns.blocks()[0]).unwrap();
        assert_eq!(present.len(), 2);
        assert_eq!(present[0].start.residue().insertion_code(), None);
        assert_eq!(present[0].end.residue().insertion_code(), None);
        assert_eq!(present[0].pdb_helix_class, BioHelixClass::RAlpha);
        assert_eq!(present[0].length, -1);
        assert_eq!(present[1].start.residue().sequence_number(), Some(14));
        assert_eq!(present[1].start.residue().insertion_code(), Some(0));
        assert_eq!(present[1].end.residue().sequence_number(), Some(15));
        assert_eq!(present[1].end.residue().insertion_code(), Some(b'B'));
        assert_eq!(present[1].pdb_helix_class, BioHelixClass::UnknownHelix);
        assert_eq!(present[1].length, 8);
    }

    #[test]
    fn bio_mmcif_c22_helices_propagate_typed_sequence_width_and_integer_errors() {
        use super::MmcifSequenceIdError;

        let header = concat!(
            "data_c22_helix_errors\n",
            "loop_\n",
            "_struct_conf.conf_type_id\n",
            "_struct_conf.beg_auth_asym_id\n",
            "_struct_conf.beg_label_comp_id\n",
            "_struct_conf.beg_auth_seq_id\n",
            "_struct_conf.pdbx_beg_PDB_ins_code\n",
            "_struct_conf.end_auth_asym_id\n",
            "_struct_conf.end_label_comp_id\n",
            "_struct_conf.end_auth_seq_id\n",
            "_struct_conf.pdbx_end_PDB_ins_code\n",
            "_struct_conf.pdbx_PDB_helix_class\n",
            "_struct_conf.pdbx_PDB_helix_length\n",
        );
        let parse_row = |row: &str| {
            let text = format!("{header}{row}\n");
            let doc = document(&text, "c22-helix-error.cif");
            parse_mmcif_helices(&doc.blocks()[0])
        };

        let insertion_conflict = parse_row("HELX_P A ALA 15A B B GLY 16 . 1 7").unwrap_err();
        assert!(matches!(
            insertion_conflict,
            MmcifHelixError::SequenceId {
                endpoint: "start",
                error: MmcifSequenceIdError::InconsistentInsertionCode { .. }
            }
        ));

        let bad_class = parse_row("HELX_P A ALA 1 . B GLY 2 . nope 5").unwrap_err();
        assert!(matches!(
            bad_class,
            MmcifHelixError::HelixClass(error)
                if error.kind() == CifReadErrorKind::InvalidValue
        ));
        let bad_length = parse_row("HELX_P A ALA 1 . B GLY 2 . . nope").unwrap_err();
        assert!(matches!(
            bad_length,
            MmcifHelixError::HelixLength(error)
                if error.kind() == CifReadErrorKind::InvalidValue
        ));

        assert!(matches!(
            parse_row("HELX_P ABCDE ALA 1 . B GLY 2 . 1 7").unwrap_err(),
            MmcifHelixError::ChainNameNotRepresentable {
                endpoint: "start",
                ..
            }
        ));
        assert!(matches!(
            parse_row("HELX_P A ABCDE 1 . B GLY 2 . 1 7").unwrap_err(),
            MmcifHelixError::ResidueNameNotRepresentable {
                endpoint: "start",
                ..
            }
        ));
    }

    #[test]
    fn bio_mmcif_c23_sheets_preserve_source_order_duplicates_sense_and_foreign_ids() {
        let doc = document(
            concat!(
                "data_c23_sheet_order\n",
                "loop_\n",
                "_struct_sheet.id\n",
                "S1\n",
                "S2\n",
                "S1\n",
                "loop_\n",
                "_struct_sheet_range.sheet_id\n",
                "_struct_sheet_range.id\n",
                "_struct_sheet_range.beg_auth_asym_id\n",
                "_struct_sheet_range.beg_label_comp_id\n",
                "_struct_sheet_range.beg_auth_seq_id\n",
                "_struct_sheet_range.pdbx_beg_PDB_ins_code\n",
                "_struct_sheet_range.end_auth_asym_id\n",
                "_struct_sheet_range.end_label_comp_id\n",
                "_struct_sheet_range.end_auth_seq_id\n",
                "_struct_sheet_range.pdbx_end_PDB_ins_code\n",
                "S1 A A ALA 10A . B GLY 20 .\n",
                "S1 B C SER 30 . D THR 40 .\n",
                "S1 A E ASN 50 . F VAL 60 .\n",
                "S3 D G CYS 70 . H LEU 80 .\n",
                "loop_\n",
                "_struct_sheet_order.sheet_id\n",
                "_struct_sheet_order.range_id_2\n",
                "_struct_sheet_order.sense\n",
                "S1 A P\n",
                "S1 A x\n",
                "S1 B a\n",
                "S1 D ?\n",
                "ghost missing P\n",
                "S3 D p\n",
            ),
            "c23-sheet-order.cif",
        );

        let sheets = parse_mmcif_sheets(&doc.blocks()[0]).unwrap();
        assert_eq!(
            sheets
                .iter()
                .map(|sheet| sheet.name.as_str())
                .collect::<Vec<_>>(),
            ["S1", "S2", "S1", "S3"]
        );
        assert_eq!(sheets[0].strands.len(), 3);
        assert_eq!(sheets[0].strands[0].name, "A");
        assert_eq!(sheets[0].strands[1].name, "B");
        assert_eq!(sheets[0].strands[2].name, "A");
        assert_eq!(sheets[0].strands[0].sense, 1);
        assert_eq!(sheets[0].strands[1].sense, -1);
        assert_eq!(sheets[0].strands[2].sense, 0);
        assert_eq!(sheets[0].strands[0].start.chain_name().as_str(), "A");
        assert_eq!(sheets[0].strands[0].start.residue().name().as_str(), "ALA");
        assert_eq!(
            sheets[0].strands[0].start.residue().sequence_number(),
            Some(10)
        );
        assert_eq!(
            sheets[0].strands[0].start.residue().insertion_code(),
            Some(b'A')
        );
        assert_eq!(sheets[0].strands[0].end.chain_name().as_str(), "B");
        assert_eq!(sheets[0].strands[0].end.residue().name().as_str(), "GLY");
        assert_eq!(
            sheets[0].strands[0].end.residue().sequence_number(),
            Some(20)
        );
        assert_eq!(sheets[0].strands[0].end.residue().insertion_code(), None);
        assert!(sheets[1].strands.is_empty());
        assert!(sheets[2].strands.is_empty());
        assert_eq!(sheets[3].strands.len(), 1);
        assert_eq!(sheets[3].strands[0].name, "D");
        assert_eq!(sheets[3].strands[0].sense, 1);
    }

    #[test]
    fn bio_mmcif_c23_sheet_hbonds_resolve_first_ids_and_retain_full_addresses() {
        let doc = document(
            concat!(
                "data_c23_sheet_hbond\n",
                "loop_\n",
                "_struct_sheet.id\n",
                "S\n",
                "S\n",
                "loop_\n",
                "_struct_sheet_range.sheet_id\n",
                "_struct_sheet_range.id\n",
                "_struct_sheet_range.beg_auth_asym_id\n",
                "_struct_sheet_range.beg_label_comp_id\n",
                "_struct_sheet_range.beg_auth_seq_id\n",
                "_struct_sheet_range.end_auth_asym_id\n",
                "_struct_sheet_range.end_label_comp_id\n",
                "_struct_sheet_range.end_auth_seq_id\n",
                "S R R ALA 1 X GLY 2\n",
                "loop_\n",
                "_pdbx_struct_sheet_hbond.sheet_id\n",
                "_pdbx_struct_sheet_hbond.range_id_2\n",
                "_pdbx_struct_sheet_hbond.range_1_auth_asym_id\n",
                "_pdbx_struct_sheet_hbond.range_1_label_comp_id\n",
                "_pdbx_struct_sheet_hbond.range_1_auth_seq_id\n",
                "_pdbx_struct_sheet_hbond.range_1_PDB_ins_code\n",
                "_pdbx_struct_sheet_hbond.range_1_label_atom_id\n",
                "_pdbx_struct_sheet_hbond.range_2_auth_asym_id\n",
                "_pdbx_struct_sheet_hbond.range_2_label_comp_id\n",
                "_pdbx_struct_sheet_hbond.range_2_auth_seq_id\n",
                "_pdbx_struct_sheet_hbond.range_2_PDB_ins_code\n",
                "_pdbx_struct_sheet_hbond.range_2_label_atom_id\n",
                "ghost R ABCDE WRONG invalid AB BADATOM ABCDE WRONG invalid AB BADATOM\n",
                "S R A ASN 11A . N B GLY 12 ? O\n",
                "S R C SER 21 . LONGATOM D THR 22 . O2\n",
            ),
            "c23-sheet-hbond.cif",
        );

        let sheets = parse_mmcif_sheets(&doc.blocks()[0]).unwrap();
        assert_eq!(sheets.len(), 2);
        assert_eq!(sheets[0].strands.len(), 1);
        assert!(sheets[1].strands.is_empty());
        let strand = &sheets[0].strands[0];
        assert_eq!(strand.hbond_atom1.chain_name().as_str(), "C");
        assert_eq!(strand.hbond_atom1.residue().name().as_str(), "SER");
        assert_eq!(strand.hbond_atom1.residue().sequence_number(), Some(21));
        assert_eq!(strand.hbond_atom1.residue().insertion_code(), None);
        assert_eq!(strand.hbond_atom1.logical_atom_name(), "LONGATOM");
        assert_eq!(strand.hbond_atom1.altloc(), 0);
        assert_eq!(strand.hbond_atom2.chain_name().as_str(), "D");
        assert_eq!(strand.hbond_atom2.residue().name().as_str(), "THR");
        assert_eq!(strand.hbond_atom2.residue().sequence_number(), Some(22));
        assert_eq!(strand.hbond_atom2.residue().insertion_code(), None);
        assert_eq!(strand.hbond_atom2.logical_atom_name(), "O2");
        assert_eq!(strand.hbond_atom2.altloc(), 0);
    }

    #[test]
    fn bio_mmcif_c23_sheets_match_pair_null_optional_and_missing_category_paths() {
        let pair_values = document(
            concat!(
                "data_c23_sheet_pairs\n",
                "_struct_sheet.id S\n",
                "_struct_sheet_range.sheet_id S\n",
                "_struct_sheet_range.id R\n",
                "_struct_sheet_range.beg_auth_asym_id A\n",
                "_struct_sheet_range.beg_label_comp_id ALA\n",
                "_struct_sheet_range.beg_auth_seq_id 15A\n",
                "_struct_sheet_range.end_auth_asym_id B\n",
                "_struct_sheet_range.end_label_comp_id .\n",
                "_struct_sheet_range.end_auth_seq_id ?\n",
                "_struct_sheet_order.sheet_id S\n",
                "_struct_sheet_order.range_id_2 R\n",
                "_struct_sheet_order.sense p\n",
                "_pdbx_struct_sheet_hbond.sheet_id S\n",
                "_pdbx_struct_sheet_hbond.range_id_2 R\n",
                "_pdbx_struct_sheet_hbond.range_1_auth_asym_id X\n",
                "_pdbx_struct_sheet_hbond.range_1_label_comp_id ASN\n",
                "_pdbx_struct_sheet_hbond.range_1_auth_seq_id 21\n",
                "_pdbx_struct_sheet_hbond.range_1_label_atom_id N\n",
                "_pdbx_struct_sheet_hbond.range_2_auth_asym_id Y\n",
                "_pdbx_struct_sheet_hbond.range_2_label_comp_id .\n",
                "_pdbx_struct_sheet_hbond.range_2_auth_seq_id ?\n",
                "_pdbx_struct_sheet_hbond.range_2_label_atom_id ?\n",
            ),
            "c23-sheet-pairs.cif",
        );
        let sheets = parse_mmcif_sheets(&pair_values.blocks()[0]).unwrap();
        assert_eq!(sheets.len(), 1);
        assert_eq!(sheets[0].strands.len(), 1);
        let strand = &sheets[0].strands[0];
        assert_eq!(strand.sense, 1);
        assert_eq!(strand.start.residue().sequence_number(), Some(15));
        assert_eq!(strand.start.residue().insertion_code(), Some(b'A'));
        assert_eq!(strand.end.residue().name().as_str(), "");
        assert_eq!(strand.end.residue().sequence_number(), None);
        assert_eq!(strand.hbond_atom1.residue().sequence_number(), Some(21));
        assert_eq!(strand.hbond_atom1.logical_atom_name(), "N");
        assert_eq!(strand.hbond_atom2.chain_name().as_str(), "Y");
        assert_eq!(strand.hbond_atom2.residue().name().as_str(), "");
        assert_eq!(strand.hbond_atom2.residue().sequence_number(), None);
        assert_eq!(strand.hbond_atom2.logical_atom_name(), "");

        let empty = document("data_c23_sheet_empty\n", "c23-sheet-empty.cif");
        assert!(parse_mmcif_sheets(&empty.blocks()[0]).unwrap().is_empty());

        let missing_required = document(
            concat!(
                "data_c23_sheet_missing_required\n",
                "_struct_sheet.id S\n",
                "loop_\n",
                "_struct_sheet_range.sheet_id\n",
                "_struct_sheet_range.id\n",
                "S R\n",
            ),
            "c23-sheet-missing-required.cif",
        );
        let sheets = parse_mmcif_sheets(&missing_required.blocks()[0]).unwrap();
        assert_eq!(sheets.len(), 1);
        assert!(sheets[0].strands.is_empty());
    }

    #[test]
    fn bio_mmcif_c23_sheet_address_failures_remain_typed() {
        let header = concat!(
            "data_c23_sheet_errors\n",
            "loop_\n",
            "_struct_sheet_range.sheet_id\n",
            "_struct_sheet_range.id\n",
            "_struct_sheet_range.beg_auth_asym_id\n",
            "_struct_sheet_range.beg_label_comp_id\n",
            "_struct_sheet_range.beg_auth_seq_id\n",
            "_struct_sheet_range.pdbx_beg_PDB_ins_code\n",
            "_struct_sheet_range.end_auth_asym_id\n",
            "_struct_sheet_range.end_label_comp_id\n",
            "_struct_sheet_range.end_auth_seq_id\n",
            "_struct_sheet_range.pdbx_end_PDB_ins_code\n",
        );
        let parse_range = |row: &str| {
            let text = format!("{header}{row}\n");
            let doc = document(&text, "c23-sheet-address-error.cif");
            parse_mmcif_sheets(&doc.blocks()[0])
        };

        assert!(matches!(
            parse_range("S R ABCDE ALA 1 . B GLY 2 .").unwrap_err(),
            MmcifSheetError::ChainNameNotRepresentable {
                address: "start",
                value
            } if value == "ABCDE"
        ));
        assert!(matches!(
            parse_range("S R A ALA 1 . B ABCDE 2 .").unwrap_err(),
            MmcifSheetError::ResidueNameNotRepresentable {
                address: "end",
                value
            } if value == "ABCDE"
        ));
        assert!(matches!(
            parse_range("S R A ALA invalid . B GLY 2 .").unwrap_err(),
            MmcifSheetError::SequenceId {
                address: "start",
                error: MmcifSequenceIdError::InvalidSequenceNumber(error)
            } if error.kind() == CifReadErrorKind::InvalidValue
        ));
        assert!(matches!(
            parse_range("S R A ALA 1 AB B GLY 2 .").unwrap_err(),
            MmcifSheetError::SequenceId {
                address: "start",
                error: MmcifSequenceIdError::InvalidInsertionCode(error)
            } if error.kind() == CifReadErrorKind::InvalidValue
        ));

        let hbond = document(
            concat!(
                "data_c23_sheet_hbond_error\n",
                "loop_\n",
                "_struct_sheet_range.sheet_id\n",
                "_struct_sheet_range.id\n",
                "_struct_sheet_range.beg_auth_asym_id\n",
                "_struct_sheet_range.beg_label_comp_id\n",
                "_struct_sheet_range.beg_auth_seq_id\n",
                "_struct_sheet_range.end_auth_asym_id\n",
                "_struct_sheet_range.end_label_comp_id\n",
                "_struct_sheet_range.end_auth_seq_id\n",
                "S R A ALA 1 B GLY 2\n",
                "loop_\n",
                "_pdbx_struct_sheet_hbond.sheet_id\n",
                "_pdbx_struct_sheet_hbond.range_id_2\n",
                "_pdbx_struct_sheet_hbond.range_1_auth_asym_id\n",
                "_pdbx_struct_sheet_hbond.range_1_label_comp_id\n",
                "_pdbx_struct_sheet_hbond.range_1_auth_seq_id\n",
                "_pdbx_struct_sheet_hbond.range_1_label_atom_id\n",
                "_pdbx_struct_sheet_hbond.range_2_auth_asym_id\n",
                "_pdbx_struct_sheet_hbond.range_2_label_comp_id\n",
                "_pdbx_struct_sheet_hbond.range_2_auth_seq_id\n",
                "_pdbx_struct_sheet_hbond.range_2_label_atom_id\n",
                "S R ABCDE ALA 3 N B GLY 4 O\n",
            ),
            "c23-sheet-hbond-error.cif",
        );
        assert!(matches!(
            parse_mmcif_sheets(&hbond.blocks()[0]).unwrap_err(),
            MmcifSheetError::ChainNameNotRepresentable {
                address: "hbond_atom1",
                value
            } if value == "ABCDE"
        ));
    }

    #[test]
    fn bio_mmcif_c24_label_address_follows_chain_and_contiguous_span_order() {
        let rows = [
            c11_atom_site_row("1", "x", "missing_entity", "1", "10A", "A", "1"),
            c11_atom_site_row("2", "x", "missing_entity", "1", "11", "A", "1"),
            c11_atom_site_row("3", "y", "missing_entity", "2", "20", "A", "1"),
            c11_atom_site_row("4", "x", "missing_entity", "3", "30", "A", "1"),
            c11_atom_site_row("5", "x", "missing_entity", "3", "40", "B", "1"),
            c11_atom_site_row("6", "x", "missing_entity", "3", "41", "B", "1"),
        ];
        // No Entity table is present: the pinned helper searches Model chains
        // and residue spans only.
        let doc = c05_scalar_document(&rows, &[]);
        let grouping = group_mmcif_atom_sites(&atom_site_table(&doc.blocks()[0])).unwrap();
        let model = &grouping.models[0];

        assert_eq!(
            resolve_mmcif_label_address(model, "x", "1"),
            Ok(Some((
                PdbChainId::from_ascii(b"A").unwrap(),
                Some(PdbSeqId::new(10, Some(b'A'))),
            )))
        );
        // The second "x" span in chain A must not be searched by get_subchain;
        // after chain A's first span misses label 3, source traversal continues
        // to chain B and takes its first matching author sequence.
        assert_eq!(
            resolve_mmcif_label_address(model, "x", "3"),
            Ok(Some((
                PdbChainId::from_ascii(b"B").unwrap(),
                Some(PdbSeqId::new(40, None)),
            )))
        );
        // A non-null sequence may have Gemmi's accepted C-locale whitespace
        // and sign syntax. A matching subchain with no matching label sequence
        // retains its last source chain assignment but has no sequence update.
        assert_eq!(
            resolve_mmcif_label_address(model, "x", "\t+2\n"),
            Ok(Some((PdbChainId::from_ascii(b"B").unwrap(), None)))
        );
        assert_eq!(resolve_mmcif_label_address(model, "absent", "1"), Ok(None));
    }

    #[test]
    fn bio_mmcif_c24_label_address_maps_null_sequences_to_the_source_sentinel() {
        let rows = [
            c11_atom_site_row("1", "x", "e1", ".", "50", "A", "1"),
            c11_atom_site_row("2", "x", "e1", "?", "60", "A", "1"),
        ];
        let doc = c05_scalar_document(&rows, &[]);
        let grouping = group_mmcif_atom_sites(&atom_site_table(&doc.blocks()[0])).unwrap();
        let model = &grouping.models[0];
        let expected = Ok(Some((
            PdbChainId::from_ascii(b"A").unwrap(),
            Some(PdbSeqId::new(50, None)),
        )));

        // Gemmi maps both null tokens to OptionalNum's INT_MIN value, so the
        // first matching null-labeled residue wins for either raw spelling.
        assert_eq!(resolve_mmcif_label_address(model, "x", "."), expected);
        assert_eq!(resolve_mmcif_label_address(model, "x", "?"), expected);
    }

    #[test]
    fn bio_mmcif_c24_label_address_preserves_typed_conversion_and_width_errors() {
        let rows = [c11_atom_site_row("1", "x", "e1", "1", "10", "A", "1")];
        let doc = c05_scalar_document(&rows, &[]);
        let mut grouping = group_mmcif_atom_sites(&atom_site_table(&doc.blocks()[0])).unwrap();
        let model = &mut grouping.models[0];

        assert!(matches!(
            resolve_mmcif_label_address(model, "x", "1tail"),
            Err(MmcifLabelAddressError::LabelSequenceId(error))
                if error.kind() == CifReadErrorKind::InvalidValue
        ));

        model.chains[0].source_name = "ABCDE".to_owned();
        assert_eq!(
            resolve_mmcif_label_address(model, "x", "1"),
            Err(MmcifLabelAddressError::ChainNameNotRepresentable {
                value: "ABCDE".to_owned(),
            })
        );
    }

    #[test]
    fn bio_mmcif_c25_connection_lowercase_disulf() {
        // Gemmi enumstr.hpp::connection_type_from_string compares exact text;
        // read_connectivity retains source row and partner order.
        let mut disulf = c25_connection_row();
        disulf[0] = "lowercase";
        disulf[1] = "disulf";
        let mut uppercase = c25_connection_row();
        uppercase[0] = "uppercase";
        uppercase[1] = "DISULF";
        let doc = c25_connection_document(&[disulf, uppercase], &[]);
        let connections =
            parse_mmcif_connections(&doc.blocks()[0], &MmcifAtomSiteGrouping::default()).unwrap();
        assert_eq!(connections.len(), 2);
        assert_eq!(connections[0].name, "lowercase");
        assert_eq!(connections[0].kind, BioConnectionKind::Disulf);
        assert_eq!(connections[1].name, "uppercase");
        assert_eq!(connections[1].kind, BioConnectionKind::Unknown);
        for connection in &connections {
            assert_eq!(
                connection.partner1.chain_name(),
                PdbChainId::from_ascii(b"A1").unwrap()
            );
            assert_eq!(
                connection.partner1.residue().name(),
                ResidueName::from_ascii(b"LIG").unwrap()
            );
            assert_eq!(connection.partner1.residue().sequence_number(), Some(15));
            assert_eq!(connection.partner1.residue().insertion_code(), Some(b'A'));
            assert_eq!(connection.partner1.logical_atom_name(), "C1");
            assert_eq!(connection.partner1.altloc(), b'A');
            assert_eq!(
                connection.partner2.chain_name(),
                PdbChainId::from_ascii(b"B2").unwrap()
            );
            assert_eq!(
                connection.partner2.residue().name(),
                ResidueName::from_ascii(b"GLY").unwrap()
            );
            assert_eq!(connection.partner2.residue().sequence_number(), Some(16));
            assert_eq!(connection.partner2.residue().insertion_code(), None);
            assert_eq!(connection.partner2.logical_atom_name(), "N");
            assert_eq!(connection.partner2.altloc(), 0);
        }
    }

    #[test]
    fn bio_mmcif_c25_connection_kind_asu_symmetry_distance_and_defaults_follow_source() {
        let mut covalent = c25_connection_row();
        covalent[0] = "covalent";
        let mut uppercase_kind = c25_connection_row();
        uppercase_kind[0] = "uppercase";
        uppercase_kind[1] = "DISULF";
        uppercase_kind[18] = ".";
        uppercase_kind[19] = ".";
        uppercase_kind[20] = ".";
        uppercase_kind[21] = ".";
        let mut equal_symmetry = c25_connection_row();
        equal_symmetry[0] = "equal";
        equal_symmetry[1] = "hydrog";
        equal_symmetry[18] = "1_555";
        equal_symmetry[19] = "1_555";
        equal_symmetry[20] = "bad";
        equal_symmetry[21] = "link-equal";
        let mut source_duplicate_lookup = c25_connection_row();
        source_duplicate_lookup[0] = "duplicate-lookup";
        source_duplicate_lookup[1] = "metalc";
        source_duplicate_lookup[18] = "A555X";
        source_duplicate_lookup[19] = "2_565";
        source_duplicate_lookup[20] = "0";
        source_duplicate_lookup[21] = "link-duplicate";
        let mut malformed_symmetry = c25_connection_row();
        malformed_symmetry[0] = "malformed-symmetry";
        malformed_symmetry[1] = "unrecognized";
        malformed_symmetry[18] = "1_55";
        malformed_symmetry[19] = "2_565";
        malformed_symmetry[20] = ".";
        malformed_symmetry[21] = ".";

        let rows = [
            covalent,
            uppercase_kind,
            equal_symmetry,
            source_duplicate_lookup,
            malformed_symmetry,
        ];
        let doc = c25_connection_document(&rows, &[]);
        let connections =
            parse_mmcif_connections(&doc.blocks()[0], &MmcifAtomSiteGrouping::default()).unwrap();
        assert_eq!(connections.len(), 5);

        let first = &connections[0];
        assert_eq!(first.name, "covalent");
        assert_eq!(first.link_id, "LINK");
        assert_eq!(first.kind, BioConnectionKind::Covale);
        assert_eq!(first.asu, BioAsu::Different);
        assert_eq!(first.reported_sym, [2, 0, 1, 0]);
        assert_eq!(first.reported_distance, 1.25);
        assert_eq!(
            first.partner1.chain_name(),
            PdbChainId::from_ascii(b"A1").unwrap()
        );
        assert_eq!(
            first.partner1.residue().name(),
            ResidueName::from_ascii(b"LIG").unwrap()
        );
        assert_eq!(first.partner1.residue().sequence_number(), Some(15));
        assert_eq!(first.partner1.residue().insertion_code(), Some(b'A'));
        assert_eq!(first.partner1.logical_atom_name(), "C1");
        assert_eq!(first.partner1.altloc(), b'A');
        assert_eq!(
            first.partner2.chain_name(),
            PdbChainId::from_ascii(b"B2").unwrap()
        );
        assert_eq!(
            first.partner2.residue().name(),
            ResidueName::from_ascii(b"GLY").unwrap()
        );
        assert_eq!(first.partner2.residue().sequence_number(), Some(16));
        assert_eq!(first.partner2.residue().insertion_code(), None);
        assert_eq!(first.partner2.logical_atom_name(), "N");
        assert_eq!(first.partner2.altloc(), 0);

        let uppercase = &connections[1];
        assert_eq!(uppercase.kind, BioConnectionKind::Unknown);
        assert_eq!(uppercase.asu, BioAsu::Any);
        assert_eq!(uppercase.reported_sym, [0; 4]);
        assert_eq!(uppercase.reported_distance, 0.0);
        assert!(uppercase.link_id.is_empty());

        let equal = &connections[2];
        assert_eq!(equal.kind, BioConnectionKind::Hydrog);
        assert_eq!(equal.asu, BioAsu::Same);
        assert_eq!(equal.reported_sym, [0; 4]);
        assert!(equal.reported_distance.is_nan());
        assert_eq!(equal.link_id, "link-equal");

        let duplicate = &connections[3];
        assert_eq!(duplicate.kind, BioConnectionKind::MetalC);
        assert_eq!(duplicate.asu, BioAsu::Different);
        assert_eq!(duplicate.reported_sym, [99, 0, 1, -35]);
        assert_eq!(duplicate.reported_distance, 0.0);

        let malformed = &connections[4];
        assert_eq!(malformed.kind, BioConnectionKind::Unknown);
        assert_eq!(malformed.asu, BioAsu::Different);
        assert_eq!(malformed.reported_sym, [0; 4]);

        let mut absent_optional = c25_connection_row();
        absent_optional[0] = "absent-optional";
        let omitted = [10, 11, 16, 17, 18, 19, 20, 21];
        let omitted_doc = c25_connection_document(&[absent_optional], &omitted);
        let omitted_connection =
            parse_mmcif_connections(&omitted_doc.blocks()[0], &MmcifAtomSiteGrouping::default())
                .unwrap()
                .remove(0);
        assert_eq!(omitted_connection.asu, BioAsu::Any);
        assert_eq!(omitted_connection.reported_sym, [0; 4]);
        assert_eq!(omitted_connection.reported_distance, 0.0);
        assert!(omitted_connection.link_id.is_empty());
        assert_eq!(omitted_connection.partner1.altloc(), 0);
    }

    #[test]
    fn bio_mmcif_c25_connection_address_presence_precedence_and_label_mapping_follow_source() {
        let mut null_auth = c25_connection_row();
        null_auth[0] = "null-auth-precedence";
        null_auth[2] = ".";
        null_auth[12] = ".";
        let null_doc = c25_connection_document(&[null_auth], &[]);
        let null_address =
            parse_mmcif_connections(&null_doc.blocks()[0], &MmcifAtomSiteGrouping::default())
                .unwrap()
                .remove(0)
                .partner1;
        // Column presence selects the auth branch even though the value is
        // null; valid label identifiers do not replace it.
        assert_eq!(
            null_address.chain_name(),
            PdbChainId::from_ascii(b"").unwrap()
        );
        assert_eq!(null_address.residue().sequence_number(), None);
        assert_eq!(null_address.residue().insertion_code(), None);
        assert_eq!(
            null_address.residue().name(),
            ResidueName::from_ascii(b"LIG").unwrap()
        );

        let mut label_row = c25_connection_row();
        label_row[0] = "label-addresses";
        label_row[15] = "3";
        let omitted_auth_columns = [2, 3, 12, 13];
        let label_doc = c25_connection_document(&[label_row], &omitted_auth_columns);
        let grouping = c11_atom_site_grouping();
        let label_connection = parse_mmcif_connections(&label_doc.blocks()[0], &grouping)
            .unwrap()
            .remove(0);
        assert_eq!(
            label_connection.partner1.chain_name(),
            PdbChainId::from_ascii(b"A").unwrap()
        );
        assert_eq!(
            label_connection.partner1.residue().sequence_number(),
            Some(1)
        );
        assert_eq!(label_connection.partner1.residue().insertion_code(), None);
        assert_eq!(label_connection.partner1.logical_atom_name(), "C1");
        assert_eq!(label_connection.partner1.altloc(), b'A');
        assert_eq!(
            label_connection.partner2.chain_name(),
            PdbChainId::from_ascii(b"A").unwrap()
        );
        assert_eq!(
            label_connection.partner2.residue().sequence_number(),
            Some(3)
        );
        assert_eq!(
            label_connection.partner2.residue().name(),
            ResidueName::from_ascii(b"GLY").unwrap()
        );
        assert_eq!(label_connection.partner2.logical_atom_name(), "N");
        assert_eq!(label_connection.partner2.altloc(), 0);
    }

    #[test]
    fn bio_mmcif_c25_connection_unmatched_label() {
        // Gemmi mmcif.cpp::set_part_of_address_from_label retains the most
        // recently matching subchain's chain name even when no sequence matches;
        // read_connectivity still fills the other address fields and row.
        let atom_rows = [
            c11_atom_site_row("1", "x", "e1", "1", "1", "A", "1"),
            c11_atom_site_row("2", "y", "e1", "3", "3", "A", "1"),
            c11_atom_site_row("3", "x", "e2", "1", "1", "B", "1"),
        ];
        let atom_doc = c05_scalar_document(&atom_rows, &[]);
        let grouping = group_mmcif_atom_sites(&atom_site_table(&atom_doc.blocks()[0])).unwrap();
        assert_eq!(grouping.models.len(), 1);

        let mut absent_first = c25_connection_row();
        absent_first[0] = "absent-first";
        absent_first[4] = "z";
        absent_first[14] = "99";
        absent_first[5] = "y";
        absent_first[15] = "3";
        let mut absent_second = c25_connection_row();
        absent_second[0] = "absent-second";
        absent_second[4] = "y";
        absent_second[14] = "3";
        absent_second[5] = "z";
        absent_second[15] = "99";
        let mut present_first = c25_connection_row();
        present_first[0] = "present-first";
        present_first[4] = "x";
        present_first[14] = "99";
        present_first[5] = "y";
        present_first[15] = "3";
        let mut present_second = c25_connection_row();
        present_second[0] = "present-second";
        present_second[4] = "y";
        present_second[14] = "3";
        present_second[5] = "x";
        present_second[15] = "99";
        // Omit the auth columns entirely: present-but-null auth columns select
        // Gemmi's auth branch and do not exercise label lookup.
        let doc = c25_connection_document(
            &[absent_first, absent_second, present_first, present_second],
            &[2, 3, 12, 13],
        );
        let connections = parse_mmcif_connections(&doc.blocks()[0], &grouping).unwrap();
        assert_eq!(connections.len(), 4);
        for (index, expected_name) in [
            "absent-first",
            "absent-second",
            "present-first",
            "present-second",
        ]
        .into_iter()
        .enumerate()
        {
            assert_eq!(connections[index].name, expected_name);
            let addresses = [&connections[index].partner1, &connections[index].partner2];
            for (partner, address) in addresses.into_iter().enumerate() {
                let unmatched = (index % 2) == partner;
                let expected_chain = if unmatched {
                    if index < 2 {
                        b"".as_slice()
                    } else {
                        b"B".as_slice()
                    }
                } else {
                    b"A".as_slice()
                };
                assert_eq!(
                    address.chain_name(),
                    PdbChainId::from_ascii(expected_chain).unwrap()
                );
                assert_eq!(
                    address.residue().sequence_number(),
                    (!unmatched).then_some(3)
                );
                assert_eq!(address.residue().insertion_code(), None);
                assert_eq!(
                    address.residue().name(),
                    ResidueName::from_ascii(if partner == 0 { b"LIG" } else { b"GLY" }).unwrap()
                );
                assert_eq!(
                    address.logical_atom_name(),
                    if partner == 0 { "C1" } else { "N" }
                );
                assert_eq!(address.altloc(), if partner == 0 { b'A' } else { 0 });
            }
        }
    }

    #[test]
    fn bio_mmcif_c25_connection_reports_source_address_and_representation_errors() {
        let mut missing_address = c25_connection_row();
        missing_address[0] = "missing-address";
        let missing_address_columns = [2, 3, 4, 5, 12, 13, 14, 15];
        let missing_doc = c25_connection_document(&[missing_address], &missing_address_columns);
        assert!(matches!(
            parse_mmcif_connections(&missing_doc.blocks()[0], &MmcifAtomSiteGrouping::default()),
            Err(MmcifConnectionError::MissingAddressIdentifiers { partner: 1 })
        ));

        let label_columns = [2, 3, 12, 13];
        let no_model_doc = c25_connection_document(&[c25_connection_row()], &label_columns);
        assert!(matches!(
            parse_mmcif_connections(&no_model_doc.blocks()[0], &MmcifAtomSiteGrouping::default()),
            Err(MmcifConnectionError::NoStructuralModels)
        ));

        let mut long_chain = c25_connection_row();
        long_chain[0] = "long-chain";
        long_chain[2] = "ABCDE";
        let long_chain_doc = c25_connection_document(&[long_chain], &[]);
        assert!(matches!(
            parse_mmcif_connections(
                &long_chain_doc.blocks()[0],
                &MmcifAtomSiteGrouping::default()
            ),
            Err(MmcifConnectionError::ChainNameNotRepresentable {
                partner: 1,
                value
            }) if value == "ABCDE"
        ));

        let mut long_residue = c25_connection_row();
        long_residue[0] = "long-residue";
        long_residue[6] = "ABCDE";
        let long_residue_doc = c25_connection_document(&[long_residue], &[]);
        assert!(matches!(
            parse_mmcif_connections(
                &long_residue_doc.blocks()[0],
                &MmcifAtomSiteGrouping::default()
            ),
            Err(MmcifConnectionError::ResidueNameNotRepresentable {
                partner: 1,
                value
            }) if value == "ABCDE"
        ));

        let mut invalid_sequence = c25_connection_row();
        invalid_sequence[0] = "invalid-sequence";
        invalid_sequence[12] = "12?";
        let invalid_sequence_doc = c25_connection_document(&[invalid_sequence], &[]);
        assert!(matches!(
            parse_mmcif_connections(
                &invalid_sequence_doc.blocks()[0],
                &MmcifAtomSiteGrouping::default()
            ),
            Err(MmcifConnectionError::SequenceId(
                MmcifSequenceIdError::InvalidSequenceNumber(error)
            )) if error.kind() == CifReadErrorKind::InvalidValue
        ));

        let mut invalid_altloc = c25_connection_row();
        invalid_altloc[0] = "invalid-altloc";
        invalid_altloc[10] = "AB";
        let invalid_altloc_doc = c25_connection_document(&[invalid_altloc], &[]);
        assert!(matches!(
            parse_mmcif_connections(
                &invalid_altloc_doc.blocks()[0],
                &MmcifAtomSiteGrouping::default()
            ),
            Err(MmcifConnectionError::Cif(error))
                if error.kind() == CifReadErrorKind::InvalidValue
        ));

        let required_id_omitted = c25_connection_document(&[c25_connection_row()], &[0]);
        assert!(
            parse_mmcif_connections(
                &required_id_omitted.blocks()[0],
                &MmcifAtomSiteGrouping::default()
            )
            .unwrap()
            .is_empty()
        );
    }

    #[test]
    fn bio_mmcif_c26_cis_peptide_rows_preserve_model_endpoints_altloc_and_angle() {
        // Gemmi mmcif.cpp::read_prot_cis appends in table order, prefers auth
        // component names, and leaves optional second-partner fields defaulted.
        let first = c26_cis_peptide_row();
        let mut second = c26_cis_peptide_row();
        second[0] = ".";
        second[2] = "7";
        second[5] = ".";
        second[6] = ".";
        second[7] = ".";
        second[10] = ".";
        second[11] = ".";
        second[12] = ".";
        let doc = c26_cis_peptide_document(&[first, second], &[]);
        let cis_peptides = parse_mmcif_cis_peptides(&doc.blocks()[0]).unwrap();
        assert_eq!(cis_peptides.len(), 2);

        let first = &cis_peptides[0];
        assert_eq!(first.model_num, 2);
        assert_eq!(
            first.partner_c.chain_name(),
            PdbChainId::from_ascii(b"A").unwrap()
        );
        assert_eq!(first.partner_c.residue().sequence_number(), Some(15));
        assert_eq!(first.partner_c.residue().insertion_code(), Some(b'A'));
        assert_eq!(
            first.partner_c.residue().name(),
            ResidueName::from_ascii(b"CYS").unwrap()
        );
        assert_eq!(first.partner_c.logical_atom_name(), "");
        assert_eq!(first.partner_c.altloc(), 0);
        assert_eq!(
            first.partner_n.chain_name(),
            PdbChainId::from_ascii(b"B").unwrap()
        );
        assert_eq!(first.partner_n.residue().sequence_number(), Some(16));
        assert_eq!(first.partner_n.residue().insertion_code(), Some(b'C'));
        assert_eq!(
            first.partner_n.residue().name(),
            ResidueName::from_ascii(b"ALA").unwrap()
        );
        assert_eq!(first.partner_n.logical_atom_name(), "");
        assert_eq!(first.partner_n.altloc(), 0);
        assert_eq!(first.only_altloc, b'X');
        assert_eq!(first.reported_angle, -2.5);

        let second = &cis_peptides[1];
        assert_eq!(second.model_num, 0);
        assert_eq!(second.partner_c.residue().sequence_number(), Some(7));
        assert_eq!(second.partner_c.residue().insertion_code(), None);
        assert_eq!(
            second.partner_c.residue().name(),
            ResidueName::from_ascii(b"SER").unwrap()
        );
        assert_eq!(
            second.partner_n.chain_name(),
            PdbChainId::from_ascii(b"").unwrap()
        );
        assert_eq!(second.partner_n.residue().sequence_number(), None);
        assert_eq!(second.partner_n.residue().insertion_code(), Some(b'C'));
        assert_eq!(
            second.partner_n.residue().name(),
            ResidueName::from_ascii(b"GLY").unwrap()
        );
        assert_eq!(second.only_altloc, 0);
        assert!(second.reported_angle.is_nan());
    }

    #[test]
    fn bio_mmcif_c26_cis_peptide_optional_columns_and_errors_follow_source() {
        // The optional tag columns are absent, not present-but-null. Gemmi's
        // row.has guards keep the second sequence, altloc, and angle defaults.
        let doc = c26_cis_peptide_document(&[c26_cis_peptide_row()], &[3, 5, 6, 7, 8, 10, 11, 12]);
        let row = parse_mmcif_cis_peptides(&doc.blocks()[0])
            .unwrap()
            .remove(0);
        assert_eq!(row.model_num, 2);
        assert_eq!(row.partner_c.residue().sequence_number(), Some(15));
        assert_eq!(row.partner_c.residue().insertion_code(), Some(b'A'));
        assert_eq!(
            row.partner_c.residue().name(),
            ResidueName::from_ascii(b"SER").unwrap()
        );
        assert_eq!(
            row.partner_n.chain_name(),
            PdbChainId::from_ascii(b"").unwrap()
        );
        assert_eq!(row.partner_n.residue().sequence_number(), None);
        assert_eq!(
            row.partner_n.residue().name(),
            ResidueName::from_ascii(b"GLY").unwrap()
        );
        assert_eq!(row.only_altloc, 0);
        assert!(row.reported_angle.is_nan());

        let mut bad_model = c26_cis_peptide_row();
        bad_model[0] = "bad";
        let doc = c26_cis_peptide_document(&[bad_model], &[]);
        assert!(matches!(parse_mmcif_cis_peptides(&doc.blocks()[0]),
            Err(MmcifCisPepError::Cif(error)) if error.kind() == CifReadErrorKind::InvalidValue));

        for (column, partner) in [(2, "C"), (7, "N")] {
            let mut bad_sequence = c26_cis_peptide_row();
            bad_sequence[column] = "12?";
            let doc = c26_cis_peptide_document(&[bad_sequence], &[]);
            assert!(matches!(parse_mmcif_cis_peptides(&doc.blocks()[0]),
                Err(MmcifCisPepError::SequenceId { partner: actual, .. }) if actual == partner));
        }

        let mut bad_insertion = c26_cis_peptide_row();
        bad_insertion[3] = "B";
        let doc = c26_cis_peptide_document(&[bad_insertion], &[]);
        assert!(matches!(
            parse_mmcif_cis_peptides(&doc.blocks()[0]),
            Err(MmcifCisPepError::SequenceId { partner: "C", .. })
        ));

        let mut bad_altloc = c26_cis_peptide_row();
        bad_altloc[11] = "AB";
        let doc = c26_cis_peptide_document(&[bad_altloc], &[]);
        assert!(matches!(parse_mmcif_cis_peptides(&doc.blocks()[0]),
            Err(MmcifCisPepError::Cif(error)) if error.kind() == CifReadErrorKind::InvalidValue));

        for (column, partner, chain) in [
            (1, "C", true),
            (6, "N", true),
            (5, "C", false),
            (10, "N", false),
        ] {
            let mut too_wide = c26_cis_peptide_row();
            too_wide[column] = "ABCDE";
            let doc = c26_cis_peptide_document(&[too_wide], &[]);
            assert!(
                matches!(parse_mmcif_cis_peptides(&doc.blocks()[0]),
                Err(MmcifCisPepError::ChainNameNotRepresentable { partner: actual, .. }) if chain && actual == partner)
                    || matches!(parse_mmcif_cis_peptides(&doc.blocks()[0]),
                    Err(MmcifCisPepError::ResidueNameNotRepresentable { partner: actual, .. }) if !chain && actual == partner)
            );
        }

        let missing_required = c26_cis_peptide_document(&[c26_cis_peptide_row()], &[0]);
        assert!(
            parse_mmcif_cis_peptides(&missing_required.blocks()[0])
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn bio_mmcif_c27_modified_residues_preserve_source_order_and_optional_text() {
        // Gemmi mmcif.cpp::read_struct_mod_residue assigns row.one_of(3,4)
        // RAW to the residue name, but row.str() decodes optional text fields.
        let first = c27_modified_residue_row();
        let mut second = c27_modified_residue_row();
        second[0] = "B";
        second[1] = "7";
        second[3] = ".";
        second[4] = "SEP";
        second[5] = ".";
        second[6] = "?";
        second[7] = ".";
        let doc = c27_modified_residue_document(&[first, second], &[]);
        let rows = parse_mmcif_modified_residues(&doc.blocks()[0]).unwrap();
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].chain_name, PdbChainId::from_ascii(b"A").unwrap());
        assert_eq!(rows[0].res_id.sequence_number(), Some(15));
        assert_eq!(rows[0].res_id.insertion_code(), Some(b'A'));
        assert_eq!(
            rows[0].res_id.name(),
            ResidueName::from_ascii(b"MSE").unwrap()
        );
        assert_eq!(rows[0].parent_comp_id, "MET");
        assert_eq!(rows[0].details, "oxidized sulfur");
        assert_eq!(rows[0].mod_id, "REFMAC1");
        assert_eq!(rows[1].chain_name, PdbChainId::from_ascii(b"B").unwrap());
        assert_eq!(rows[1].res_id.sequence_number(), Some(7));
        assert_eq!(rows[1].res_id.insertion_code(), None);
        assert_eq!(
            rows[1].res_id.name(),
            ResidueName::from_ascii(b"SEP").unwrap()
        );
        assert_eq!(rows[1].parent_comp_id, "");
        assert_eq!(rows[1].details, "");
        assert_eq!(rows[1].mod_id, "");

        let mut quoted = c27_modified_residue_row();
        quoted[3] = "'ME'";
        let quoted_doc = c27_modified_residue_document(&[quoted], &[]);
        let quoted_row = parse_mmcif_modified_residues(&quoted_doc.blocks()[0])
            .unwrap()
            .remove(0);
        assert_eq!(
            quoted_row.res_id.name(),
            ResidueName::from_ascii(b"'ME'").unwrap()
        );
    }

    #[test]
    fn bio_mmcif_c27_modified_residue_optional_columns_and_boundaries() {
        let absent_optional =
            c27_modified_residue_document(&[c27_modified_residue_row()], &[2, 3, 5, 6, 7]);
        let row = parse_mmcif_modified_residues(&absent_optional.blocks()[0])
            .unwrap()
            .remove(0);
        assert_eq!(row.res_id.sequence_number(), Some(15));
        assert_eq!(row.res_id.insertion_code(), Some(b'A'));
        assert_eq!(row.res_id.name(), ResidueName::from_ascii(b"SEC").unwrap());
        assert_eq!(row.parent_comp_id, "");
        assert_eq!(row.details, "");
        assert_eq!(row.mod_id, "");

        let mut raw_null_name = c27_modified_residue_row();
        raw_null_name[3] = ".";
        raw_null_name[4] = ".";
        let doc = c27_modified_residue_document(&[raw_null_name], &[]);
        let row = parse_mmcif_modified_residues(&doc.blocks()[0])
            .unwrap()
            .remove(0);
        assert_eq!(row.res_id.name(), ResidueName::from_ascii(b".").unwrap());

        let mut invalid_seq = c27_modified_residue_row();
        invalid_seq[1] = "12?";
        let doc = c27_modified_residue_document(&[invalid_seq], &[]);
        assert!(matches!(
            parse_mmcif_modified_residues(&doc.blocks()[0]),
            Err(MmcifModResError::SequenceId(_))
        ));

        let mut conflicting_code = c27_modified_residue_row();
        conflicting_code[2] = "B";
        let doc = c27_modified_residue_document(&[conflicting_code], &[]);
        assert!(matches!(
            parse_mmcif_modified_residues(&doc.blocks()[0]),
            Err(MmcifModResError::SequenceId(_))
        ));

        let mut long_chain = c27_modified_residue_row();
        long_chain[0] = "ABCDE";
        let doc = c27_modified_residue_document(&[long_chain], &[]);
        assert!(matches!(parse_mmcif_modified_residues(&doc.blocks()[0]),
            Err(MmcifModResError::ChainNameNotRepresentable { value }) if value == "ABCDE"));

        let mut long_name = c27_modified_residue_row();
        long_name[3] = "ABCDE";
        let doc = c27_modified_residue_document(&[long_name], &[]);
        assert!(matches!(parse_mmcif_modified_residues(&doc.blocks()[0]),
            Err(MmcifModResError::ResidueNameNotRepresentable { value }) if value == "ABCDE"));

        let missing_required = c27_modified_residue_document(&[c27_modified_residue_row()], &[0]);
        assert!(
            parse_mmcif_modified_residues(&missing_required.blocks()[0])
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn bio_mmcif_c28_operation_expr_preserves_ranges_order_and_first_group_only() {
        // Gemmi mmcif.cpp::parse_operation_expr expands numeric ranges but
        // explicitly ignores products after the first parenthesized group.
        let cases: [(&str, &[&str]); 8] = [
            ("3", &["3"]),
            ("1,3,5", &["1", "3", "5"]),
            ("one,two", &["one", "two"]),
            ("(3)", &["3"]),
            ("(a)", &["a"]),
            ("(1-3)", &["1", "2", "3"]),
            ("(2,3-4,XY)", &["2", "3", "4", "XY"]),
            ("(2,3)(4,5)", &["2", "3"]),
        ];
        for (input, expected) in cases {
            assert_eq!(
                parse_mmcif_operation_expr(input).unwrap(),
                expected
                    .iter()
                    .map(|text| (*text).to_owned())
                    .collect::<Vec<_>>(),
                "{input:?}"
            );
        }
    }

    #[test]
    fn bio_mmcif_c28_operation_expr_malformed_prefixes_follow_pinned_source() {
        // no_sign_atoi accepts whitespace/digit prefixes and produces zero
        // when no digits convert; parse_operation_expr does not validate
        // parentheses or reject empty tokens.
        let cases: [(&str, &[&str]); 9] = [
            ("", &[""]),
            ("(1,)", &["1", ""]),
            ("(1,2", &["1", "2"]),
            ("(3-1)", &[]),
            ("(1-)", &[]),
            ("(-2)", &["0", "1", "2"]),
            ("(x-y)", &["0"]),
            (" 2-4", &["2", "3", "4"]),
            ("1x-3", &["1", "2", "3"]),
        ];
        for (input, expected) in cases {
            assert_eq!(
                parse_mmcif_operation_expr(input).unwrap(),
                expected
                    .iter()
                    .map(|text| (*text).to_owned())
                    .collect::<Vec<_>>(),
                "{input:?}"
            );
        }
        assert_eq!(
            parse_mmcif_operation_expr("1-2147483647"),
            Err(MmcifOperationExprError::SourceUndefinedRangeIncrement),
        );
        assert_eq!(
            parse_mmcif_operation_expr("2147483648-2"),
            Err(MmcifOperationExprError::SourceUndefinedIntegerOverflow),
        );
    }

    #[test]
    fn bio_mmcif_c29_assemblies_preserve_operator_order_first_match_and_subchains() {
        // Gemmi read_assemblies reads operators first, scans each matching
        // generator in source order, and uses find_or_null's first ID match.
        let doc = c29_assembly_document(
            &[[
                "asm",
                "author_and_software_defined_assembly",
                "program",
                "tetramer",
                "4",
            ]],
            &[
                ["asm", "'ABSA (A^2)'", "10"],
                ["asm", "'SSA (A^2)'", "20"],
                ["asm", "MORE", "30"],
                ["asm", "'ABSA (A^2)'", "40"],
                ["other", "MORE", "999"],
            ],
            &[
                ["asm", "(2,1,missing,1)(x)", "A,B,,C,"],
                ["asm", "3", "Z"],
                ["other", "1", "Q"],
            ],
            &[
                ("1", "identity", ["1", "2", "3"]),
                ("2", "rotation", ["4", "5", "6"]),
                ("1", "duplicate", ["9", "9", "9"]),
            ],
        );
        let assemblies = parse_mmcif_assemblies(&doc.blocks()[0]).unwrap();
        assert_eq!(assemblies.len(), 1);
        let assembly = &assemblies[0];
        assert_eq!(assembly.name, "asm");
        assert!(assembly.author_determined);
        assert!(assembly.software_determined);
        assert_eq!(assembly.special_kind, BioAssemblySpecialKind::NotApplicable);
        assert_eq!(assembly.software_name, "program");
        assert_eq!(assembly.oligomeric_details, "tetramer");
        assert_eq!(assembly.oligomeric_count, 4);
        assert_eq!(
            (
                assembly.buried_surface_area,
                assembly.surface_area,
                assembly.solvent_free_energy_change
            ),
            (40.0, 20.0, 30.0)
        );
        assert_eq!(assembly.generators.len(), 2);
        let first = &assembly.generators[0];
        assert!(first.chains.is_empty());
        assert_eq!(first.subchains, ["A", "B", "", "C", ""]);
        assert_eq!(first.operators.len(), 3);
        assert_eq!(
            first
                .operators
                .iter()
                .map(|op| op.name.as_deref())
                .collect::<Vec<_>>(),
            [Some("2"), Some("1"), Some("1")]
        );
        assert_eq!(
            first
                .operators
                .iter()
                .map(|op| op.operator_type.as_deref())
                .collect::<Vec<_>>(),
            [Some("rotation"), Some("identity"), Some("identity")]
        );
        assert_eq!(*first.operators[0].transform.translation(), [4.0, 5.0, 6.0]);
        assert_eq!(*first.operators[1].transform.translation(), [1.0, 2.0, 3.0]);
        assert_eq!(*first.operators[2].transform.translation(), [1.0, 2.0, 3.0]);
        assert_eq!(assembly.generators[1].subchains, ["Z"]);
        assert!(assembly.generators[1].operators.is_empty());
    }

    #[test]
    fn bio_mmcif_c29_assemblies_filter_details_before_conversion_and_keep_defaults() {
        let doc = c29_assembly_document(
            &[
                ["skip", "unsupported_kind", ".", ".", "bad"],
                ["point", "'complete point assembly'", ".", ".", "2"],
                ["empty", ".", ".", ".", "."],
            ],
            &[],
            &[["point", "absent", "X,,"]],
            &[],
        );
        let rows = parse_mmcif_assemblies(&doc.blocks()[0]).unwrap();
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].name, "point");
        assert_eq!(rows[0].special_kind, BioAssemblySpecialKind::CompletePoint);
        assert!(!rows[0].author_determined);
        assert!(!rows[0].software_determined);
        assert_eq!(rows[0].oligomeric_count, 2);
        assert_eq!(rows[0].generators.len(), 1);
        assert_eq!(rows[0].generators[0].subchains, ["X", "", ""]);
        assert!(rows[0].generators[0].operators.is_empty());
        assert_eq!(rows[1].name, "empty");
        assert_eq!(rows[1].special_kind, BioAssemblySpecialKind::NotApplicable);
        assert_eq!(rows[1].oligomeric_count, 0);
        assert!(rows[1].buried_surface_area.is_nan());
        assert!(rows[1].surface_area.is_nan());
        assert!(rows[1].solvent_free_energy_change.is_nan());

        let bad_count = c29_assembly_document(
            &[["author", "author_defined_assembly", ".", ".", "bad"]],
            &[],
            &[],
            &[],
        );
        assert!(matches!(parse_mmcif_assemblies(&bad_count.blocks()[0]),
            Err(MmcifAssemblyError::Cif(error)) if error.kind() == CifReadErrorKind::InvalidValue));

        let bad_expression = c29_assembly_document(
            &[["author", "author_defined_assembly", ".", ".", "1"]],
            &[],
            &[["author", "1-2147483647", "A"]],
            &[],
        );
        assert!(matches!(
            parse_mmcif_assemblies(&bad_expression.blocks()[0]),
            Err(MmcifAssemblyError::OperationExpr(
                MmcifOperationExprError::SourceUndefinedRangeIncrement
            ))
        ));
    }
}
