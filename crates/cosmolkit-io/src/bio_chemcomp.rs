//! Private Gemmi-shaped chemical-component read helpers.

use crate::cif::{CifBlock, CifDocument, CifReadError, CifTable, cif_as_f64, cif_as_string};
use cosmolkit_bio::{
    AtomName, AtomSourceIds, BioAtomRow, BioCalcFlag, BioChainId, BioChainRow, BioCoordinateBlock,
    BioCoordinateFormat, BioMetadata, BioModelId, BioModelRow, BioResidueId, BioResidueRow,
    BioRowSpan, BioSiftsUnpResidue, BioStructureData, BioStructureError, BioStructureParts,
    BioStructureSourceState, ChainKind, ChainSourceIds, EntityKind, PdbAtomSerial, PdbSeqId,
    ResidueName, ResidueSourceIds, find_residue_info,
};
use std::error::Error;
use std::fmt;

#[derive(Debug)]
enum ChemCompAtomError {
    MissingValue { row: usize, column: usize },
    NameOutsideCanonicalBoundary { row: usize, value: String },
    ChargeOutsideSourceDefinedRange { row: usize, value: String },
    Number(CifReadError),
}

impl fmt::Display for ChemCompAtomError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingValue { row, column } => {
                write!(f, "chemcomp row {row} is missing column {column}")
            }
            Self::NameOutsideCanonicalBoundary { row, value } => write!(
                f,
                "chemcomp row {row} atom name is not representable: {value:?}"
            ),
            Self::ChargeOutsideSourceDefinedRange { row, value } => write!(
                f,
                "chemcomp row {row} charge is not representable: {value:?}"
            ),
            Self::Number(error) => error.fmt(f),
        }
    }
}

impl Error for ChemCompAtomError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::Number(error) => Some(error),
            _ => None,
        }
    }
}

#[derive(Debug)]
pub(super) enum ChemCompModelError {
    NotChemComp,
    Table(CifReadError),
    Atom(ChemCompAtomError),
    Hierarchy(BioStructureError),
    ResidueName(String),
}

impl fmt::Display for ChemCompModelError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NotChemComp => write!(f, "Not a chem_comp format."),
            Self::Table(error) => error.fmt(f),
            Self::Atom(error) => error.fmt(f),
            Self::Hierarchy(error) => error.fmt(f),
            Self::ResidueName(name) => {
                write!(f, "chemcomp residue name is not representable: {name:?}")
            }
        }
    }
}

impl Error for ChemCompModelError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Self::NotChemComp => None,
            Self::Table(error) => Some(error),
            Self::Atom(error) => Some(error),
            Self::Hierarchy(error) => Some(error),
            Self::ResidueName(_) => None,
        }
    }
}

fn append_chemcomp_model(
    parts: &mut BioStructureParts,
    block: &CifBlock,
    kind: ChemCompModel,
) -> Result<(), ChemCompModelError> {
    // Gemmi✔️✔️: inline Model make_model_from_chemcomp_block(const cif::Block& block, ChemCompModel kind) {
    // Gemmi✔️✔️:   Model model;
    // Gemmi✔️✔️:   model.chains.emplace_back("");
    // Gemmi✔️✔️:   model.chains[0].residues.push_back(make_residue_from_chemcomp_block(block, kind));
    // Gemmi✔️✔️:   return model;
    // Gemmi✔️✔️: }
    // Behavior: the source's default model number is zero until H06 renumbers
    // all models. The one empty chain and one residue are retained even if
    // the atom table is empty. All row IDs and coordinate rows are aligned.
    // Complexity: the atom table is parsed once and appended once. Existing
    // coordinate storage is moved, not cloned or revalidated on each append.
    let setup = chemcomp_residue_setup(block, kind).map_err(ChemCompModelError::Table)?;
    let parsed_atoms = chemcomp_atom_rows(&setup).map_err(ChemCompModelError::Atom)?;
    let residue_name = ResidueName::from_ascii(setup.name.as_bytes())
        .ok_or_else(|| ChemCompModelError::ResidueName(setup.name.clone()))?;
    let model_id = BioModelId::new(u32::try_from(parts.models.len()).map_err(|_| {
        ChemCompModelError::Hierarchy(BioStructureError::RowIndexTooLarge {
            value: parts.models.len(),
        })
    })?);
    let chain_id = BioChainId::new(u32::try_from(parts.chains.len()).map_err(|_| {
        ChemCompModelError::Hierarchy(BioStructureError::RowIndexTooLarge {
            value: parts.chains.len(),
        })
    })?);
    let residue_id = BioResidueId::new(u32::try_from(parts.residues.len()).map_err(|_| {
        ChemCompModelError::Hierarchy(BioStructureError::RowIndexTooLarge {
            value: parts.residues.len(),
        })
    })?);
    let atom_span = BioRowSpan::from_usize(parts.atoms.len(), parsed_atoms.len())
        .map_err(ChemCompModelError::Hierarchy)?;
    let residue_span =
        BioRowSpan::from_usize(parts.residues.len(), 1).map_err(ChemCompModelError::Hierarchy)?;
    let chain_span =
        BioRowSpan::from_usize(parts.chains.len(), 1).map_err(ChemCompModelError::Hierarchy)?;
    let source = ResidueSourceIds::new(
        Some(PdbSeqId::new(setup.sequence_number, None)),
        None,
        None,
        None,
        None,
    )
    .expect("empty segment ID is valid");

    let mut positions = std::mem::take(&mut parts.coordinates).into_positions();
    positions.reserve(parsed_atoms.len());
    for (atom, position) in parsed_atoms {
        parts.atoms.push(BioAtomRow::new(
            residue_id,
            atom.name(),
            atom.element(),
            atom.isotope_mass_number(),
            atom.altloc(),
            atom.formal_charge(),
            atom.calc_flag(),
            atom.occupancy(),
            atom.b_iso(),
            *atom.anisou(),
            atom.tls_group_id(),
            atom.fraction(),
            *atom.source(),
        ));
        positions.push(position);
    }
    parts.coordinates = BioCoordinateBlock::new(positions);
    parts.residues.push(BioResidueRow::new(
        chain_id,
        atom_span,
        residue_name,
        find_residue_info(residue_name.as_str()).kind,
        EntityKind::Unknown,
        None,
        None,
        source,
        BioSiftsUnpResidue::default(),
    ));
    parts.chains.push(BioChainRow::new(
        model_id,
        None,
        residue_span,
        ChainKind::Unknown,
        ChainSourceIds::default(),
    ));
    parts.models.push(BioModelRow::new(chain_span, Some(0)));
    Ok(())
}

pub(super) fn make_structure_from_chemcomp_block(
    block: &CifBlock,
    which: i32,
) -> Result<BioStructureData, ChemCompModelError> {
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
    // Behavior: mask and value-presence checks retain the source's fixed model
    // order. Source text is retained raw; only the completed detached value is
    // validated, once after all selected models have been materialized.
    // Complexity: three bounded tag lookups and one append per selected model;
    // no repeated whole-structure validation or atom-table re-scan is added.
    let mut source_state = BioStructureSourceState::default();
    if let Some(name) = block.find_value("_chem_comp.id") {
        source_state.name = name.raw().to_owned();
    }
    let mut parts = BioStructureParts {
        input_format: BioCoordinateFormat::ChemComp,
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
        metadata: BioMetadata::default(),
        source_state,
        coordinates: BioCoordinateBlock::default(),
        crystal: None,
        ncs_operators: Vec::new(),
        assemblies: Vec::new(),
    };
    for (bit, tag, kind) in [
        (1, "_chem_comp_atom.x", ChemCompModel::Xyz),
        (2, "_chem_comp_atom.model_Cartn_x", ChemCompModel::Example),
        (
            4,
            "_chem_comp_atom.pdbx_model_Cartn_x_ideal",
            ChemCompModel::Ideal,
        ),
    ] {
        if which & bit != 0 && block.has_any_value(tag) {
            append_chemcomp_model(&mut parts, block, kind)?;
        }
    }
    for (index, model) in parts.models.iter_mut().enumerate() {
        *model = BioModelRow::new(model.chain_span(), Some(index as i32 + 1));
    }
    BioStructureData::from_parts(parts).map_err(ChemCompModelError::Hierarchy)
}

fn chemcomp_atom_rows(
    setup: &ChemCompResidueSetup<'_>,
) -> Result<Vec<(BioAtomRow, [f64; 3])>, ChemCompAtomError> {
    // Gemmi✔️✔️: res.atoms.resize(table.length());
    // Gemmi✔️✔️: int n = 0;
    // Gemmi✔️✔️: for (auto row : table) {
    // Gemmi✔️✔️:   Atom& atom = res.atoms[n++];
    // Gemmi✔️✔️:   atom.name = row.str(0);
    // Gemmi✔️✔️:   atom.element = Element(row.str(1));
    // Gemmi✔️✔️:   if (row.has2(2))
    // Gemmi✔️✔️:     // Charge is defined as integer, but some cif files in the wild have
    // Gemmi✔️✔️:     // trailing '.000', so we read it as floating-point number.
    // Gemmi❗❌:     atom.charge = (signed char) std::round(cif::as_number(row[2]));
    // Gemmi❗❌:   atom.pos = Position(cif::as_number(row[3]),
    // Gemmi❗❌:                       cif::as_number(row[4]),
    // Gemmi❗❌:                       cif::as_number(row[5]));
    // Gemmi✔️✔️: }
    // Gemmi✔️✔️: char altloc = '\0'; // 0 if not set
    // Gemmi✔️✔️: signed char charge = 0;  // [-8, +8]
    // Gemmi✔️✔️: Element element = El::X;
    // Gemmi✔️✔️: CalcFlag calc_flag = CalcFlag::NotSet;  // mmCIF _atom_site.calc_flag
    // Gemmi✔️✔️: short tls_group_id = -1;
    // Gemmi✔️✔️: int serial = 0;
    // Gemmi✔️✔️: float fraction = 0.f;  // custom value, one use is Refmac's ccp4_deuterium_fraction
    // Gemmi✔️✔️: Position pos;
    // Gemmi✔️✔️: float occ = 1.0f;
    // Gemmi✔️✔️: float b_iso = 20.0f; // arbitrary default value
    // Gemmi✔️✔️: SMat33<float> aniso = {0, 0, 0, 0, 0, 0};
    // Behavior: typed rows retain exact decoded names, element/deuterium identity,
    // source defaults and coordinate numbers within the existing cif_as_f64
    // acceptance scope. Its ordinary finite bit parity is not yet proven;
    // out-of-range C++ float-to-signed-char conversion is undefined and the
    // canonical i8 boundary reports it explicitly.
    // Complexity: one pass and preallocated output, same O(row count) as source;
    // fixed-size row construction does not clone the surrounding structure.
    // Numeric conversion inherits cif_as_f64's known extra scans.
    let mut atoms = Vec::with_capacity(setup.atoms.len());
    for (index, row) in setup.atoms.iter().enumerate() {
        let raw_name = row
            .get(0)
            .ok_or(ChemCompAtomError::MissingValue {
                row: index,
                column: 0,
            })?
            .raw();
        let decoded_name = cif_as_string(raw_name);
        let name = AtomName::from_ascii(decoded_name.as_bytes()).ok_or_else(|| {
            ChemCompAtomError::NameOutsideCanonicalBoundary {
                row: index,
                value: decoded_name.clone(),
            }
        })?;
        let symbol = row.get(1).ok_or(ChemCompAtomError::MissingValue {
            row: index,
            column: 1,
        })?;
        let decoded_symbol = cif_as_string(symbol.raw());
        let symbol_bytes = decoded_symbol.as_bytes();
        let (element, isotope_mass_number) = super::bio_pdb::gemmi_find_element([
            symbol_bytes.first().copied().unwrap_or_default(),
            symbol_bytes.get(1).copied().unwrap_or_default(),
        ]);
        let charge = if row.has_value(2) {
            let raw = row
                .get(2)
                .expect("has_value ensures an existing value")
                .raw();
            let rounded = cif_as_f64(raw, f64::NAN)
                .map_err(ChemCompAtomError::Number)?
                .round();
            if !rounded.is_finite() || rounded < i8::MIN as f64 || rounded > i8::MAX as f64 {
                return Err(ChemCompAtomError::ChargeOutsideSourceDefinedRange {
                    row: index,
                    value: raw.to_owned(),
                });
            }
            rounded as i8
        } else {
            0
        };
        let mut position = [0.0; 3];
        for (axis, slot) in position.iter_mut().enumerate() {
            let raw = row
                .get(axis + 3)
                .ok_or(ChemCompAtomError::MissingValue {
                    row: index,
                    column: axis + 3,
                })?
                .raw();
            *slot = cif_as_f64(raw, f64::NAN).map_err(ChemCompAtomError::Number)?;
        }
        atoms.push((
            BioAtomRow::new(
                BioResidueId::new(0),
                name,
                element,
                isotope_mass_number,
                None,
                charge,
                BioCalcFlag::NotSet,
                1.0,
                20.0,
                [0.0; 6],
                -1,
                0.0,
                AtomSourceIds::new(Some(PdbAtomSerial::new(0))),
            ),
            position,
        ));
    }
    Ok(atoms)
}

struct ChemCompResidueSetup<'a> {
    sequence_number: i32,
    name: String,
    atoms: CifTable<'a>,
}

fn chemcomp_residue_setup(
    block: &CifBlock,
    kind: ChemCompModel,
) -> Result<ChemCompResidueSetup<'_>, CifReadError> {
    // Gemmi✔️✔️: Residue res;
    // Gemmi✔️✔️: res.seqid.num = 1;
    // Gemmi✔️✔️: cif::Column col =
    // Gemmi✔️✔️:   const_cast<cif::Block&>(block).find_values("_chem_comp_atom.comp_id");
    // Gemmi✔️✔️: if (col && col.length() > 0)
    // Gemmi✔️✔️:   res.name = col[0];
    // Gemmi✔️✔️: else
    // Gemmi✔️✔️:   res.name = block.name.substr(starts_with(block.name, "comp_") ? 5 : 0);
    // Gemmi✔️✔️: cif::Table table = const_cast<cif::Block&>(block).find("_chem_comp_atom.",
    // Gemmi✔️✔️:         {"atom_id", "type_symbol", "?charge",
    // Gemmi✔️✔️:          xyz_tags[0], xyz_tags[1], xyz_tags[2]});
    // Behavior: comp_id is copied as its raw CIF value, including quoted or
    // null spelling; only absence/zero length selects the block-name path.
    // Complexity: one column lookup, one table lookup and at most one name
    // allocation, matching the source's scan and string construction shape.
    let name = if let Some(column) = block.find_values("_chem_comp_atom.comp_id") {
        if let Some(first) = column.get(0) {
            first.raw().to_owned()
        } else {
            block
                .name()
                .strip_prefix("comp_")
                .unwrap_or(block.name())
                .to_owned()
        }
    } else {
        block
            .name()
            .strip_prefix("comp_")
            .unwrap_or(block.name())
            .to_owned()
    };
    let tags = chemcomp_xyz_tags(block, kind);
    let atoms = block.find(
        "_chem_comp_atom.",
        &[
            "atom_id",
            "type_symbol",
            "?charge",
            tags[0],
            tags[1],
            tags[2],
        ],
    )?;
    Ok(ChemCompResidueSetup {
        sequence_number: 1,
        name,
        atoms,
    })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ChemCompModel {
    Xyz,
    Example,
    Ideal,
    First,
}

fn chemcomp_xyz_tags(block: &CifBlock, mut kind: ChemCompModel) -> [&'static str; 3] {
    // Gemmi✔️✔️: std::array<std::string, 3> xyz_tags;
    // Gemmi✔️✔️: if (kind == ChemCompModel::First) {
    // Gemmi✔️✔️:   const cif::Item* loop_item =
    // Gemmi✔️✔️:       const_cast<cif::Block&>(block).find_loop("_chem_comp_atom.atom_id").item();
    // Gemmi✔️✔️:   if (loop_item && loop_item->type == cif::ItemType::Loop)
    // Gemmi✔️✔️:     for (const std::string& tag : loop_item->loop.tags) {
    // Gemmi✔️✔️:       std::string lctag = gemmi::to_lower(tag);
    // Gemmi✔️✔️:       if (lctag == "_chem_comp_atom.x") {
    // Gemmi✔️✔️:         kind = ChemCompModel::Xyz;
    // Gemmi✔️✔️:         break;
    // Gemmi✔️✔️:       } else if (lctag == "_chem_comp_atom.model_cartn_x") {
    // Gemmi✔️✔️:         kind = ChemCompModel::Example;
    // Gemmi✔️✔️:         break;
    // Gemmi✔️✔️:       } else if (lctag == "_chem_comp_atom.pdbx_model_cartn_x_ideal") {
    // Gemmi✔️✔️:         kind = ChemCompModel::Ideal;
    // Gemmi✔️✔️:         break;
    // Gemmi✔️✔️:       }
    // Gemmi✔️✔️:     }
    // Gemmi✔️✔️: }
    // Gemmi✔️✔️: switch (kind) {
    // Gemmi✔️✔️:   case ChemCompModel::Xyz:
    // Gemmi✔️✔️:     xyz_tags = {{"x", "y", "z"}};
    // Gemmi✔️✔️:     break;
    // Gemmi✔️✔️:   case ChemCompModel::Example:
    // Gemmi✔️✔️:     xyz_tags = {{"model_Cartn_x", "model_Cartn_y", "model_Cartn_z"}};
    // Gemmi✔️✔️:     break;
    // Gemmi✔️✔️:   case ChemCompModel::Ideal:
    // Gemmi✔️✔️:     xyz_tags = {{"pdbx_model_Cartn_x_ideal",
    // Gemmi✔️✔️:                    "pdbx_model_Cartn_y_ideal",
    // Gemmi✔️✔️:                    "pdbx_model_Cartn_z_ideal"}};
    // Gemmi✔️✔️:     break;
    // Gemmi✔️✔️:   default:
    // Gemmi✔️✔️:     break;
    // Gemmi✔️✔️: }
    // Behavior: First chooses the first recognized X tag in the atom-id loop's
    // original column order. An absent match leaves all three tags empty.
    // Complexity: one loop lookup and one linear column scan, with no extra
    // allocation or repeated search beyond the source's scan.
    if kind == ChemCompModel::First
        && let Some(loop_) = block.find_loop("_chem_comp_atom.atom_id")
    {
        for tag in loop_.tags() {
            if tag.eq_ignore_ascii_case("_chem_comp_atom.x") {
                kind = ChemCompModel::Xyz;
                break;
            } else if tag.eq_ignore_ascii_case("_chem_comp_atom.model_cartn_x") {
                kind = ChemCompModel::Example;
                break;
            } else if tag.eq_ignore_ascii_case("_chem_comp_atom.pdbx_model_cartn_x_ideal") {
                kind = ChemCompModel::Ideal;
                break;
            }
        }
    }
    match kind {
        ChemCompModel::Xyz => ["x", "y", "z"],
        ChemCompModel::Example => ["model_Cartn_x", "model_Cartn_y", "model_Cartn_z"],
        ChemCompModel::Ideal => [
            "pdbx_model_Cartn_x_ideal",
            "pdbx_model_Cartn_y_ideal",
            "pdbx_model_Cartn_z_ideal",
        ],
        ChemCompModel::First => ["", "", ""],
    }
}

pub(super) fn check_chemcomp_block_number(document: &CifDocument) -> Option<usize> {
    // Gemmi✔️✔️: inline int check_chemcomp_block_number(const cif::Document& doc) {
    // Gemmi✔️✔️:   // monomer library file without global_
    // Gemmi✔️✔️:   if (doc.blocks.size() == 2 && doc.blocks[0].name == "comp_list")
    // Gemmi✔️✔️:     return 1;
    // Gemmi✔️✔️:   // monomer library file with global_
    // Gemmi✔️✔️:   if (doc.blocks.size() == 3 && doc.blocks[0].name.empty() &&
    // Gemmi✔️✔️:       doc.blocks[1].name == "comp_list")
    // Gemmi✔️✔️:     return 2;
    // Gemmi✔️✔️:   // CCD file
    // Gemmi✔️✔️:   if (doc.blocks.size() == 1 &&
    // Gemmi✔️✔️:       !doc.blocks[0].has_tag("_atom_site.id") &&
    // Gemmi✔️✔️:       !doc.blocks[0].has_tag("_cell.length_a") &&
    // Gemmi✔️✔️:       doc.blocks[0].has_tag("_chem_comp_atom.atom_id"))
    // Gemmi✔️✔️:     return 0;
    // Gemmi✔️✔️:   return -1;
    // Gemmi✔️✔️: }
    // Behavior: exact block-count/name branches before CCD tag screening;
    // Option::None is the private equivalent of the source's -1 sentinel.
    // Complexity: fixed number of indexed block/tag lookups, as in source.
    let blocks = document.blocks();
    if blocks.len() == 2 && blocks[0].name() == "comp_list" {
        return Some(1);
    }
    if blocks.len() == 3 && blocks[0].name().is_empty() && blocks[1].name() == "comp_list" {
        return Some(2);
    }
    if blocks.len() == 1
        && !blocks[0].has_tag("_atom_site.id")
        && !blocks[0].has_tag("_cell.length_a")
        && blocks[0].has_tag("_chem_comp_atom.atom_id")
    {
        return Some(0);
    }
    None
}

pub(super) fn make_structure_from_chemcomp_doc(
    document: &CifDocument,
    which: i32,
) -> Result<BioStructureData, ChemCompModelError> {
    // Gemmi✔️✔️: inline Structure make_structure_from_chemcomp_doc(const cif::Document& doc,
    // Gemmi✔️✔️:                                                   cif::Document* save_doc=nullptr,
    // Gemmi✔️✔️:                                                   int which=7) {
    // Gemmi✔️✔️:   int n = check_chemcomp_block_number(doc);
    // Gemmi✔️✔️:   if (n == -1)
    // Gemmi✔️✔️:     fail("Not a chem_comp format.");
    // Gemmi✔️✔️:   Structure st = make_structure_from_chemcomp_block(doc.blocks[n], which);
    // Gemmi❌❌:   if (save_doc)
    // Gemmi❌❌:     *save_doc = std::move(doc);
    // Gemmi✔️✔️:   return st;
    // Gemmi✔️✔️: }
    // Behavior: the private reader returns a structured not-chemcomp error
    // and reuses the exact H01-selected block; save_doc is not requested by
    // this Rust value contract and remains unmodeled independently.
    // Complexity: fixed branch selection plus the H06 model construction;
    // the borrowed document is not cloned or re-parsed.
    let index = check_chemcomp_block_number(document).ok_or(ChemCompModelError::NotChemComp)?;
    make_structure_from_chemcomp_block(&document.blocks()[index], which)
}

#[cfg(test)]
mod tests {
    use super::{
        ChemCompAtomError, ChemCompModel, ChemCompModelError, append_chemcomp_model,
        check_chemcomp_block_number, chemcomp_atom_rows, chemcomp_residue_setup, chemcomp_xyz_tags,
        make_structure_from_chemcomp_block, make_structure_from_chemcomp_doc,
    };
    use crate::cif::{CifCheckLevel, read_cif_document};
    use cosmolkit_bio::{
        BioCalcFlag, BioCoordinateBlock, BioCoordinateFormat, BioMetadata, BioStructureParts,
        BioStructureSourceState,
    };
    use cosmolkit_types::Element;

    fn selected(text: &str) -> Option<usize> {
        let document = read_cif_document(text, "h01.cif", CifCheckLevel::Syntax).unwrap();
        check_chemcomp_block_number(&document)
    }

    #[test]
    fn bio_read_h01_exact_one_two_three_block_shapes() {
        assert_eq!(selected("data_CCD\n_chem_comp_atom.atom_id C1\n"), Some(0));
        assert_eq!(selected("data_comp_list\ndata_CCD\n"), Some(1));
        assert_eq!(selected("global_\ndata_comp_list\ndata_CCD\n"), Some(2));
    }

    #[test]
    fn bio_read_h01_exclusions_and_first_blocks_are_source_ordered() {
        assert_eq!(
            selected("data_CCD\n_chem_comp_atom.atom_id C1\n_atom_site.id 1\n"),
            None
        );
        assert_eq!(
            selected("data_CCD\n_chem_comp_atom.atom_id C1\n_cell.length_a 1\n"),
            None
        );
        assert_eq!(selected("data_CCD\n_chem_comp.id CCD\n"), None);
        assert_eq!(
            selected("data_other\ndata_CCD\n_chem_comp_atom.atom_id C1\n"),
            None
        );
        assert_eq!(selected("data_comp_list\ndata_CCD\ndata_extra\n"), None);
        assert_eq!(
            selected("global_\ndata_other\ndata_CCD\n_chem_comp_atom.atom_id C1\n"),
            None
        );
    }

    fn tags(text: &str, kind: ChemCompModel) -> [&'static str; 3] {
        let document = read_cif_document(text, "h02.cif", CifCheckLevel::Syntax).unwrap();
        chemcomp_xyz_tags(&document.blocks()[0], kind)
    }

    #[test]
    fn bio_read_h02_explicit_kind_selects_its_exact_triplet() {
        let text = "data_comp_CMP\n_chem_comp_atom.atom_id A1\n";
        assert_eq!(tags(text, ChemCompModel::Xyz), ["x", "y", "z"]);
        assert_eq!(
            tags(text, ChemCompModel::Example),
            ["model_Cartn_x", "model_Cartn_y", "model_Cartn_z"]
        );
        assert_eq!(
            tags(text, ChemCompModel::Ideal),
            [
                "pdbx_model_Cartn_x_ideal",
                "pdbx_model_Cartn_y_ideal",
                "pdbx_model_Cartn_z_ideal"
            ]
        );
    }

    #[test]
    fn bio_read_h02_first_uses_actual_atom_id_loop_column_order_case_insensitively() {
        let first_example = "data_comp_CMP\nloop_\n_CHEM_COMP_ATOM.ATOM_ID\n_CHEM_COMP_ATOM.MODEL_CARTN_X\n_chem_comp_atom.x\nA1 1 2\n";
        assert_eq!(
            tags(first_example, ChemCompModel::First),
            ["model_Cartn_x", "model_Cartn_y", "model_Cartn_z"]
        );
        let first_xyz = "data_comp_CMP\nloop_\n_chem_comp_atom.atom_id\n_CHEM_COMP_ATOM.X\n_chem_comp_atom.model_cartn_x\nA1 1 2\n";
        assert_eq!(tags(first_xyz, ChemCompModel::First), ["x", "y", "z"]);
        let first_ideal = "data_comp_CMP\nloop_\n_chem_comp_atom.atom_id\n_CHEM_COMP_ATOM.PDBX_MODEL_CARTN_X_IDEAL\n_chem_comp_atom.x\nA1 1 2\n";
        assert_eq!(
            tags(first_ideal, ChemCompModel::First),
            [
                "pdbx_model_Cartn_x_ideal",
                "pdbx_model_Cartn_y_ideal",
                "pdbx_model_Cartn_z_ideal"
            ]
        );
    }

    #[test]
    fn bio_read_h02_first_without_matching_atom_id_loop_retains_empty_tags() {
        let pair = "data_comp_CMP\n_chem_comp_atom.atom_id A1\n_chem_comp_atom.x 1\n";
        assert_eq!(tags(pair, ChemCompModel::First), ["", "", ""]);
        let unrelated = "data_comp_CMP\nloop_\n_chem_comp_atom.atom_id\n_chem_comp_atom.y\nA1 1\n";
        assert_eq!(tags(unrelated, ChemCompModel::First), ["", "", ""]);
    }

    #[test]
    fn bio_read_h03_name_uses_raw_comp_id_then_source_block_fallback() {
        let with_comp = "data_comp_CMP\n_chem_comp_atom.comp_id 'RAW'\n";
        let document = read_cif_document(with_comp, "h03.cif", CifCheckLevel::Syntax).unwrap();
        let setup = chemcomp_residue_setup(&document.blocks()[0], ChemCompModel::Xyz).unwrap();
        assert_eq!(setup.sequence_number, 1);
        assert_eq!(setup.name, "'RAW'");

        let empty_value = "data_comp_CMP\n_chem_comp_atom.comp_id ''\n";
        let document = read_cif_document(empty_value, "h03.cif", CifCheckLevel::Syntax).unwrap();
        assert_eq!(
            chemcomp_residue_setup(&document.blocks()[0], ChemCompModel::Xyz)
                .unwrap()
                .name,
            "''"
        );

        let without_comp = "data_comp_CMP\n_chem_comp.id CMP\n";
        let document = read_cif_document(without_comp, "h03.cif", CifCheckLevel::Syntax).unwrap();
        assert_eq!(
            chemcomp_residue_setup(&document.blocks()[0], ChemCompModel::Xyz)
                .unwrap()
                .name,
            "CMP"
        );
        let document = read_cif_document("data_CMP\n", "h03.cif", CifCheckLevel::Syntax).unwrap();
        assert_eq!(
            chemcomp_residue_setup(&document.blocks()[0], ChemCompModel::Xyz)
                .unwrap()
                .name,
            "CMP"
        );
    }

    #[test]
    fn bio_read_h03_pair_and_loop_tables_require_fields_but_allow_optional_charge() {
        let pair = "data_comp_CMP\n_chem_comp_atom.atom_id 'A 1'\n_chem_comp_atom.type_symbol C\n_chem_comp_atom.x 1\n_chem_comp_atom.y 2\n_chem_comp_atom.z 3\n";
        let document = read_cif_document(pair, "h03.cif", CifCheckLevel::Syntax).unwrap();
        let setup = chemcomp_residue_setup(&document.blocks()[0], ChemCompModel::Xyz).unwrap();
        assert_eq!(setup.atoms.len(), 1);
        assert!(!setup.atoms.has_column(2));
        let row = setup.atoms.row(0).unwrap();
        assert_eq!(row.get(0).unwrap().raw(), "'A 1'");
        assert_eq!(row.decoded(0).as_deref(), Some("A 1"));

        let loop_text = "data_comp_CMP\nloop_\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n_chem_comp_atom.charge\n_chem_comp_atom.x\n_chem_comp_atom.y\n_chem_comp_atom.z\nA1 C 2 1 2 3\nA2 N -1 4 5 6\n";
        let document = read_cif_document(loop_text, "h03.cif", CifCheckLevel::Syntax).unwrap();
        let setup = chemcomp_residue_setup(&document.blocks()[0], ChemCompModel::Xyz).unwrap();
        assert_eq!(setup.atoms.len(), 2);
        assert!(setup.atoms.has_column(2));
        assert_eq!(setup.atoms.row(1).unwrap().get(2).unwrap().raw(), "-1");

        let missing = pair.replace("_chem_comp_atom.type_symbol C\n", "");
        let document = read_cif_document(&missing, "h03.cif", CifCheckLevel::Syntax).unwrap();
        let setup = chemcomp_residue_setup(&document.blocks()[0], ChemCompModel::Xyz).unwrap();
        assert_eq!(setup.atoms.len(), 0);
        assert!(!setup.atoms.is_present());
    }

    fn atom_rows(
        text: &str,
        kind: ChemCompModel,
    ) -> Result<Vec<(cosmolkit_bio::BioAtomRow, [f64; 3])>, ChemCompAtomError> {
        let document = read_cif_document(text, "h04.cif", CifCheckLevel::Syntax).unwrap();
        let setup = chemcomp_residue_setup(&document.blocks()[0], kind).unwrap();
        chemcomp_atom_rows(&setup)
    }

    #[test]
    fn bio_read_h04_all_coordinate_kinds_and_source_atom_defaults() {
        for (kind, columns) in [
            (ChemCompModel::Xyz, ["x", "y", "z"]),
            (
                ChemCompModel::Example,
                ["model_Cartn_x", "model_Cartn_y", "model_Cartn_z"],
            ),
            (
                ChemCompModel::Ideal,
                [
                    "pdbx_model_Cartn_x_ideal",
                    "pdbx_model_Cartn_y_ideal",
                    "pdbx_model_Cartn_z_ideal",
                ],
            ),
        ] {
            let text = format!(
                "data_comp_CMP\n_chem_comp_atom.atom_id CA\n_chem_comp_atom.type_symbol C\n_chem_comp_atom.{} 1.25\n_chem_comp_atom.{} -0.0\n_chem_comp_atom.{} 3e2\n",
                columns[0], columns[1], columns[2]
            );
            let rows = atom_rows(&text, kind).unwrap();
            assert_eq!(rows.len(), 1);
            let (atom, position) = &rows[0];
            assert_eq!(atom.name().as_str(), "CA");
            assert_eq!(atom.element(), Element::C);
            assert_eq!(atom.isotope_mass_number(), None);
            assert_eq!(atom.altloc(), None);
            assert_eq!(atom.formal_charge(), 0);
            assert_eq!(atom.calc_flag(), BioCalcFlag::NotSet);
            assert_eq!(atom.occupancy(), 1.0);
            assert_eq!(atom.b_iso(), 20.0);
            assert_eq!(*atom.anisou(), [0.0; 6]);
            assert_eq!(atom.tls_group_id(), -1);
            assert_eq!(atom.fraction(), 0.0);
            assert_eq!(atom.source().serial().unwrap().value(), 0);
            assert_eq!(position[0], 1.25);
            assert_eq!(position[1].to_bits(), (-0.0_f64).to_bits());
            assert_eq!(position[2], 300.0);
        }
    }

    #[test]
    fn bio_read_h04_charge_rounding_missing_markers_and_deuterium() {
        let text = "data_comp_CMP\nloop_\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n_chem_comp_atom.charge\n_chem_comp_atom.x\n_chem_comp_atom.y\n_chem_comp_atom.z\nH1 D 1.5 . 2 3\nC1 C -1.5 4 ? 6\nN1 N . 7 8 9\nO1 O 2.000 1 2 3\n";
        let rows = atom_rows(text, ChemCompModel::Xyz).unwrap();
        assert_eq!(rows.len(), 4);
        assert_eq!(rows[0].0.element(), Element::H);
        assert_eq!(rows[0].0.isotope_mass_number(), Some(2));
        assert_eq!(rows[0].0.formal_charge(), 2);
        assert!(rows[0].1[0].is_nan());
        assert_eq!(rows[1].0.formal_charge(), -2);
        assert!(rows[1].1[1].is_nan());
        assert_eq!(rows[2].0.formal_charge(), 0);
        assert_eq!(rows[3].0.formal_charge(), 2);
    }

    #[test]
    fn bio_read_h04_checked_name_and_charge_boundaries() {
        let long_name = "data_comp_CMP\n_chem_comp_atom.atom_id LONG5\n_chem_comp_atom.type_symbol C\n_chem_comp_atom.x 1\n_chem_comp_atom.y 2\n_chem_comp_atom.z 3\n";
        assert!(matches!(
            atom_rows(long_name, ChemCompModel::Xyz),
            Err(ChemCompAtomError::NameOutsideCanonicalBoundary { row: 0, .. })
        ));
        let charge = "data_comp_CMP\n_chem_comp_atom.atom_id C1\n_chem_comp_atom.type_symbol C\n_chem_comp_atom.charge 128\n_chem_comp_atom.x 1\n_chem_comp_atom.y 2\n_chem_comp_atom.z 3\n";
        assert!(matches!(
            atom_rows(charge, ChemCompModel::Xyz),
            Err(ChemCompAtomError::ChargeOutsideSourceDefinedRange { row: 0, .. })
        ));
    }

    fn empty_chemcomp_parts() -> BioStructureParts {
        BioStructureParts {
            input_format: BioCoordinateFormat::ChemComp,
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
            metadata: BioMetadata::default(),
            source_state: BioStructureSourceState::default(),
            coordinates: BioCoordinateBlock::default(),
            crystal: None,
            ncs_operators: Vec::new(),
            assemblies: Vec::new(),
        }
    }

    #[test]
    fn bio_read_h05_empty_atom_table_still_appends_one_empty_chain_and_residue() {
        let document =
            read_cif_document("data_comp_CMP\n", "h05.cif", CifCheckLevel::Syntax).unwrap();
        let mut parts = empty_chemcomp_parts();
        append_chemcomp_model(&mut parts, &document.blocks()[0], ChemCompModel::Xyz).unwrap();
        assert_eq!(
            (
                parts.models.len(),
                parts.chains.len(),
                parts.residues.len(),
                parts.atoms.len()
            ),
            (1, 1, 1, 0)
        );
        assert_eq!(
            (
                parts.models[0].chain_span().start(),
                parts.models[0].chain_span().len()
            ),
            (0, 1)
        );
        assert_eq!(
            (
                parts.chains[0].residue_span().start(),
                parts.chains[0].residue_span().len()
            ),
            (0, 1)
        );
        assert_eq!(
            (
                parts.residues[0].atom_span().start(),
                parts.residues[0].atom_span().len()
            ),
            (0, 0)
        );
        assert_eq!(parts.residues[0].name().as_str(), "CMP");
        assert_eq!(parts.residues[0].source().seq_id().unwrap().seq_num(), 1);
        assert!(parts.coordinates.is_empty());
    }

    #[test]
    fn bio_read_h05_multiple_appends_keep_order_ids_spans_and_coordinates_aligned() {
        let text = "data_comp_CMP\nloop_\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n_chem_comp_atom.x\n_chem_comp_atom.y\n_chem_comp_atom.z\nC1 C 1 2 3\nO1 O 4 5 6\n";
        let document = read_cif_document(text, "h05.cif", CifCheckLevel::Syntax).unwrap();
        let mut parts = empty_chemcomp_parts();
        for _ in 0..2 {
            append_chemcomp_model(&mut parts, &document.blocks()[0], ChemCompModel::Xyz).unwrap();
        }
        assert_eq!(
            (
                parts.models.len(),
                parts.chains.len(),
                parts.residues.len(),
                parts.atoms.len()
            ),
            (2, 2, 2, 4)
        );
        assert_eq!(
            parts
                .models
                .iter()
                .map(|row| (
                    row.chain_span().start(),
                    row.chain_span().len(),
                    row.source_model_number()
                ))
                .collect::<Vec<_>>(),
            vec![(0, 1, Some(0)), (1, 1, Some(0))]
        );
        assert_eq!(
            parts
                .chains
                .iter()
                .map(|row| (
                    row.model_id().value(),
                    row.residue_span().start(),
                    row.residue_span().len()
                ))
                .collect::<Vec<_>>(),
            vec![(0, 0, 1), (1, 1, 1)]
        );
        assert_eq!(
            parts
                .residues
                .iter()
                .map(|row| (
                    row.chain_id().value(),
                    row.atom_span().start(),
                    row.atom_span().len()
                ))
                .collect::<Vec<_>>(),
            vec![(0, 0, 2), (1, 2, 2)]
        );
        assert_eq!(
            parts
                .atoms
                .iter()
                .map(|row| (row.residue_id().value(), row.name().as_str().to_owned()))
                .collect::<Vec<_>>(),
            vec![
                (0, "C1".to_owned()),
                (0, "O1".to_owned()),
                (1, "C1".to_owned()),
                (1, "O1".to_owned())
            ]
        );
        assert_eq!(
            parts.coordinates.positions(),
            &[
                [1.0, 2.0, 3.0],
                [4.0, 5.0, 6.0],
                [1.0, 2.0, 3.0],
                [4.0, 5.0, 6.0]
            ]
        );
    }

    #[test]
    fn bio_read_h06_all_masks_cross_all_coordinate_set_presence() {
        // The pinned source checks each bit and has_any_value in fixed Xyz,
        // Example, Ideal order, independent of input tag order.
        for present in 0..8 {
            let x = if present & 1 != 0 { "1" } else { "." };
            let example = if present & 2 != 0 { "4" } else { "?" };
            let ideal = if present & 4 != 0 { "7" } else { "." };
            let text = format!(
                "data_comp_CMP\n_chem_comp.id 'RAW'\nloop_\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n_chem_comp_atom.pdbx_model_Cartn_x_ideal\n_chem_comp_atom.pdbx_model_Cartn_y_ideal\n_chem_comp_atom.pdbx_model_Cartn_z_ideal\n_chem_comp_atom.model_Cartn_x\n_chem_comp_atom.model_Cartn_y\n_chem_comp_atom.model_Cartn_z\n_chem_comp_atom.x\n_chem_comp_atom.y\n_chem_comp_atom.z\nC1 C {ideal} 8 9 {example} 5 6 {x} 2 3\n"
            );
            let document = read_cif_document(&text, "h06.cif", CifCheckLevel::Syntax).unwrap();
            for mask in 0..8 {
                let result =
                    make_structure_from_chemcomp_block(&document.blocks()[0], mask).unwrap();
                let selected = [
                    (1, [1.0, 2.0, 3.0]),
                    (2, [4.0, 5.0, 6.0]),
                    (4, [7.0, 8.0, 9.0]),
                ]
                .into_iter()
                .filter(|(bit, _)| (present & mask & bit) != 0)
                .map(|(_, xyz)| xyz)
                .collect::<Vec<_>>();
                assert_eq!(result.input_format(), BioCoordinateFormat::ChemComp);
                assert_eq!(result.source_state().name, "'RAW'");
                assert_eq!(
                    result.models().len(),
                    selected.len(),
                    "present={present}, mask={mask}"
                );
                assert_eq!(result.chains().len(), selected.len());
                assert_eq!(result.residues().len(), selected.len());
                assert_eq!(result.atoms().len(), selected.len());
                assert_eq!(result.coordinates().positions(), selected.as_slice());
                for (index, model) in result.models().iter().enumerate() {
                    assert_eq!(model.source_model_number(), Some(index as i32 + 1));
                    assert_eq!(
                        (model.chain_span().start(), model.chain_span().len()),
                        (index as u32, 1)
                    );
                    assert_eq!(result.atoms()[index].residue_id().value(), index as u32);
                }
            }
        }
    }

    #[test]
    fn bio_read_h06_null_only_tags_and_empty_result_keep_raw_source_name() {
        let text = "data_comp_CMP\n_chem_comp.id 'AB C'\n_chem_comp_atom.x .\n_chem_comp_atom.model_Cartn_x ?\n_chem_comp_atom.pdbx_model_Cartn_x_ideal .\n";
        let document = read_cif_document(text, "h06.cif", CifCheckLevel::Syntax).unwrap();
        let result = make_structure_from_chemcomp_block(&document.blocks()[0], 7).unwrap();
        assert!(result.models().is_empty());
        assert!(result.coordinates().is_empty());
        assert_eq!(result.source_state().name, "'AB C'");
    }

    #[test]
    fn bio_read_h07_uses_exact_chemcomp_block_in_all_three_document_shapes() {
        let atom_block = "data_comp_CMP\n_chem_comp.id 'SELECTED'\n_chem_comp_atom.atom_id C1\n_chem_comp_atom.type_symbol C\n_chem_comp_atom.x 1\n_chem_comp_atom.y 2\n_chem_comp_atom.z 3\n";
        for prefix in [
            "",
            "data_comp_list\n_chem_comp.id WRONG\n",
            "global_\ndata_comp_list\n_chem_comp.id WRONG\n",
        ] {
            let text = format!("{prefix}{atom_block}");
            let document = read_cif_document(&text, "h07.cif", CifCheckLevel::Syntax).unwrap();
            let result = make_structure_from_chemcomp_doc(&document, 1).unwrap();
            assert_eq!(result.source_state().name, "'SELECTED'");
            assert_eq!(result.input_format(), BioCoordinateFormat::ChemComp);
            assert_eq!(result.models().len(), 1);
            assert_eq!(result.residues()[0].name().as_str(), "CMP");
            assert_eq!(result.coordinates().positions(), &[[1.0, 2.0, 3.0]]);
        }
    }

    #[test]
    fn bio_read_h07_rejects_all_nonchemcomp_document_shapes_structurally() {
        for text in [
            "data_other\n_chem_comp.id OTHER\n",
            "data_comp_CMP\n_chem_comp_atom.atom_id C1\n_atom_site.id 1\n",
            "data_comp_CMP\n_chem_comp_atom.atom_id C1\n_cell.length_a 1\n",
            "data_other\ndata_comp_CMP\n_chem_comp_atom.atom_id C1\n",
            "data_comp_list\ndata_comp_CMP\ndata_extra\n",
        ] {
            let document = read_cif_document(text, "h07.cif", CifCheckLevel::Syntax).unwrap();
            assert!(matches!(
                make_structure_from_chemcomp_doc(&document, 7),
                Err(ChemCompModelError::NotChemComp)
            ));
        }
    }
}
