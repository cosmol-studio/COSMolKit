use std::collections::BTreeMap;

use cosmolkit_core::{RingInfo, fast_find_rings_from_parts};
use cosmolkit_model::{
    AtomId, Bond, BondId, Conformer2D, Conformer3D, SGroupConnection, StereoGroup, StereoGroupKind,
    SubstanceGroup, SubstanceGroupKind,
};
use cosmolkit_types::{BondDirection, BondOrder, BondStereo};

use crate::{SmilesParseError, SmilesRecord, SmilesWriteParams, writer::write_smiles_for_cx};

/// CXSMILES extension families selected for writing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CxSmilesFields(u32);

impl CxSmilesFields {
    pub const NONE: Self = Self(0);
    pub const ATOM_LABELS: Self = Self(1 << 0);
    pub const MOLFILE_VALUES: Self = Self(1 << 1);
    pub const COORDS: Self = Self(1 << 2);
    pub const RADICALS: Self = Self(1 << 3);
    pub const ATOM_PROPS: Self = Self(1 << 4);
    pub const LINKNODES: Self = Self(1 << 5);
    pub const ENHANCED_STEREO: Self = Self(1 << 6);
    pub const SGROUPS: Self = Self(1 << 7);
    pub const POLYMER: Self = Self(1 << 8);
    pub const BOND_CFG: Self = Self(1 << 9);
    pub const BOND_ATROPISOMER: Self = Self(1 << 10);
    pub const COORDINATE_BONDS: Self = Self(1 << 11);
    pub const HYDROGEN_BONDS: Self = Self(1 << 12);
    pub const ZERO_BONDS: Self = Self(1 << 13);
    pub const ALL: Self = Self(0x7fff_ffff);
    pub const ALL_BUT_COORDS: Self = Self(Self::ALL.0 ^ Self::COORDS.0);

    #[must_use]
    pub const fn bits(self) -> u32 {
        self.0
    }

    #[must_use]
    pub const fn contains(self, other: Self) -> bool {
        self.0 & other.0 == other.0
    }
}

impl std::ops::BitOr for CxSmilesFields {
    type Output = Self;

    fn bitor(self, rhs: Self) -> Self {
        Self(self.0 | rhs.0)
    }
}

/// Options for detached CXSMILES writing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CxSmilesWriteParams {
    pub smiles: SmilesWriteParams,
    pub fields: CxSmilesFields,
}

impl Default for CxSmilesWriteParams {
    fn default() -> Self {
        Self {
            smiles: SmilesWriteParams::default(),
            fields: CxSmilesFields::ALL,
        }
    }
}

/// Writes canonical CXSMILES with every modeled extension family enabled.
pub fn write_cx_smiles(record: &SmilesRecord) -> Result<String, SmilesParseError> {
    write_cx_smiles_with_params(record, &CxSmilesWriteParams::default())
}

/// Writes CXSMILES from detached values with explicit traversal and field policy.
pub fn write_cx_smiles_with_params(
    record: &SmilesRecord,
    params: &CxSmilesWriteParams,
) -> Result<String, SmilesParseError> {
    record
        .coordinates
        .validate_for_atom_count(record.topology.atoms.len())
        .map_err(|error| SmilesParseError::Model(error.to_string()))?;
    let output = write_smiles_for_cx(record, &params.smiles)?;
    let extension = write_cx_extensions(
        record,
        params.fields,
        &output.atom_order,
        &output.bond_order,
    )?;
    if extension.is_empty() {
        Ok(output.text)
    } else if output.text.is_empty() {
        Ok(extension)
    } else {
        Ok(format!("{} {}", output.text, extension))
    }
}

fn append_extension(addition: String, output: &mut String) {
    // RDKit✔️✔️: void appendToCXExtension(const std::string &addition, std::string &base) {
    // RDKit✔️✔️:   if (!addition.empty()) {
    // RDKit✔️✔️:     if (base.size() > 1) { base += ","; }
    // RDKit✔️✔️:     base += addition;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    if !addition.is_empty() {
        if output.len() > 1 {
            output.push(',');
        }
        output.push_str(&addition);
    }
}

fn write_cx_extensions(
    record: &SmilesRecord,
    fields: CxSmilesFields,
    atom_order: &[AtomId],
    bond_order: &[BondId],
) -> Result<String, SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION SmilesWrite::getCXExtensions field order
    // RDKit✔️✔️: if ((flags & SmilesWrite::CXSmilesFields::CX_COORDS) &&
    // RDKit✔️✔️:     mol.getNumConformers()) {
    // RDKit✔️✔️:   res += "(" + get_coords_block(mol, atomOrder) + ")";
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if ((flags & SmilesWrite::CXSmilesFields::CX_ATOM_LABELS) && needLabels) {
    // RDKit✔️✔️:   auto lbls = get_atomlabel_block(mol, atomOrder);
    // RDKit✔️✔️:   if (!lbls.empty()) {
    // RDKit✔️✔️:     if (res.size() > 1) { res += ","; }
    // RDKit✔️✔️:     res += "$" + lbls + "$";
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if ((flags & SmilesWrite::CXSmilesFields::CX_MOLFILE_VALUES) && needValues) {
    // RDKit✔️✔️:   if (res.size() > 1) { res += ","; }
    // RDKit✔️✔️:   res += "$_AV:" +
    // RDKit✔️✔️:          get_value_block(mol, atomOrder, common_properties::molFileValue) + "$";
    // RDKit✔️✔️: }
    // RDKit✔️✔️: auto radblock = get_radical_block(mol, atomOrder);
    // RDKit✔️✔️: if ((flags & SmilesWrite::CXSmilesFields::CX_RADICALS) && radblock.size()) {
    // RDKit✔️✔️:   if (res.size() > 1) { res += ","; }
    // RDKit✔️✔️:   res += radblock;
    // RDKit✔️✔️:   if (res.back() == ',') { res.erase(res.size() - 1); }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (flags & SmilesWrite::CXSmilesFields::CX_ATOM_PROPS) {
    // RDKit✔️✔️:   const auto atomblock = get_atom_props_block(mol, atomOrder);
    // RDKit✔️✔️:   appendToCXExtension(atomblock, res);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (flags & SmilesWrite::CXSmilesFields::CX_BOND_CFG) {
    // RDKit✔️✔️:   bool includeCoords = flags & SmilesWrite::CX_COORDS && mol.getNumConformers();
    // RDKit✔️✔️:   const auto cfgblock = get_bond_config_block(mol, atomOrder, bondOrder,
    // RDKit✔️✔️:                                               includeCoords, wedgeBonds);
    // RDKit✔️✔️:   appendToCXExtension(cfgblock, res);
    // RDKit✔️✔️:   const auto cistransblock = get_ringbond_cistrans_block(mol, atomOrder, bondOrder);
    // RDKit✔️✔️:   appendToCXExtension(cistransblock, res);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (flags & SmilesWrite::CXSmilesFields::CX_COORDINATE_BONDS) {
    // RDKit✔️✔️:   const auto block = get_coord_or_hydrogen_bonds_block(
    // RDKit✔️✔️:       mol, Bond::BondType::DATIVE, "C", atomOrder, bondOrder);
    // RDKit✔️✔️:   appendToCXExtension(block, res);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (flags & SmilesWrite::CXSmilesFields::CX_HYDROGEN_BONDS) {
    // RDKit✔️✔️:   const auto block = get_coord_or_hydrogen_bonds_block(
    // RDKit✔️✔️:       mol, Bond::BondType::HYDROGEN, "H", atomOrder, bondOrder);
    // RDKit✔️✔️:   appendToCXExtension(block, res);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (flags & SmilesWrite::CXSmilesFields::CX_ZERO_BONDS) {
    // RDKit✔️✔️:   const auto block = get_zerobonds_block(mol, atomOrder, bondOrder);
    // RDKit✔️✔️:   appendToCXExtension(block, res);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (flags & SmilesWrite::CXSmilesFields::CX_LINKNODES) {
    // RDKit✔️✔️:   const auto linknodeblock = get_linknodes_block(mol, atomOrder);
    // RDKit✔️✔️:   appendToCXExtension(linknodeblock, res);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (flags & SmilesWrite::CXSmilesFields::CX_ENHANCEDSTEREO) {
    // RDKit✔️✔️:   const auto stereoblock = get_enhanced_stereo_block(mol, atomOrder, wedgeBonds);
    // RDKit✔️✔️:   appendToCXExtension(stereoblock, res);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (flags & SmilesWrite::CXSmilesFields::CX_SGROUPS) {
    // RDKit✔️✔️:   const auto sgroupdatablock = get_sgroup_data_block(mol, atomOrder);
    // RDKit✔️✔️:   appendToCXExtension(sgroupdatablock, res);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (flags & SmilesWrite::CXSmilesFields::CX_POLYMER) {
    // RDKit✔️✔️:   const auto sgrouppolyblock = get_sgroup_polymer_block(mol, atomOrder, bondOrder);
    // RDKit✔️✔️:   appendToCXExtension(sgrouppolyblock, res);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (flags & (SmilesWrite::CXSmilesFields::CX_SGROUPS |
    // RDKit✔️✔️:              SmilesWrite::CXSmilesFields::CX_POLYMER)) {
    // RDKit✔️✔️:   const auto sgrouphierarchyblock = get_sgroup_hierarchy_block(mol);
    // RDKit✔️✔️:   appendToCXExtension(sgrouphierarchyblock, res);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (res.size() > 1) { res += "|"; } else { res = ""; }
    // END RDKIT CPP FUNCTION SmilesWrite::getCXExtensions field order
    let mut result = String::from("|");

    if fields.contains(CxSmilesFields::COORDS)
        && let Some(coords) = write_coordinates(record, atom_order)
    {
        result.push('(');
        result.push_str(&coords);
        result.push(')');
    }
    if fields.contains(CxSmilesFields::ATOM_LABELS) {
        let labels = write_atom_labels(record, atom_order);
        if !labels.is_empty() {
            append_extension(format!("${labels}$"), &mut result);
        }
    }
    if fields.contains(CxSmilesFields::MOLFILE_VALUES) {
        let values = write_atom_values(record, atom_order);
        if !values.is_empty() {
            append_extension(format!("$_AV:{values}$"), &mut result);
        }
    }
    if fields.contains(CxSmilesFields::RADICALS) {
        append_extension(write_radicals(record, atom_order)?, &mut result);
    }
    if fields.contains(CxSmilesFields::ATOM_PROPS) {
        append_extension(write_atom_properties(record, atom_order), &mut result);
    }

    let ring_info = if fields.contains(CxSmilesFields::BOND_CFG) {
        Some(
            fast_find_rings_from_parts(
                record.topology.atoms.len(),
                &record.topology.bonds,
                &record.topology.adjacency,
            )
            .map_err(|error| SmilesParseError::WriterStereo(error.to_string()))?,
        )
    } else {
        None
    };
    if fields.contains(CxSmilesFields::BOND_CFG) {
        append_extension(
            write_bond_config(
                record,
                atom_order,
                bond_order,
                has_coordinates(record),
                false,
            ),
            &mut result,
        );
        append_extension(
            write_ring_bond_stereo(
                record,
                atom_order,
                bond_order,
                ring_info.as_ref().expect("ring info was computed"),
            ),
            &mut result,
        );
    } else if fields.contains(CxSmilesFields::BOND_ATROPISOMER) {
        append_extension(
            write_bond_config(
                record,
                atom_order,
                bond_order,
                has_coordinates(record),
                true,
            ),
            &mut result,
        );
    }
    if fields.contains(CxSmilesFields::COORDINATE_BONDS) {
        append_extension(
            write_coordinate_bonds(record, atom_order, bond_order, BondOrder::Dative, "C"),
            &mut result,
        );
    }
    if fields.contains(CxSmilesFields::HYDROGEN_BONDS) {
        append_extension(
            write_coordinate_bonds(record, atom_order, bond_order, BondOrder::Hydrogen, "H"),
            &mut result,
        );
    }
    if fields.contains(CxSmilesFields::ZERO_BONDS) {
        append_extension(write_zero_bonds(record, bond_order), &mut result);
    }
    if fields.contains(CxSmilesFields::LINKNODES) {
        append_extension(write_link_nodes(record, atom_order), &mut result);
    }
    if fields.contains(CxSmilesFields::ENHANCED_STEREO) {
        append_extension(write_enhanced_stereo(record, atom_order), &mut result);
    }
    if fields.contains(CxSmilesFields::SGROUPS) {
        append_extension(write_data_sgroups(record, atom_order), &mut result);
    }
    if fields.contains(CxSmilesFields::POLYMER) {
        append_extension(
            write_polymer_sgroups(record, atom_order, bond_order),
            &mut result,
        );
    }
    if fields.contains(CxSmilesFields::SGROUPS) || fields.contains(CxSmilesFields::POLYMER) {
        append_extension(
            write_sgroup_hierarchy(
                record,
                fields.contains(CxSmilesFields::SGROUPS),
                fields.contains(CxSmilesFields::POLYMER),
            ),
            &mut result,
        );
    }

    if result.len() == 1 {
        Ok(String::new())
    } else {
        result.push('|');
        Ok(result)
    }
}

fn atom_positions(atom_order: &[AtomId], atom_count: usize) -> Vec<Option<usize>> {
    let mut positions = vec![None; atom_count];
    for (position, atom) in atom_order.iter().copied().enumerate() {
        positions[atom.index()] = Some(position);
    }
    positions
}

fn bond_positions(bond_order: &[BondId], bond_count: usize) -> Vec<Option<usize>> {
    let mut positions = vec![None; bond_count];
    for (position, bond) in bond_order.iter().copied().enumerate() {
        positions[bond.index()] = Some(position);
    }
    positions
}

fn zero_small(value: f64) -> f64 {
    if value.abs() < 1e-4 { 0.0 } else { value }
}

// `boost::format("%g")` uses six significant digits. This implements the
// same fixed/scientific selection and trimming for finite coordinate values.
fn format_general(value: f64) -> String {
    if value == 0.0 {
        return "0".to_owned();
    }
    if !value.is_finite() {
        return value.to_string();
    }
    let exponent = value.abs().log10().floor() as i32;
    let scale = 10_f64.powi(5 - exponent);
    let rounded = (value * scale).round() / scale;
    let exponent = rounded.abs().log10().floor() as i32;
    if !(-4..6).contains(&exponent) {
        let mut text = format!("{rounded:.5e}");
        let exponent_at = text.find('e').expect("scientific format contains exponent");
        let exponent_text = text.split_off(exponent_at);
        while text.ends_with('0') {
            text.pop();
        }
        if text.ends_with('.') {
            text.pop();
        }
        let exponent_value = exponent_text[1..].parse::<i32>().unwrap_or(exponent);
        format!("{text}e{exponent_value:+03}")
    } else {
        let decimals = usize::try_from((5 - exponent).max(0)).unwrap_or(0);
        let mut text = format!("{rounded:.decimals$}");
        if text.contains('.') {
            while text.ends_with('0') {
                text.pop();
            }
            if text.ends_with('.') {
                text.pop();
            }
        }
        text
    }
}

enum CoordinateSource<'a> {
    ThreeD(&'a Conformer3D),
    TwoD(&'a Conformer2D),
}

fn coordinate_source(record: &SmilesRecord) -> Option<CoordinateSource<'_>> {
    record
        .coordinates
        .conformers_3d
        .first()
        .map(CoordinateSource::ThreeD)
        .or_else(|| {
            record
                .coordinates
                .conformers_2d
                .first()
                .map(CoordinateSource::TwoD)
        })
}

fn has_coordinates(record: &SmilesRecord) -> bool {
    coordinate_source(record).is_some()
}

fn write_coordinates(record: &SmilesRecord, atom_order: &[AtomId]) -> Option<String> {
    // BEGIN RDKIT CPP FUNCTION get_coords_block
    // RDKit✔️✔️: const auto &conf = mol.getConformer();
    // RDKit✔️✔️: for (auto idx : atomOrder) {
    // RDKit✔️✔️:   const auto &pt = conf.getAtomPos(idx);
    // RDKit✔️✔️:   res += boost::str(boost::format("%g,%g,") % zero_small_vals(pt.x) %
    // RDKit✔️✔️:                     zero_small_vals(pt.y));
    // RDKit✔️✔️:   if (conf.is3D()) {
    // RDKit✔️✔️:     auto zc = boost::str(boost::format("%g") % zero_small_vals(pt.z));
    // RDKit✔️✔️:     if (zc != "0") { res += zc; }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION get_coords_block
    let source = coordinate_source(record)?;
    let mut points = Vec::with_capacity(atom_order.len());
    for atom in atom_order {
        let point = match source {
            CoordinateSource::ThreeD(conformer) => {
                let point = conformer.coordinates()[atom.index()];
                let mut text = format!(
                    "{},{},",
                    format_general(zero_small(point[0])),
                    format_general(zero_small(point[1]))
                );
                if conformer.is_3d() {
                    let z = format_general(zero_small(point[2]));
                    if z != "0" {
                        text.push_str(&z);
                    }
                }
                text
            }
            CoordinateSource::TwoD(conformer) => {
                let point = conformer.coordinates()[atom.index()];
                format!(
                    "{},{},",
                    format_general(zero_small(point[0])),
                    format_general(zero_small(point[1]))
                )
            }
        };
        points.push(point);
    }
    Some(points.join(";"))
}

fn write_atom_labels(record: &SmilesRecord, atom_order: &[AtomId]) -> String {
    // BEGIN RDKIT CPP FUNCTION get_atomlabel_block
    // RDKit✔️✔️: if (atom->getPropIfPresent(common_properties::_QueryAtomGenericLabel, lbl)) {
    // RDKit✔️✔️:   res += quote_string(lbl + "_p");
    // RDKit✔️✔️: } else if (!atom->getAtomicNum() &&
    // RDKit✔️✔️:            atom->getPropIfPresent(common_properties::dummyLabel, lbl) &&
    // RDKit✔️✔️:            std::find(SmilesParseOps::pseudoatoms.begin(),
    // RDKit✔️✔️:                      SmilesParseOps::pseudoatoms.end(), lbl) !=
    // RDKit✔️✔️:                SmilesParseOps::pseudoatoms.end()) {
    // RDKit✔️✔️:   res += quote_string(lbl + "_p");
    // RDKit✔️✔️: } else if (!atom->getAtomicNum() &&
    // RDKit✔️✔️:            atom->getPropIfPresent(common_properties::_fromAttachPoint,
    // RDKit✔️✔️:                                       val) &&
    // RDKit✔️✔️:            (val == 1 || val == 2)) {
    // RDKit✔️✔️:   res += quote_string("_AP" + std::to_string(val));
    // RDKit✔️✔️: } else if (atom->getPropIfPresent(common_properties::atomLabel, lbl)) {
    // RDKit✔️✔️:   res += quote_string(lbl);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION get_atomlabel_block
    const PSEUDOATOMS: [&str; 8] = ["Pol", "Mod", "Het", "Any", "A", "Q", "X", "*"];
    let labels = atom_order
        .iter()
        .map(|atom| {
            let atom = &record.topology.atoms[atom.index()];
            if let Some(label) = atom.prop("_QueryAtomGenericLabel") {
                format!("{label}_p")
            } else if atom.atomic_number() == 0
                && atom
                    .prop("dummyLabel")
                    .is_some_and(|label| PSEUDOATOMS.contains(&label))
            {
                format!("{}_p", atom.prop("dummyLabel").unwrap_or_default())
            } else if atom.atomic_number() == 0
                && atom
                    .prop("_fromAttachPoint")
                    .is_some_and(|value| matches!(value, "1" | "2"))
            {
                format!("_AP{}", atom.prop("_fromAttachPoint").unwrap_or_default())
            } else {
                atom.prop("atomLabel").unwrap_or_default().to_owned()
            }
        })
        .collect::<Vec<_>>();
    if labels.iter().all(String::is_empty) {
        String::new()
    } else {
        labels.join(";")
    }
}

fn write_atom_values(record: &SmilesRecord, atom_order: &[AtomId]) -> String {
    let values = atom_order
        .iter()
        .map(|atom| {
            record.topology.atoms[atom.index()]
                .prop("molFileValue")
                .unwrap_or_default()
                .to_owned()
        })
        .collect::<Vec<_>>();
    if values.iter().all(String::is_empty) {
        String::new()
    } else {
        values.join(";")
    }
}

fn write_radicals(
    record: &SmilesRecord,
    atom_order: &[AtomId],
) -> Result<String, SmilesParseError> {
    // BEGIN RDKIT CPP FUNCTION get_radical_block
    // RDKit✔️✔️: std::map<unsigned int, std::vector<unsigned int>> rads;
    // RDKit✔️✔️: for (unsigned int i = 0; i < atomOrder.size(); ++i) {
    // RDKit✔️✔️:   auto nrad = mol.getAtomWithIdx(atomOrder[i])->getNumRadicalElectrons();
    // RDKit✔️✔️:   if (nrad) { rads[nrad].push_back(i); }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: switch (pr.first) { case 1: res += "^1:"; break;
    // RDKit✔️✔️:   case 2: res += "^2:"; break; case 3: res += "^5:"; break; }
    // END RDKIT CPP FUNCTION get_radical_block
    let mut radicals = BTreeMap::<u8, Vec<usize>>::new();
    for (position, atom) in atom_order.iter().copied().enumerate() {
        let count = record.topology.atoms[atom.index()].radical_electrons();
        if count > 0 {
            radicals.entry(count).or_default().push(position);
        }
    }
    let mut blocks = Vec::new();
    for (count, atoms) in radicals {
        let marker = match count {
            1 => "^1:",
            2 => "^2:",
            3 => "^5:",
            _ => {
                return Err(SmilesParseError::UnsupportedWriter(
                    "CXSMILES supports at most three radical electrons per atom",
                ));
            }
        };
        blocks.push(format!(
            "{marker}{}",
            atoms
                .iter()
                .map(usize::to_string)
                .collect::<Vec<_>>()
                .join(",")
        ));
    }
    Ok(blocks.join(","))
}

fn quote_atom_property(text: &str) -> String {
    text.replace('.', "&#46;")
}

fn write_atom_properties(record: &SmilesRecord, atom_order: &[AtomId]) -> String {
    // BEGIN RDKIT CPP FUNCTION get_atom_props_block
    // RDKit✔️✔️: constexpr std::array<std::string_view, 7> skip = {
    // RDKit✔️✔️:     atomLabel, molFileValue, molParity, molAtomMapNumber,
    // RDKit✔️✔️:     molStereoCare, molRxnExactChange, molInversionFlag};
    // RDKit✔️✔️: for (const auto &pn : atom->getPropList(false, false)) {
    // RDKit✔️✔️:   if (std::find(skip.begin(), skip.end(), pn) == skip.end()) {
    // RDKit✔️✔️:     if (pn == "dummyLabel" &&
    // RDKit✔️✔️:         (isAttachmentPoint || pv == "*" ||
    // RDKit✔️✔️:          std::find(SmilesParseOps::pseudoatoms.begin(),
    // RDKit✔️✔️:                    SmilesParseOps::pseudoatoms.end(),
    // RDKit✔️✔️:                    pv) != SmilesParseOps::pseudoatoms.end())) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     res +=
    // RDKit✔️✔️:         boost::str(boost::format(":%d.%s.%s") % which %
    // RDKit✔️✔️:                    quote_atomprop_string(pn) % quote_atomprop_string(pv));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION get_atom_props_block
    const SKIP: [&str; 7] = [
        "atomLabel",
        "molFileValue",
        "molParity",
        "molAtomMapNumber",
        "molStereoCare",
        "molRxnExactChange",
        "molInversionFlag",
    ];
    const PSEUDOATOMS: [&str; 8] = ["Pol", "Mod", "Het", "Any", "A", "Q", "X", "*"];
    let mut entries = Vec::new();
    for (position, atom_id) in atom_order.iter().copied().enumerate() {
        let atom = &record.topology.atoms[atom_id.index()];
        let attachment = atom.atomic_number() == 0 && atom.prop("_fromAttachPoint").is_some();
        for (name, value) in atom.props() {
            if name.starts_with('_') || atom.is_prop_computed(name) || SKIP.contains(&name.as_str())
            {
                continue;
            }
            if name == "dummyLabel"
                && (attachment || value == "*" || PSEUDOATOMS.contains(&value.as_str()))
            {
                continue;
            }
            entries.push(format!(
                "{position}.{}.{}",
                quote_atom_property(name),
                quote_atom_property(value)
            ));
        }
    }
    if entries.is_empty() {
        String::new()
    } else {
        format!("atomProp:{}", entries.join(":"))
    }
}

fn can_have_direction(bond: &Bond) -> bool {
    matches!(bond.order(), BondOrder::Single | BondOrder::Aromatic)
}

fn normalized_wedge_direction(bond: &Bond) -> BondDirection {
    match bond.direction() {
        BondDirection::BeginDash | BondDirection::BeginWedge | BondDirection::Unknown => {
            bond.direction()
        }
        _ => BondDirection::None,
    }
}

fn write_bond_config(
    record: &SmilesRecord,
    atom_order: &[AtomId],
    bond_order: &[BondId],
    coords_included: bool,
    atropisomer_only: bool,
) -> String {
    // BEGIN RDKIT CPP FUNCTION get_bond_config_block core emission
    // RDKit✔️✔️: for (unsigned int i = 0; i < bondOrder.size(); ++i) {
    // RDKit✔️✔️:   const auto bond = mol.getBondWithIdx(bondOrder[i]);
    // RDKit✔️✔️:   unsigned int wedgeStartAtomIdx = bond->getBeginAtomIdx();
    // RDKit✔️✔️:   if (!canHaveDirection(*bond)) { continue; }
    // RDKit✔️✔️:   Bond::BondDir bd = bond->getBondDir();
    // RDKit✔️✔️:   switch (bd) {
    // RDKit✔️✔️:     case Bond::BondDir::BEGINDASH:
    // RDKit✔️✔️:     case Bond::BondDir::BEGINWEDGE:
    // RDKit✔️✔️:     case Bond::BondDir::UNKNOWN:
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     default:
    // RDKit✔️✔️:       bd = Bond::BondDir::NONE;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!atropisomerOnly && bd == Bond::BondDir::NONE &&
    // RDKit✔️✔️:       bond->getPropIfPresent(common_properties::_MolFileBondCfg, cfg)) {
    // RDKit✔️✔️:     switch (cfg) {
    // RDKit✔️✔️:       case 1:
    // RDKit✔️✔️:         bd = Bond::BondDir::BEGINWEDGE;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       case 2:
    // RDKit✔️✔️:         bd = Bond::BondDir::UNKNOWN;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       case 3:
    // RDKit✔️✔️:         bd = Bond::BondDir::BEGINDASH;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       default:
    // RDKit✔️✔️:         bd = Bond::BondDir::NONE;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (bd == Bond::BondDir::UNKNOWN) {
    // RDKit✔️✔️:     wType = "w";
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   else if (coordsIncluded || isAnAtropisomer) {
    // RDKit✔️✔️:     if (bd == Bond::BondDir::BEGINWEDGE) {
    // RDKit✔️✔️:       wType = "wU";
    // RDKit✔️✔️:     } else if (bd == Bond::BondDir::BEGINDASH) {
    // RDKit✔️✔️:       wType = "wD";
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   wParts[wType].push_back(format("%d.%d", begAtomOrder, i));
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION get_bond_config_block core emission
    let positions = atom_positions(atom_order, record.topology.atoms.len());
    let mut parts = BTreeMap::<&'static str, Vec<String>>::new();
    for (bond_position, bond_id) in bond_order.iter().copied().enumerate() {
        let bond = &record.topology.bonds[bond_id.index()];
        if !can_have_direction(bond) {
            continue;
        }
        let mut direction = normalized_wedge_direction(bond);
        if matches!(
            direction,
            BondDirection::BeginWedge | BondDirection::BeginDash
        ) && !matches!(
            record.topology.atoms[bond.begin().index()].chiral_tag(),
            cosmolkit_types::ChiralTag::TetrahedralCw | cosmolkit_types::ChiralTag::TetrahedralCcw
        ) {
            // `MolToSmiles()` cleans stereochemistry before
            // `get_bond_config_block()`, which removes non-stereogenic
            // wedge/dash directions such as `CC |wU:0.0|`.
            direction = BondDirection::None;
        }
        if !atropisomer_only && direction == BondDirection::None {
            direction = match bond
                .prop("_MolFileBondCfg")
                .and_then(|value| value.parse::<u8>().ok())
            {
                Some(1) => BondDirection::BeginWedge,
                Some(2) => BondDirection::Unknown,
                Some(3) => BondDirection::BeginDash,
                _ => BondDirection::None,
            };
        }
        if matches!(
            direction,
            BondDirection::BeginWedge | BondDirection::BeginDash
        ) && !matches!(
            record.topology.atoms[bond.begin().index()].chiral_tag(),
            cosmolkit_types::ChiralTag::TetrahedralCw | cosmolkit_types::ChiralTag::TetrahedralCcw
        ) {
            direction = BondDirection::None;
        }
        // Atropisomeric bond stereochemistry is rejected by the base detached
        // writer, so there is no modeled atropisomer-derived wedge here.
        if atropisomer_only {
            continue;
        }
        let kind = match direction {
            BondDirection::Unknown => Some("w"),
            BondDirection::BeginWedge if coords_included => Some("wU"),
            BondDirection::BeginDash if coords_included => Some("wD"),
            _ => None,
        };
        let (Some(kind), Some(begin_position)) = (kind, positions[bond.begin().index()]) else {
            continue;
        };
        parts
            .entry(kind)
            .or_default()
            .push(format!("{begin_position}.{bond_position}"));
    }
    parts
        .into_iter()
        .map(|(kind, entries)| format!("{kind}:{}", entries.join(",")))
        .collect::<Vec<_>>()
        .join(",")
}

fn write_coordinate_bonds(
    record: &SmilesRecord,
    atom_order: &[AtomId],
    bond_order: &[BondId],
    order: BondOrder,
    symbol: &str,
) -> String {
    // BEGIN RDKIT CPP FUNCTION get_coord_or_hydrogen_bonds_block
    // RDKit✔️✔️: for (unsigned int i = 0; i < bondOrder.size(); ++i) {
    // RDKit✔️✔️:   const auto bond = mol.getBondWithIdx(bondOrder[i]);
    // RDKit✔️✔️:   if (bond->getBondType() != bondType) { continue; }
    // RDKit✔️✔️:   auto begAtomOrder = std::find(atomOrder.begin(), atomOrder.end(),
    // RDKit✔️✔️:                                  bond->getBeginAtomIdx()) - atomOrder.begin();
    // RDKit✔️✔️:   res += format("%d.%d", begAtomOrder, i);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION get_coord_or_hydrogen_bonds_block
    let positions = atom_positions(atom_order, record.topology.atoms.len());
    let entries = bond_order
        .iter()
        .copied()
        .enumerate()
        .filter_map(|(position, bond_id)| {
            let bond = &record.topology.bonds[bond_id.index()];
            let matches = if order == BondOrder::Dative {
                matches!(bond.order(), BondOrder::Dative | BondOrder::DativeOne)
            } else {
                bond.order() == order
            };
            matches.then(|| {
                positions[bond.begin().index()].map(|begin| format!("{begin}.{position}"))
            })?
        })
        .collect::<Vec<_>>();
    if entries.is_empty() {
        String::new()
    } else {
        format!("{symbol}:{}", entries.join(","))
    }
}

fn write_zero_bonds(record: &SmilesRecord, bond_order: &[BondId]) -> String {
    // RDKit✔️✔️: if (bond->getBondType() != Bond::BondType::ZERO) { continue; }
    // RDKit✔️✔️: res += boost::str(boost::format("%d") % i);
    let entries = bond_order
        .iter()
        .copied()
        .enumerate()
        .filter_map(|(position, bond)| {
            (record.topology.bonds[bond.index()].order() == BondOrder::Zero)
                .then(|| position.to_string())
        })
        .collect::<Vec<_>>();
    if entries.is_empty() {
        String::new()
    } else {
        format!("Z:{}", entries.join(","))
    }
}

fn other_atom(bond: &Bond, atom: AtomId) -> Option<AtomId> {
    if bond.begin() == atom {
        Some(bond.end())
    } else if bond.end() == atom {
        Some(bond.begin())
    } else {
        None
    }
}

fn write_ring_bond_stereo(
    record: &SmilesRecord,
    atom_order: &[AtomId],
    bond_order: &[BondId],
    rings: &RingInfo,
) -> String {
    // BEGIN RDKIT CPP FUNCTION get_ringbond_cistrans_block
    // RDKit✔️✔️: if (!rinfo->numBondRings(idx) ||
    // RDKit✔️✔️:     rinfo->minBondRingSize(idx) < Chirality::minRingSizeForDoubleBondStereo) continue;
    // RDKit✔️✔️: if (bond->getBondType() != Bond::BondType::DOUBLE &&
    // RDKit✔️✔️:     bond->getBondType() != Bond::BondType::AROMATIC) {
    // RDKit✔️✔️:   continue;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (bstereo == Bond::BondStereo::STEREOANY) {
    // RDKit✔️✔️:   if (ctu.empty()) { ctu += "ctu:"; } else { ctu += ","; }
    // RDKit✔️✔️:   ctu += label;
    // RDKit✔️✔️: } else {
    // RDKit✔️✔️:   bool needSwap = false;
    // RDKit✔️✔️:   if (bstereo == Bond::BondStereo::STEREOCIS || needSwap) {
    // RDKit✔️✔️:     if (c.empty()) { c += "c:"; } else { c += ","; }
    // RDKit✔️✔️:     c += label;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     if (t.empty()) { t += "t:"; } else { t += ","; }
    // RDKit✔️✔️:     t += label;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: return c + t + ctu;
    // END RDKIT CPP FUNCTION get_ringbond_cistrans_block
    const MIN_RING_SIZE: usize = 8;
    let positions = atom_positions(atom_order, record.topology.atoms.len());
    let mut cis = Vec::new();
    let mut trans = Vec::new();
    let mut unknown = Vec::new();
    for (bond_position, bond_id) in bond_order.iter().copied().enumerate() {
        if rings.num_bond_rings(bond_id) == 0 || rings.min_bond_ring_size(bond_id) < MIN_RING_SIZE {
            continue;
        }
        let bond = &record.topology.bonds[bond_id.index()];
        if !matches!(bond.order(), BondOrder::Double | BondOrder::Aromatic)
            || !matches!(
                bond.stereo(),
                BondStereo::Any | BondStereo::Cis | BondStereo::Trans
            )
        {
            continue;
        }
        if bond.stereo() == BondStereo::Any {
            unknown.push(bond_position.to_string());
            continue;
        }
        let Some([begin_ref, end_ref]) = bond.stereo_atoms() else {
            continue;
        };
        let mut swap = false;
        for (center, opposite, reference) in [
            (bond.begin(), bond.end(), begin_ref),
            (bond.end(), bond.begin(), end_ref),
        ] {
            let neighbors = record.topology.adjacency.neighbors_of(center.index());
            if neighbors.len() <= 2 {
                continue;
            }
            let Some(reference_position) = positions[reference.index()] else {
                continue;
            };
            for neighbor in neighbors {
                let neighbor_atom =
                    other_atom(&record.topology.bonds[neighbor.bond.index()], center)
                        .expect("adjacency endpoint");
                if neighbor_atom != opposite
                    && neighbor_atom != reference
                    && positions[neighbor_atom.index()]
                        .is_some_and(|position| position < reference_position)
                {
                    swap = !swap;
                }
            }
        }
        if bond.stereo() == BondStereo::Cis || swap {
            cis.push(bond_position.to_string());
        } else {
            trans.push(bond_position.to_string());
        }
    }
    [
        (!cis.is_empty()).then(|| format!("c:{}", cis.join(","))),
        (!trans.is_empty()).then(|| format!("t:{}", trans.join(","))),
        (!unknown.is_empty()).then(|| format!("ctu:{}", unknown.join(","))),
    ]
    .into_iter()
    .flatten()
    .collect::<Vec<_>>()
    .join(",")
}

fn stereo_kind_order(kind: StereoGroupKind) -> u8 {
    match kind {
        StereoGroupKind::Absolute => 0,
        StereoGroupKind::Or => 1,
        StereoGroupKind::And => 2,
    }
}

fn assign_stereo_group_ids(groups: &[StereoGroup]) -> Vec<Option<u32>> {
    // BEGIN RDKIT CPP FUNCTION assignStereoGroupIds helpers
    // RDKit✔️✔️: if (groupId == 0) { return; }
    // RDKit✔️✔️: if (ids[groupId]) { sg.setWriteId(0); } else { ids[groupId] = true; }
    // RDKit✔️✔️: if (sg.getWriteId() == 0) {
    // RDKit✔️✔️:   ++nextId; while (nextId < ids.size() && ids[nextId]) { ++nextId; }
    // RDKit✔️✔️:   sg.setWriteId(nextId);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION assignStereoGroupIds helpers
    // The detached model's ID is the source/read ID. RDKit does not forward
    // read IDs into write IDs during CX parsing, so every parsed OR/AND group
    // enters `assignStereoGroupIds()` with write ID zero.
    let mut ids = vec![None; groups.len()];
    let mut next_or = 0;
    let mut next_and = 0;
    for (index, group) in groups.iter().enumerate() {
        if group.kind() == StereoGroupKind::Absolute {
            continue;
        }
        let next = match group.kind() {
            StereoGroupKind::Or => &mut next_or,
            StereoGroupKind::And => &mut next_and,
            StereoGroupKind::Absolute => unreachable!(),
        };
        *next += 1;
        ids[index] = Some(*next);
    }
    ids
}

fn write_enhanced_stereo(record: &SmilesRecord, atom_order: &[AtomId]) -> String {
    // BEGIN RDKIT CPP FUNCTION getSortedStereoGroupsAndIndices/get_enhanced_stereo_block
    // RDKit✔️✔️: const auto newAtomIndexes = getSortedMappedIndexes(atomIds, revOrder);
    // RDKit✔️✔️: if (!newAtomIndexes.empty()) sortingGroups.emplace_back(sg, newAtomIndexes);
    // RDKit✔️✔️: // sort by 1) StereoGroup type; 2) StereoGroup id; 3) atom indexes
    // RDKit✔️✔️: assignStereoGroupIds(groups);
    // RDKit✔️✔️: switch (sgItr->getGroupType()) {
    // RDKit✔️✔️:   case StereoGroupType::STEREO_ABSOLUTE:
    // RDKit✔️✔️:     res << "a:";
    // RDKit✔️✔️:     break;
    // RDKit✔️✔️:   case StereoGroupType::STEREO_OR:
    // RDKit✔️✔️:     res << "o" << sgItr->getWriteId() << ":";
    // RDKit✔️✔️:     break;
    // RDKit✔️✔️:   case StereoGroupType::STEREO_AND:
    // RDKit✔️✔️:     res << "&" << sgItr->getWriteId() << ":";
    // RDKit✔️✔️:     break;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: for (const auto &aid : *grpAtomsItr) {
    // RDKit✔️✔️:   res << aid << ",";
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION getSortedStereoGroupsAndIndices/get_enhanced_stereo_block
    let positions = atom_positions(atom_order, record.topology.atoms.len());
    let mut groups = record
        .topology
        .stereo_groups
        .iter()
        .filter_map(|group| {
            let mut atoms = group
                .atoms()
                .iter()
                .filter_map(|atom| positions[atom.index()])
                .collect::<Vec<_>>();
            atoms.sort_unstable();
            (!atoms.is_empty()).then(|| (group.clone(), atoms))
        })
        .collect::<Vec<_>>();
    groups.sort_by(|(left_group, left_atoms), (right_group, right_atoms)| {
        stereo_kind_order(left_group.kind())
            .cmp(&stereo_kind_order(right_group.kind()))
            .then_with(|| left_atoms.cmp(right_atoms))
    });
    let sorted_groups = groups
        .iter()
        .map(|(group, _)| group.clone())
        .collect::<Vec<_>>();
    let ids = assign_stereo_group_ids(&sorted_groups);
    groups
        .into_iter()
        .zip(ids)
        .map(|((group, atoms), id)| {
            let prefix = match group.kind() {
                StereoGroupKind::Absolute => "a".to_owned(),
                StereoGroupKind::Or => format!("o{}", id.unwrap_or(1)),
                StereoGroupKind::And => format!("&{}", id.unwrap_or(1)),
            };
            format!(
                "{prefix}:{}",
                atoms
                    .iter()
                    .map(usize::to_string)
                    .collect::<Vec<_>>()
                    .join(",")
            )
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn write_link_nodes(record: &SmilesRecord, atom_order: &[AtomId]) -> String {
    // BEGIN RDKIT CPP FUNCTION get_linknodes_block
    // RDKit✔️✔️: auto linkNodes = MolEnumerator::utils::getMolLinkNodes(mol, false);
    // RDKit✔️✔️: for (const auto &ln : linkNodes) {
    // RDKit✔️✔️:   unsigned int atomIdx = atomOrder[ln.bondAtoms[0].first];
    // RDKit✔️✔️:   res << atomIdx << ":" << ln.minRep << "." << ln.maxRep;
    // RDKit✔️✔️:   if (mol.getAtomWithIdx(ln.bondAtoms[0].first)->getDegree() > 2) {
    // RDKit✔️✔️:     res << "." << atomOrder[ln.bondAtoms[0].second] << "."
    // RDKit✔️✔️:         << atomOrder[ln.bondAtoms[1].second];
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION get_linknodes_block
    let Some(raw) = record.properties.prop("_MolFileLinkNodes") else {
        return String::new();
    };
    let mut entries = Vec::new();
    for item in raw.split('|').filter(|item| !item.trim().is_empty()) {
        let values = item
            .split_whitespace()
            .map(str::parse::<usize>)
            .collect::<Result<Vec<_>, _>>();
        let Ok(values) = values else { continue };
        if values.len() < 7 || values[2] < 2 {
            continue;
        }
        let Some(center) = values[3].checked_sub(1) else {
            continue;
        };
        let Some(center_output) = atom_order.get(center).map(|atom| atom.index()) else {
            continue;
        };
        let mut entry = format!("{center_output}:{}.{}", values[0], values[1]);
        if record.topology.adjacency.neighbors_of(center).len() > 2 {
            let (Some(first), Some(second)) = (values[4].checked_sub(1), values[6].checked_sub(1))
            else {
                continue;
            };
            let (Some(first_output), Some(second_output)) = (
                atom_order.get(first).map(|atom| atom.index()),
                atom_order.get(second).map(|atom| atom.index()),
            ) else {
                continue;
            };
            entry.push_str(&format!(".{first_output}.{second_output}"));
        }
        entries.push(entry);
    }
    if entries.is_empty() {
        String::new()
    } else {
        format!("LN:{}", entries.join(","))
    }
}

fn is_data_sgroup(group: &SubstanceGroup) -> bool {
    matches!(group.kind(), SubstanceGroupKind::Data)
        || group
            .props()
            .get("TYPE")
            .is_some_and(|value| value == "DAT")
}

fn data_sgroup_value(group: &SubstanceGroup, key: &str) -> String {
    group.props().get(key).cloned().unwrap_or_default()
}

fn write_data_sgroups(record: &SmilesRecord, atom_order: &[AtomId]) -> String {
    // BEGIN RDKIT CPP FUNCTION get_sgroup_data_block
    // RDKit✔️✔️: for (const auto &sg : sgs) {
    // RDKit✔️✔️:   if (sg.hasProp("TYPE") && sg.getProp<std::string>("TYPE") == "DAT") {
    // RDKit✔️✔️:     res << "SgD:";
    // RDKit✔️✔️:     for (const auto oaid : sg.getAtoms()) { res << revOrder[oaid] << ","; }
    // RDKit✔️✔️:     res << ":";
    // RDKit✔️✔️:     std::string prop;
    // RDKit✔️✔️:     if (sg.getPropIfPresent("FIELDNAME", prop) && !prop.empty()) {
    // RDKit✔️✔️:       res << prop;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     res << ":";
    // RDKit✔️✔️:     std::vector<std::string> vprop;
    // RDKit✔️✔️:     if (sg.getPropIfPresent("DATAFIELDS", vprop) && !vprop.empty()) {
    // RDKit✔️✔️:       for (const auto &pv : vprop) { res << pv << ","; }
    // RDKit✔️✔️:       res.seekp(-1, res.cur);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     res << ":";
    // RDKit✔️✔️:     if (sg.getPropIfPresent("QUERYOP", prop) && !prop.empty()) { res << prop; }
    // RDKit✔️✔️:     res << ":";
    // RDKit✔️✔️:     if (sg.getPropIfPresent("FIELDINFO", prop) && !prop.empty()) { res << prop; }
    // RDKit✔️✔️:     res << ":";
    // RDKit✔️✔️:     if (sg.getPropIfPresent("FIELDTAG", prop) && !prop.empty()) { res << prop; }
    // RDKit✔️✔️:     res << ":";
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION get_sgroup_data_block
    let positions = atom_positions(atom_order, record.topology.atoms.len());
    record
        .topology
        .substance_groups
        .iter()
        .filter(|group| is_data_sgroup(group))
        .filter_map(|group| {
            let atoms = group
                .atoms()
                .iter()
                .filter_map(|atom| positions[atom.index()])
                .map(|position| position.to_string())
                .collect::<Vec<_>>();
            if atoms.is_empty() {
                return None;
            }
            let data = if group.data_fields().is_empty() {
                data_sgroup_value(group, "DATAFIELDS")
            } else {
                group.data_fields().join(",")
            };
            Some(format!(
                "SgD:{}:{}:{}:{}:{}:{}:",
                atoms.join(","),
                data_sgroup_value(group, "FIELDNAME"),
                data,
                data_sgroup_value(group, "QUERYOP"),
                data_sgroup_value(group, "FIELDINFO"),
                data_sgroup_value(group, "FIELDTAG")
            ))
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn polymer_type(group: &SubstanceGroup) -> Option<&'static str> {
    match group.kind() {
        SubstanceGroupKind::StructuralRepeatUnit => Some("n"),
        SubstanceGroupKind::Monomer => Some("mon"),
        SubstanceGroupKind::Mer => Some("mer"),
        SubstanceGroupKind::Copolymer => match group
            .subtype()
            .or_else(|| group.props().get("SUBTYPE").map(String::as_str))
        {
            Some("ALT") => Some("alt"),
            Some("RAN") => Some("ran"),
            Some("BLO") => Some("blk"),
            _ => Some("co"),
        },
        SubstanceGroupKind::Crosslink => Some("xl"),
        SubstanceGroupKind::Modification => Some("mod"),
        SubstanceGroupKind::MixtureComponent => Some("mix"),
        SubstanceGroupKind::Formulation => Some("f"),
        SubstanceGroupKind::AnyPolymer => Some("any"),
        SubstanceGroupKind::Graft => Some("grf"),
        SubstanceGroupKind::Generic(value) if value == "GEN" => Some("gen"),
        SubstanceGroupKind::Generic(value) if value == "COM" => Some("c"),
        _ => None,
    }
}

fn parse_crossings(group: &SubstanceGroup, primary: &str, fallback: &str) -> Vec<usize> {
    group
        .props()
        .get(primary)
        .or_else(|| group.props().get(fallback))
        .map(|value| {
            value
                .split(',')
                .filter_map(|part| part.trim().parse::<usize>().ok())
                .collect()
        })
        .unwrap_or_default()
}

fn connection_text(group: &SubstanceGroup) -> String {
    group
        .connection()
        .map(|connection| match connection {
            SGroupConnection::HeadToHead => "hh".to_owned(),
            SGroupConnection::HeadToTail => "ht".to_owned(),
            SGroupConnection::Either => "eu".to_owned(),
            SGroupConnection::Unknown(value) => value.to_ascii_lowercase(),
        })
        .or_else(|| {
            group
                .props()
                .get("CONNECT")
                .map(|value| value.to_ascii_lowercase())
        })
        .unwrap_or_default()
}

fn write_polymer_sgroups(
    record: &SmilesRecord,
    atom_order: &[AtomId],
    bond_order: &[BondId],
) -> String {
    // BEGIN RDKIT CPP FUNCTION get_sgroup_polymer_block
    // RDKit✔️✔️: for (const auto &sg : sgs) {
    // RDKit✔️✔️:   std::string typ;
    // RDKit✔️✔️:   if (sg.getPropIfPresent("TYPE", typ) &&
    // RDKit✔️✔️:       reverseTypemap.find(typ) != reverseTypemap.end()) {
    // RDKit✔️✔️:     res << "Sg:" << reverse type << ":";
    // RDKit✔️✔️:     for (const auto oaid : sg.getAtoms()) { res << revAtomOrder[oaid] << ","; }
    // RDKit✔️✔️:     res << ":" << LABEL << ":" << lowercase CONNECT << ":";
    // RDKit✔️✔️:     std::vector<unsigned int> headCrossings;
    // RDKit✔️✔️:     if (sg.getPropIfPresent("XBHEAD", headCrossings) &&
    // RDKit✔️✔️:         headCrossings.size() > 1) {
    // RDKit✔️✔️:       for (auto v : headCrossings) { res << bondOrder[v] << ","; }
    // RDKit✔️✔️:       res.seekp(-1, res.cur);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     res << ":";
    // RDKit✔️✔️:     std::vector<unsigned int> tailCrossings;
    // RDKit✔️✔️:     if (sg.getPropIfPresent("XBCORR", tailCrossings) &&
    // RDKit✔️✔️:         tailCrossings.size() > 2) {
    // RDKit✔️✔️:       for (unsigned int i = 1; i < tailCrossings.size(); i += 2) {
    // RDKit✔️✔️:         res << bondOrder[tailCrossings[i]] << ",";
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       res.seekp(-1, res.cur);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     res << ":";
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION get_sgroup_polymer_block
    let atom_positions = atom_positions(atom_order, record.topology.atoms.len());
    let bond_positions = bond_positions(bond_order, record.topology.bonds.len());
    record
        .topology
        .substance_groups
        .iter()
        .filter_map(|group| {
            let kind = polymer_type(group)?;
            let atoms = group
                .atoms()
                .iter()
                .filter_map(|atom| atom_positions[atom.index()])
                .map(|position| position.to_string())
                .collect::<Vec<_>>();
            if atoms.is_empty() {
                return None;
            }
            let heads = parse_crossings(group, "_headCrossings", "XBHEAD");
            let tails = parse_crossings(group, "_tailCrossings", "XBCORR");
            let head_text = (heads.len() > 1)
                .then(|| {
                    heads
                        .iter()
                        .filter_map(|bond| bond_positions.get(*bond).copied().flatten())
                        .map(|position| position.to_string())
                        .collect::<Vec<_>>()
                        .join(",")
                })
                .unwrap_or_default();
            let tail_text = (tails.len() > 2)
                .then(|| {
                    tails
                        .iter()
                        .skip(1)
                        .step_by(2)
                        .filter_map(|bond| bond_positions.get(*bond).copied().flatten())
                        .map(|position| position.to_string())
                        .collect::<Vec<_>>()
                        .join(",")
                })
                .unwrap_or_default();
            Some(format!(
                "Sg:{kind}:{}:{}:{}:{head_text}:{tail_text}:",
                atoms.join(","),
                group
                    .label()
                    .or_else(|| group.props().get("LABEL").map(String::as_str))
                    .unwrap_or_default(),
                connection_text(group)
            ))
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn sgroup_index(group: &SubstanceGroup) -> usize {
    group
        .props()
        .get("index")
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or_else(|| group.id().index())
}

fn write_sgroup_hierarchy(
    record: &SmilesRecord,
    include_data: bool,
    include_polymer: bool,
) -> String {
    // BEGIN RDKIT CPP FUNCTION get_sgroup_hierarchy_block
    // RDKit✔️✔️: std::string get_sgroup_hierarchy_block(const ROMol &mol) {
    // RDKit✔️✔️:   const auto &sgs = getSubstanceGroups(mol);
    // RDKit✔️✔️:   if (sgs.empty()) {
    // RDKit✔️✔️:     return "";
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::stringstream res;
    // RDKit✔️✔️:   // we need a map from sgroup index to output index;
    // RDKit✔️✔️:   std::map<unsigned int, unsigned int> sgroupOrder;
    // RDKit✔️✔️:   bool parentPresent = false;
    // RDKit✔️✔️:   for (const auto &sg : sgs) {
    // RDKit✔️✔️:     if (sg.hasProp("_cxsmilesOutputIndex")) {
    // RDKit✔️✔️:       unsigned int sgidx = sg.getIndexInMol();
    // RDKit✔️✔️:       sg.getPropIfPresent("index", sgidx);
    // RDKit✔️✔️:       sgroupOrder[sgidx] = sg.getProp<unsigned int>("_cxsmilesOutputIndex");
    // RDKit✔️✔️:       sg.clearProp("_cxsmilesOutputIndex");
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (sg.hasProp("PARENT")) {
    // RDKit✔️✔️:       parentPresent = true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (parentPresent) {
    // RDKit✔️✔️:     // now loop over them and add the information
    // RDKit✔️✔️:     std::map<unsigned int, std::vector<unsigned int>> accum;
    // RDKit✔️✔️:     for (const auto &sg : sgs) {
    // RDKit✔️✔️:       unsigned pidx;
    // RDKit✔️✔️:       if (sg.getPropIfPresent("PARENT", pidx) &&
    // RDKit✔️✔️:           sgroupOrder.find(pidx) != sgroupOrder.end()) {
    // RDKit✔️✔️:         unsigned int sgidx = sg.getIndexInMol();
    // RDKit✔️✔️:         sg.getPropIfPresent("index", sgidx);
    // RDKit✔️✔️:         if (sgroupOrder.find(sgidx) != sgroupOrder.end()) {
    // RDKit✔️✔️:           accum[sgroupOrder[pidx]].push_back(sgroupOrder[sgidx]);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (!accum.empty()) {
    // RDKit✔️✔️:       res << "SgH:";
    // RDKit✔️✔️:       for (const auto &pr : accum) {
    // RDKit✔️✔️:         res << pr.first << ":";
    // RDKit✔️✔️:         for (auto v : pr.second) {
    // RDKit✔️✔️:           res << v << ".";
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         // remove the extra ".":
    // RDKit✔️✔️:         res.seekp(-1, res.cur);
    // RDKit✔️✔️:         res << ",";
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     std::string resStr = res.str();
    // RDKit✔️✔️:     while (!resStr.empty() && resStr.back() == ',') {
    // RDKit✔️✔️:       resStr.pop_back();
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return resStr;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return "";
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION get_sgroup_hierarchy_block
    let mut output_indices = BTreeMap::new();
    let mut next = 0;
    if include_data {
        for group in &record.topology.substance_groups {
            if is_data_sgroup(group) {
                output_indices.insert(sgroup_index(group), next);
                next += 1;
            }
        }
    }
    if include_polymer {
        for group in &record.topology.substance_groups {
            if polymer_type(group).is_some() {
                output_indices.insert(sgroup_index(group), next);
                next += 1;
            }
        }
    }
    // Resolve typed IDs once, not by scanning all rows for every child. Both
    // this map and RDKit's source-index/output-index map cost O(n log n) to
    // build and O(log n) per lookup, with O(n) temporary storage overall.
    let source_indices: BTreeMap<_, _> = record
        .topology
        .substance_groups
        .iter()
        .map(|group| (group.id(), sgroup_index(group)))
        .collect();
    let mut hierarchy = BTreeMap::<usize, Vec<usize>>::new();
    for group in &record.topology.substance_groups {
        let Some(child) = output_indices.get(&sgroup_index(group)).copied() else {
            continue;
        };
        // RDKit resolves the parent through the `PARENT` property, which
        // holds the parent's source `index`, and looks it up in the map keyed
        // by that same `index` space. The typed `parent` stores the parent's
        // `SubstanceGroupId`, so translate it through the parent row's own
        // `index` rather than using the row id directly. Fall back to the
        // preserved `PARENT` property when only that is available.
        let parent_key = group
            .parent()
            .and_then(|parent_id| source_indices.get(&parent_id).copied())
            .or_else(|| {
                group
                    .props()
                    .get("PARENT")
                    .and_then(|value| value.parse::<usize>().ok())
            });
        let Some(parent) = parent_key.and_then(|parent| output_indices.get(&parent).copied())
        else {
            continue;
        };
        hierarchy.entry(parent).or_default().push(child);
    }
    if hierarchy.is_empty() {
        String::new()
    } else {
        format!(
            "SgH:{}",
            hierarchy
                .into_iter()
                .map(|(parent, children)| format!(
                    "{parent}:{}",
                    children
                        .iter()
                        .map(usize::to_string)
                        .collect::<Vec<_>>()
                        .join(".")
                ))
                .collect::<Vec<_>>()
                .join(",")
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{SmilesParseParams, parse_smiles};

    fn parse(input: &str) -> SmilesRecord {
        parse_smiles(input, &SmilesParseParams::default()).expect("parse CXSMILES")
    }

    fn write_noncanonical(input: &str) -> String {
        write_cx_smiles_with_params(
            &parse(input),
            &CxSmilesWriteParams {
                smiles: SmilesWriteParams { canonical: false },
                ..Default::default()
            },
        )
        .expect("write CXSMILES")
    }

    #[test]
    fn writes_coordinates_labels_values_radicals_and_atom_properties_in_source_order() {
        let input = "*C |$foo;$,$_AV:x;y$,atomProp:0.a&#46;b.c&#46;d,^1:1,(1,2,3;3,4,0)|";
        assert_eq!(
            write_noncanonical(input),
            "*[CH2] |(1,2,3;3,4,),$foo;$,$_AV:x;y$,^1:1,atomProp:0.a&#46;b.c&#46;d|"
        );
    }

    #[test]
    fn canonical_output_indices_follow_final_atom_and_bond_order() {
        for (input, expected) in [
            ("OC |$oxygen;carbon$|", "CO |$carbon;oxygen$|"),
            ("OC |$_AV:o;c$|", "CO |$_AV:c;o$|"),
            ("OC |atomProp:0.x.y|", "CO |atomProp:1.x.y|"),
            ("OC |^1:0|", "C[O] |^1:1|"),
            ("OC |w:0.0|", "CO |w:1.0|"),
            ("CC |C:0.0|", "C[CH4] |C:1.0|"),
            ("CC |H:0.0|", "CC |H:0.0|"),
            ("CC |Z:0|", "C~C |Z:0|"),
        ] {
            assert_eq!(write_cx_smiles(&parse(input)).unwrap(), expected, "{input}");
        }
    }

    #[test]
    fn writes_enhanced_stereo_with_source_sorting_and_fresh_write_ids() {
        let input = "F[C@H](Cl)Br.O[C@@H](N)I |&7:1,o2:5|";
        assert_eq!(
            write_cx_smiles(&parse(input)).unwrap(),
            "F[C@H](Cl)Br.N[C@H](O)I |o1:5,&1:1|"
        );
        assert_eq!(
            write_noncanonical(input),
            "F[C@H](Cl)Br.O[C@@H](N)I |o1:5,&1:1|"
        );
    }

    #[test]
    fn cx_wedges_with_2d_coordinates_assign_and_write_tetrahedral_stereo() {
        for (input, expected) in [
            (
                "CC(O)Cl |(-3.9163,5.4767,;-3.9163,3.9367,;-2.5826,3.1667,;-5.25,3.1667,),wU:1.0|",
                "C[C@@H](O)Cl |(-3.9163,5.4767,;-3.9163,3.9367,;-2.5826,3.1667,;-5.25,3.1667,),wU:1.0|",
            ),
            (
                "CC(O)Cl |(-3.9163,5.4767,;-3.9163,3.9367,;-2.5826,3.1667,;-5.25,3.1667,),wD:1.0|",
                "C[C@H](O)Cl |(-3.9163,5.4767,;-3.9163,3.9367,;-2.5826,3.1667,;-5.25,3.1667,),wD:1.0|",
            ),
        ] {
            assert_eq!(write_cx_smiles(&parse(input)).unwrap(), expected, "{input}");
        }
        assert_eq!(
            write_cx_smiles(&parse("CC |(0,0,;1,0,),wU:0.0|")).unwrap(),
            "CC |(0,0,;1,0,)|"
        );
    }

    #[test]
    fn writes_data_polymer_and_hierarchy_records() {
        for (input, expected) in [
            ("OCC |SgD:0,2:foo:bar::::|", "CCO |SgD:2,0:foo:bar::::|"),
            ("OCC |Sg:n:0,1:lab:ht:::|", "CCO |Sg:n:2,1:lab:ht:::|"),
            (
                "OCC |SgD:0:a:b::::,Sg:n:1,2::ht:::,SgH:1:0|",
                "CCO |SgD:2:a:b::::,Sg:n:1,0::ht:::,SgH:1:0|",
            ),
        ] {
            assert_eq!(write_cx_smiles(&parse(input)).unwrap(), expected, "{input}");
        }
    }

    #[test]
    fn review_hierarchy_resolves_sparse_source_indices_and_selected_outputs() {
        use cosmolkit_model::SubstanceGroupId;
        let mut record = parse("C");
        let mut child = SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data);
        child.set_prop("index", "42");
        child.set_parent(SubstanceGroupId::new(2));
        let mut sibling = SubstanceGroup::new(SubstanceGroupId::new(1), SubstanceGroupKind::Data);
        sibling.set_prop("index", "7");
        // Preserve the existing raw-property path when no typed parent exists.
        sibling.set_prop("PARENT", "900");
        let mut parent = SubstanceGroup::new(
            SubstanceGroupId::new(2),
            SubstanceGroupKind::StructuralRepeatUnit,
        );
        parent.set_prop("index", "900");
        record.topology.substance_groups = vec![child, sibling, parent];
        assert_eq!(write_sgroup_hierarchy(&record, true, true), "SgH:2:0.1");
        assert_eq!(write_sgroup_hierarchy(&record, true, false), "");
        assert_eq!(write_sgroup_hierarchy(&record, false, true), "");
        assert_eq!(write_sgroup_hierarchy(&record, false, false), "");
    }

    #[test]
    fn review_hierarchy_many_typed_children_keep_row_order() {
        use cosmolkit_model::SubstanceGroupId;
        let mut record = parse("C");
        let count = 512;
        let parent_id = SubstanceGroupId::new(count - 1);
        record.topology.substance_groups = (0..count)
            .map(|row| {
                let mut group =
                    SubstanceGroup::new(SubstanceGroupId::new(row), SubstanceGroupKind::Data);
                group.set_prop("index", (1000 + row * 7).to_string());
                if row + 1 != count {
                    group.set_parent(parent_id);
                }
                group
            })
            .collect();
        let children = (0..count - 1)
            .map(|row| row.to_string())
            .collect::<Vec<_>>()
            .join(".");
        assert_eq!(
            write_sgroup_hierarchy(&record, true, false),
            format!("SgH:{}:{children}", count - 1)
        );
    }

    #[test]
    fn writes_large_ring_cis_trans_and_unknown_blocks() {
        for (input, expected) in [
            ("C1CCCC/C=C/CCC1 |t:5|", "C1=C/CCCCCCCC/1 |t:0|"),
            ("C1CCCCC=CCCC1 |c:5|", "C1=C\\CCCCCCCC/1 |c:0|"),
            ("C1=CCCCCCCCC1 |ctu:0|", "C1=CCCCCCCCC1 |ctu:0|"),
        ] {
            assert_eq!(write_cx_smiles(&parse(input)).unwrap(), expected, "{input}");
        }
    }

    #[test]
    fn writes_link_nodes_with_rdkit_canonical_index_behavior() {
        for (input, expected) in [
            ("OC1CCC(F)C1 |LN:1:1.3.2.6|", "OC1CCC(F)C1 |LN:1:1.3.2.6|"),
            (
                "FC1CCC(O)C1 |LN:1:1.3.2.6,4:1.4.3.6|",
                "OC1CCC(F)C1 |LN:4:1.3.3.6,1:1.4.2.6|",
            ),
            ("C1OCCC1C |LN:0:1.5,1:1.3|", "CC1CCOC1 |LN:5:1.5,4:1.3|"),
        ] {
            assert_eq!(write_cx_smiles(&parse(input)).unwrap(), expected, "{input}");
        }
    }

    #[test]
    fn field_selection_can_omit_coordinates_without_losing_other_records() {
        let record = parse("CC |(0,0,;1,0,),$_AV:left;right$|");
        let result = write_cx_smiles_with_params(
            &record,
            &CxSmilesWriteParams {
                fields: CxSmilesFields::ALL_BUT_COORDS,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(result, "CC |$_AV:left;right$|");
    }

    #[test]
    fn coordinate_numbers_use_c_percent_g_shape() {
        for (value, expected) in [
            (0.0, "0"),
            (0.0001, "0.0001"),
            (0.00123456789, "0.00123457"),
            (1.23456789, "1.23457"),
            (12_345.6789, "12345.7"),
            (123_456.789, "123457"),
            (1_234_567.89, "1.23457e+06"),
            (-1_234_567.89, "-1.23457e+06"),
            (999_999.9, "1e+06"),
        ] {
            assert_eq!(format_general(value), expected, "{value:?}");
        }
    }
}
