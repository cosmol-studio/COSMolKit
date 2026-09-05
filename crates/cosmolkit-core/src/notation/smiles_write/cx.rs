use super::*;

pub fn get_cx_extensions(molecule: &Molecule, fields: CxSmilesFields) -> Result<String, SmilesWriteError> {
    get_cx_extensions_scoped(molecule, fields, &CxWriteScope::full_molecule(molecule))
}

pub(super) fn get_cx_extensions_scoped(
    molecule: &Molecule,
    fields: CxSmilesFields,
    scope: &CxWriteScope,
) -> Result<String, SmilesWriteError> {
    write_cx_smiles_fields(molecule, fields, scope)
}

pub(super) fn write_cx_smiles_fields(
    molecule: &Molecule,
    fields: CxSmilesFields,
    scope: &CxWriteScope,
) -> Result<String, SmilesWriteError> {
    let mut res = String::from("|");
    let append_to_cx = |addition: &str, buf: &mut String| {
        if !addition.is_empty() {
            if buf.len() > 1 {
                buf.push(',');
            }
            buf.push_str(addition);
        }
    };

    if fields.contains(CxSmilesFields::COORDS) {
        let coords = write_cx_coordinates(molecule, &scope.atom_order);
        if !coords.is_empty() {
            res.push('(');
            res.push_str(&coords);
            res.push(')');
        }
    }

    let need_labels = scope.atom_order.iter().any(|atom_id| {
        let atom = &molecule.atoms()[atom_id.index()];
        atom.prop("atomLabel").is_some()
            || atom.prop("_QueryAtomGenericLabel").is_some()
            || atom.prop("dummyLabel").is_some()
            || atom.prop("_fromAttachPoint").is_some()
    });
    if fields.contains(CxSmilesFields::ATOM_LABELS) && need_labels {
        let labels = write_cx_atom_labels(molecule, &scope.atom_order);
        if !labels.is_empty() {
            append_to_cx(&format!("${}$", labels), &mut res);
        }
    }

    let need_values = scope
        .atom_order
        .iter()
        .any(|atom_id| molecule.atoms()[atom_id.index()].prop("molFileValue").is_some());
    if fields.contains(CxSmilesFields::MOLFILE_VALUES) && need_values {
        let values = write_cx_atom_values(molecule, &scope.atom_order);
        if !values.is_empty() {
            append_to_cx(&format!("$_AV:{}$", values), &mut res);
        }
    }

    if fields.contains(CxSmilesFields::RADICALS) {
        let radicals = write_cx_radicals(molecule, &scope.atom_order);
        if !radicals.is_empty() {
            if res.len() > 1 {
                res.push(',');
            }
            res.push_str(&radicals);
            if res.ends_with(',') {
                res.pop();
            }
        }
    }

    if fields.contains(CxSmilesFields::ATOM_PROPS) {
        let props = write_cx_atom_props(molecule, &scope.atom_order);
        append_to_cx(&props, &mut res);
    }

    if fields.contains(CxSmilesFields::BOND_CFG) {
        let include_coords = fields.contains(CxSmilesFields::COORDS) && molecule.coordinates_2d().is_some();
        let bond_cfg =
            write_cx_bond_config_block(molecule, &scope.atom_order, &scope.bond_order, include_coords, false);
        append_to_cx(&bond_cfg, &mut res);
        let ringbond_cistrans = write_cx_ringbond_cistrans_block(molecule, &scope.atom_order, &scope.bond_order);
        append_to_cx(&ringbond_cistrans, &mut res);
    } else if fields.contains(CxSmilesFields::BOND_ATROPISOMER) {
        let include_coords = fields.contains(CxSmilesFields::COORDS) && molecule.coordinates_2d().is_some();
        let bond_cfg = write_cx_bond_config_block(molecule, &scope.atom_order, &scope.bond_order, include_coords, true);
        append_to_cx(&bond_cfg, &mut res);
    }

    if fields.contains(CxSmilesFields::ENHANCED_STEREO) {
        let stereo = write_cx_enhanced_stereo(molecule, &scope.atom_order, &scope.bond_order);
        append_to_cx(&stereo, &mut res);
    }

    if fields.contains(CxSmilesFields::SGROUPS) {
        let sgroups = write_cx_sgroups(molecule, &scope.atom_order, &scope.bond_order);
        append_to_cx(&sgroups, &mut res);
    }

    if fields.contains(CxSmilesFields::POLYMER) {
        let polymer = write_cx_polymer_sgroups(molecule, &scope.atom_order, &scope.bond_order);
        append_to_cx(&polymer, &mut res);
    }

    if fields.contains(CxSmilesFields::SGROUPS) || fields.contains(CxSmilesFields::POLYMER) {
        let hierarchy = write_cx_sgroup_hierarchy_block(
            molecule,
            &scope.atom_order,
            &scope.bond_order,
            fields.contains(CxSmilesFields::SGROUPS),
            fields.contains(CxSmilesFields::POLYMER),
        );
        append_to_cx(&hierarchy, &mut res);
    }

    if fields.contains(CxSmilesFields::COORDINATE_BONDS) {
        let coord_bonds = write_cx_coordinate_bonds(molecule, &scope.atom_order, &scope.bond_order, "C");
        append_to_cx(&coord_bonds, &mut res);
    }

    if fields.contains(CxSmilesFields::HYDROGEN_BONDS) {
        let h_bonds = write_cx_coordinate_bonds(molecule, &scope.atom_order, &scope.bond_order, "H");
        append_to_cx(&h_bonds, &mut res);
    }

    if fields.contains(CxSmilesFields::ZERO_BONDS) {
        let zero_bonds = write_cx_zero_bonds(molecule, &scope.bond_order);
        append_to_cx(&zero_bonds, &mut res);
    }

    if fields.contains(CxSmilesFields::LINKNODES) {
        let linknodes = write_cx_linknodes_block(molecule, &scope.atom_order);
        append_to_cx(&linknodes, &mut res);
    }

    // RDKit: if (res.size() > 1) { res += "|"; } else { res = ""; }
    if res.len() > 1 {
        res.push('|');
    } else {
        res.clear();
    }
    Ok(res)
}

pub(super) fn cx_atom_output_positions(atom_order: &[AtomId], atom_count: usize) -> Vec<Option<usize>> {
    let mut positions = vec![None; atom_count];
    for (position, atom_id) in atom_order.iter().copied().enumerate() {
        if atom_id.index() < positions.len() {
            positions[atom_id.index()] = Some(position);
        }
    }
    positions
}

pub(super) fn cx_bond_output_positions(bond_order: &[BondId], bond_count: usize) -> Vec<Option<usize>> {
    let mut positions = vec![None; bond_count];
    for (position, bond_id) in bond_order.iter().copied().enumerate() {
        if bond_id.index() < positions.len() {
            positions[bond_id.index()] = Some(position);
        }
    }
    positions
}

pub(super) fn zero_small_writer_coord(value: f64) -> f64 {
    if value.abs() < 1e-4 { 0.0 } else { value }
}

pub(super) fn quote_atomprop_string(text: &str) -> String {
    text.chars()
        .map(|ch| if ch == '.' { "&#46;".to_string() } else { ch.to_string() })
        .collect()
}

pub(super) fn write_cx_coordinates(molecule: &Molecule, atom_order: &[AtomId]) -> String {
    let coords = match molecule.coordinates_2d() {
        Some(c) => c,
        None => return String::new(),
    };
    let mut parts = Vec::new();
    for atom in atom_order {
        let Some(coord) = coords.get(atom.index()) else {
            continue;
        };
        parts.push(format!(
            "{}, {},",
            zero_small_writer_coord(coord[0]),
            zero_small_writer_coord(coord[1])
        ));
    }
    parts
        .into_iter()
        .map(|part| part.replace(", ", ","))
        .collect::<Vec<_>>()
        .join(";")
}

pub(super) fn write_cx_coords(molecule: &Molecule, atom_order: &[AtomId]) -> String {
    write_cx_coordinates(molecule, atom_order)
}

pub(super) fn write_cx_atom_labels(molecule: &Molecule, atom_order: &[AtomId]) -> String {
    let pseudoatoms = ["Pol", "Mod", "Het", "Any", "A", "Q", "X", "*"];
    let mut parts = Vec::new();
    for atom_id in atom_order {
        let atom = &molecule.atoms()[atom_id.index()];
        let part = if let Some(label) = atom.prop("_QueryAtomGenericLabel") {
            Some(format!("{label}_p"))
        } else if atom.atomic_number() == 0
            && atom
                .prop("dummyLabel")
                .is_some_and(|label| pseudoatoms.contains(&label))
        {
            Some(format!("{}_p", atom.prop("dummyLabel").unwrap_or_default()))
        } else if atom.atomic_number() == 0
            && atom
                .prop("_fromAttachPoint")
                .and_then(|value| value.parse::<u32>().ok())
                .is_some_and(|value| value == 1 || value == 2)
        {
            Some(format!("_AP{}", atom.prop("_fromAttachPoint").unwrap_or_default()))
        } else {
            atom.prop("atomLabel").map(str::to_string)
        };
        if let Some(part) = part {
            parts.push(part);
        } else {
            parts.push(String::new());
        }
    }
    if parts.iter().all(|part| part.is_empty()) {
        String::new()
    } else {
        parts.join(";")
    }
}

pub(super) fn write_cx_atom_values(molecule: &Molecule, atom_order: &[AtomId]) -> String {
    let mut parts = Vec::new();
    for atom_id in atom_order {
        let atom = &molecule.atoms()[atom_id.index()];
        if let Some(value) = atom.prop("molFileValue") {
            parts.push(value.to_string());
        } else {
            parts.push(String::new());
        }
    }
    parts.join(";")
}

pub(super) fn write_cx_molfile_values(molecule: &Molecule, atom_order: &[AtomId]) -> String {
    write_cx_atom_values(molecule, atom_order)
}

pub(super) fn write_cx_radicals(molecule: &Molecule, atom_order: &[AtomId]) -> String {
    let mut by_count: BTreeMap<u8, Vec<usize>> = BTreeMap::new();
    for (output_idx, atom_id) in atom_order.iter().copied().enumerate() {
        let atom = &molecule.atoms()[atom_id.index()];
        let re = atom.radical_electrons();
        if re > 0 {
            by_count.entry(re).or_default().push(output_idx);
        }
    }
    if by_count.is_empty() {
        return String::new();
    }
    let mut result = String::new();
    for (count, atoms) in by_count {
        let code = match count {
            1 => "^1:",
            2 => "^2:",
            3 => "^5:",
            _ => continue,
        };
        result.push_str(code);
        for atom in atoms {
            result.push_str(&format!("{atom},"));
        }
    }
    result
}

pub(super) fn write_cx_atom_props(molecule: &Molecule, atom_order: &[AtomId]) -> String {
    let skip = [
        "atomLabel",
        "molFileValue",
        "molParity",
        "molStereoCare",
        "molRxnExactChange",
        "molInversionFlag",
        "dummyLabel",
    ];
    let mut result = String::new();
    for (which, atom_id) in atom_order.iter().copied().enumerate() {
        let atom = &molecule.atoms()[atom_id.index()];
        let is_attachment_point = atom.atomic_number() == 0 && atom.prop("_fromAttachPoint").is_some();
        for (prop_name, prop_value) in atom.props() {
            // RDKit getPropList(includePrivate=false, includeComputed=false)
            // excludes underscore-prefixed writer-internal/cache props from
            // CX atomProp emission.
            if prop_name.starts_with('_') {
                continue;
            }
            if skip.contains(&prop_name.as_str()) || prop_name == "molAtomMapNumber" {
                continue;
            }
            if prop_name == "dummyLabel"
                && (is_attachment_point
                    || prop_value == "*"
                    || ["Pol", "Mod", "Het", "Any", "A", "Q", "X", "*"].contains(&prop_value.as_str()))
            {
                continue;
            }
            if result.is_empty() {
                result.push_str("atomProp");
            }
            result.push_str(&format!(
                ":{which}.{}.{}",
                quote_atomprop_string(prop_name),
                quote_atomprop_string(prop_value)
            ));
        }
    }
    result
}

pub(super) fn write_cx_enhanced_stereo(molecule: &Molecule, atom_order: &[AtomId], _bond_order: &[BondId]) -> String {
    use crate::stereo::StereoGroupKind;
    let atom_positions = cx_atom_output_positions(atom_order, molecule.num_atoms());
    let write_ids = assigned_writer_stereo_group_ids(molecule.stereo_groups());
    let mut parts: Vec<String> = Vec::new();
    for (group, write_id) in molecule.stereo_groups().iter().zip(write_ids) {
        let mut atom_idxs: Vec<usize> = group
            .atoms()
            .iter()
            .filter_map(|atom| atom_positions.get(atom.index()).and_then(|position| *position))
            .collect();
        if atom_idxs.is_empty() {
            continue;
        }
        atom_idxs.sort_unstable();
        let prefix = match group.kind() {
            StereoGroupKind::Absolute => "a".to_string(),
            StereoGroupKind::Or => format!("o{}", write_id.unwrap_or(1)),
            StereoGroupKind::And => format!("&{}", write_id.unwrap_or(1)),
        };
        let members = atom_idxs.into_iter().map(|idx| idx.to_string()).collect::<Vec<_>>();
        parts.push(format!("{prefix}:{}", members.join(",")));
    }
    parts.join(",")
}

pub(super) fn assigned_writer_stereo_group_ids(groups: &[crate::StereoGroup]) -> Vec<Option<u32>> {
    use crate::stereo::StereoGroupKind;

    let mut or_ids = Vec::<u32>::new();
    let mut and_ids = Vec::<u32>::new();
    let mut assigned = groups.iter().map(crate::StereoGroup::id).collect::<Vec<_>>();
    for (idx, group) in groups.iter().enumerate() {
        let Some(id) = assigned[idx] else {
            continue;
        };
        let ids = match group.kind() {
            StereoGroupKind::Or => &mut or_ids,
            StereoGroupKind::And => &mut and_ids,
            StereoGroupKind::Absolute => continue,
        };
        if id != 0 && ids.contains(&id) {
            assigned[idx] = None;
        } else if id != 0 {
            ids.push(id);
        }
    }

    let mut next_or = 0_u32;
    let mut next_and = 0_u32;
    for (idx, group) in groups.iter().enumerate() {
        if group.kind() == StereoGroupKind::Absolute || assigned[idx].is_some() {
            continue;
        }
        let (next, ids) = match group.kind() {
            StereoGroupKind::Or => (&mut next_or, &mut or_ids),
            StereoGroupKind::And => (&mut next_and, &mut and_ids),
            StereoGroupKind::Absolute => unreachable!(),
        };
        *next += 1;
        while ids.contains(&*next) {
            *next += 1;
        }
        ids.push(*next);
        assigned[idx] = Some(*next);
    }
    assigned
}

pub(super) fn write_cx_sgroups(molecule: &Molecule, atom_order: &[AtomId], bond_order: &[BondId]) -> String {
    let data = write_cx_data_sgroups(molecule, atom_order);
    let other = write_cx_non_data_sgroups(molecule, atom_order, bond_order);
    match (data.is_empty(), other.is_empty()) {
        (true, true) => String::new(),
        (false, true) => data,
        (true, false) => other,
        (false, false) => format!("{data},{other}"),
    }
}

pub(super) fn write_cx_data_sgroups(molecule: &Molecule, atom_order: &[AtomId]) -> String {
    let atom_positions = cx_atom_output_positions(atom_order, molecule.num_atoms());
    let mut parts = Vec::new();
    for sgroup in molecule.substance_groups() {
        if !writer_is_data_sgroup(sgroup) {
            continue;
        }
        let atoms = sgroup
            .atoms()
            .iter()
            .filter_map(|atom| atom_positions.get(atom.index()).and_then(|value| *value))
            .map(|idx| idx.to_string())
            .collect::<Vec<_>>();
        if atoms.is_empty() {
            continue;
        }
        let field_name = writer_data_sgroup_field_name(sgroup);
        let data_fields = writer_data_sgroup_values(sgroup).join(",");
        let query_op = writer_data_sgroup_query_op(sgroup);
        let field_info = writer_data_sgroup_field_info(sgroup);
        let field_tag = writer_data_sgroup_field_tag(sgroup);
        parts.push(format!(
            "SgD:{}:{field_name}:{data_fields}:{query_op}:{field_info}:{field_tag}:",
            atoms.join(",")
        ));
    }
    parts.join(",")
}

pub(super) fn write_cx_non_data_sgroups(molecule: &Molecule, atom_order: &[AtomId], bond_order: &[BondId]) -> String {
    use crate::sgroup::SubstanceGroupKind;
    let atom_positions = cx_atom_output_positions(atom_order, molecule.num_atoms());
    let bond_positions = cx_bond_output_positions(bond_order, molecule.num_bonds());
    let mut parts = Vec::new();
    for sgroup in molecule.substance_groups() {
        if matches!(sgroup.kind(), SubstanceGroupKind::Data)
            || sgroup.props().get("TYPE").is_some_and(|value| value == "DAT")
            || writer_polymer_sgroup_type_code(sgroup).is_some()
        {
            continue;
        }
        let kind_str = match sgroup.kind() {
            SubstanceGroupKind::Data => "DAT",
            SubstanceGroupKind::Superatom => "SUP",
            SubstanceGroupKind::MultipleGroup => "MUL",
            SubstanceGroupKind::StructuralRepeatUnit => "SRU",
            SubstanceGroupKind::Monomer => "MON",
            SubstanceGroupKind::Copolymer => "COP",
            SubstanceGroupKind::Crosslink => "CRO",
            SubstanceGroupKind::Graft => "GRA",
            SubstanceGroupKind::Modification => "MOD",
            SubstanceGroupKind::Mer => "MER",
            SubstanceGroupKind::AnyPolymer => "ANY",
            SubstanceGroupKind::MixtureComponent => "MIX",
            SubstanceGroupKind::Mixture => "MIXTURE",
            SubstanceGroupKind::Formulation => "FOR",
            SubstanceGroupKind::Generic(s) => s.as_str(),
        };
        let atom_idxs: Vec<String> = sgroup
            .atoms()
            .iter()
            .filter_map(|atom| atom_positions.get(atom.index()).and_then(|value| *value))
            .map(|idx| idx.to_string())
            .collect();
        let bond_idxs: Vec<String> = sgroup
            .bonds()
            .iter()
            .filter_map(|bond| bond_positions.get(bond.index()).and_then(|value| *value))
            .map(|idx| idx.to_string())
            .collect();
        if atom_idxs.is_empty() && bond_idxs.is_empty() {
            continue;
        }
        let mut entry = format!("_S:{}:{}:{}", kind_str, atom_idxs.join(","), bond_idxs.join(","));
        if let Some(label) = sgroup.label() {
            entry.push(':');
            entry.push_str(&label.replace(',', "\\,").replace('|', "\\|"));
        }
        if let Some(conn) = sgroup.connection() {
            let conn_str = match conn {
                crate::sgroup::SGroupConnection::HeadToHead => "HH",
                crate::sgroup::SGroupConnection::HeadToTail => "HT",
                crate::sgroup::SGroupConnection::Either => "EU",
                crate::sgroup::SGroupConnection::Unknown(s) => s,
            };
            entry.push(':');
            entry.push_str(conn_str);
        }
        parts.push(entry);
    }
    parts.join(",")
}

pub(super) fn write_cx_coordinate_bonds(
    molecule: &Molecule,
    atom_order: &[AtomId],
    bond_order: &[BondId],
    symbol: &str,
) -> String {
    let atom_positions = cx_atom_output_positions(atom_order, molecule.num_atoms());
    let target_order = match symbol {
        "C" => BondOrder::Dative,
        "H" => BondOrder::Hydrogen,
        _ => return String::new(),
    };
    let mut parts = Vec::new();
    for (bond_output_idx, bond_id) in bond_order.iter().copied().enumerate() {
        let bond = &molecule.bonds()[bond_id.index()];
        let matches = if symbol == "C" {
            matches!(bond.order(), BondOrder::Dative | BondOrder::DativeOne)
        } else {
            bond.order() == target_order
        };
        if !matches {
            continue;
        }
        let Some(begin_output_idx) = atom_positions.get(bond.begin().index()).and_then(|value| *value) else {
            continue;
        };
        parts.push(format!("{begin_output_idx}.{bond_output_idx}"));
    }
    if parts.is_empty() {
        String::new()
    } else {
        format!("{symbol}:{}", parts.join(","))
    }
}

pub(super) fn write_cx_zero_bonds(molecule: &Molecule, bond_order: &[BondId]) -> String {
    let mut parts = Vec::new();
    for (bond_output_idx, bond_id) in bond_order.iter().copied().enumerate() {
        if molecule.bonds()[bond_id.index()].order() == BondOrder::Zero {
            parts.push(bond_output_idx.to_string());
        }
    }
    if parts.is_empty() {
        String::new()
    } else {
        format!("Z:{}", parts.join(","))
    }
}

// BEGIN RDKIT CPP FUNCTION get_bond_config_block
// RDKit✔️✔️: std::string get_bond_config_block(
// RDKit✔️✔️:     const ROMol &mol, const std::vector<unsigned int> &atomOrder,
// RDKit✔️✔️:     const std::vector<unsigned int> &bondOrder, bool coordsIncluded,
// RDKit✔️✔️:     std::map<int, std::unique_ptr<RDKit::Chirality::WedgeInfoBase>> &wedgeBonds,
// RDKit✔️✔️:     bool atropisomerOnly = false) {
// RDKit✔️✔️:   std::map<std::string, std ::vector<std::string>> wParts;
// RDKit✔️✔️:   for (unsigned int i = 0; i < bondOrder.size(); ++i) {
// RDKit✔️✔️:     auto idx = bondOrder[i];
// RDKit✔️✔️:     const auto bond = mol.getBondWithIdx(idx);
// RDKit✔️✔️:     unsigned int wedgeStartAtomIdx = bond->getBeginAtomIdx();
// RDKit✔️✔️:     if (!canHaveDirection(*bond)) { continue; }
// RDKit✔️✔️:     Bond::BondDir bd = bond->getBondDir();
// RDKit✔️✔️:     switch (bd) {
// RDKit✔️✔️:       case Bond::BondDir::BEGINDASH:
// RDKit✔️✔️:       case Bond::BondDir::BEGINWEDGE:
// RDKit✔️✔️:       case Bond::BondDir::UNKNOWN:
// RDKit✔️✔️:         break;
// RDKit✔️✔️:       default:
// RDKit✔️✔️:         bd = Bond::BondDir::NONE;
// RDKit✔️✔️:     }
// RDKit✔️✔️:     if (atropisomerOnly && bd == Bond::BondDir::NONE) { continue; }
// RDKit✔️✔️:     // ... atropisomer checks and optional wedge flipping ...
// RDKit✔️✔️:     if (!atropisomerOnly) {
// RDKit✔️✔️:       unsigned int cfg = 0;
// RDKit✔️✔️:       if (bd == Bond::BondDir::NONE &&
// RDKit✔️✔️:           bond->getPropIfPresent(common_properties::_MolFileBondCfg, cfg)) {
// RDKit✔️✔️:         switch (cfg) {
// RDKit✔️✔️:           case 1: bd = Bond::BondDir::BEGINWEDGE; break;
// RDKit✔️✔️:           case 2: bd = Bond::BondDir::UNKNOWN; break;
// RDKit✔️✔️:           case 3: bd = Bond::BondDir::BEGINDASH; break;
// RDKit✔️✔️:           default: bd = Bond::BondDir::NONE;
// RDKit✔️✔️:         }
// RDKit✔️✔️:       }
// RDKit✔️✔️:     }
// RDKit✔️✔️:     auto begAtomOrder =
// RDKit✔️✔️:         std::find(atomOrder.begin(), atomOrder.end(), wedgeStartAtomIdx) -
// RDKit✔️✔️:         atomOrder.begin();
// RDKit✔️✔️:     std::string wType = "";
// RDKit✔️✔️:     if (bd == Bond::BondDir::UNKNOWN) { wType = "w"; }
// RDKit✔️✔️:     else if (coordsIncluded || isAnAtropisomer) {
// RDKit✔️✔️:       if (bd == Bond::BondDir::BEGINWEDGE) { wType = "wU"; }
// RDKit✔️✔️:       else if (bd == Bond::BondDir::BEGINDASH) { wType = "wD"; }
// RDKit✔️✔️:     }
// RDKit✔️✔️:     if (wType != "") { wParts[wType].push_back(format("%d.%d", begAtomOrder, i)); }
// RDKit✔️✔️:   }
// RDKit✔️✔️:   // join as "w:...,wD:...,wU:..."
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION get_bond_config_block
pub(super) fn write_cx_bond_config_block(
    molecule: &Molecule,
    atom_order: &[AtomId],
    bond_order: &[BondId],
    coords_included: bool,
    atropisomer_only: bool,
) -> String {
    let mut atom_order_positions: Vec<Option<usize>> = vec![None; molecule.atoms().len()];
    for (position, atom_id) in atom_order.iter().copied().enumerate() {
        if atom_id.index() < atom_order_positions.len() {
            atom_order_positions[atom_id.index()] = Some(position);
        }
    }

    let mut w_parts: BTreeMap<&'static str, Vec<String>> = BTreeMap::new();
    for (bond_output_idx, bond_id) in bond_order.iter().copied().enumerate() {
        let bond = &molecule.bonds()[bond_id.index()];
        let wedge_start_atom = bond.begin();
        if !can_have_direction_for_writer(bond.order()) {
            continue;
        }

        let mut direction = normalize_writer_cx_bond_direction(bond.direction());
        let mut is_an_atropisomer = false;

        if atropisomer_only && direction == BondDirection::None {
            continue;
        }

        if matches!(direction, BondDirection::BeginDash | BondDirection::BeginWedge) {
            for neighbor_bond_id in incident_bonds(molecule, wedge_start_atom) {
                if neighbor_bond_id == bond_id {
                    continue;
                }
                let neighbor = &molecule.bonds()[neighbor_bond_id.index()];
                if matches!(neighbor.stereo(), BondStereo::AtropCw | BondStereo::AtropCcw) {
                    is_an_atropisomer = true;
                    break;
                }
            }
        }

        if atropisomer_only {
            if !is_an_atropisomer {
                continue;
            }
        } else if matches!(direction, BondDirection::None)
            && let Some(cfg) = writer_parse_molfile_bond_cfg(bond)
        {
            direction = match cfg {
                1 => BondDirection::BeginWedge,
                2 => BondDirection::Unknown,
                3 => BondDirection::BeginDash,
                _ => BondDirection::None,
            };
        }

        // RDKit also derives direction from coordinates via
        // Chirality::GetMolFileBondStereoInfo when available; that path is
        // deferred until the full wedgeBonds/conformer parity port.
        let w_type = if direction == BondDirection::Unknown {
            Some("w")
        } else if coords_included || is_an_atropisomer {
            match direction {
                BondDirection::BeginWedge => Some("wU"),
                BondDirection::BeginDash => Some("wD"),
                _ => None,
            }
        } else {
            None
        };

        let Some(w_type) = w_type else {
            continue;
        };

        let Some(Some(begin_atom_order_idx)) = atom_order_positions.get(wedge_start_atom.index()) else {
            continue;
        };
        w_parts
            .entry(w_type)
            .or_default()
            .push(format!("{begin_atom_order_idx}.{bond_output_idx}"));
    }

    let mut parts: Vec<String> = Vec::new();
    for (w_type, entries) in w_parts {
        if entries.is_empty() {
            continue;
        }
        parts.push(format!("{w_type}:{}", entries.join(",")));
    }
    parts.join(",")
}

pub(super) fn normalize_writer_cx_bond_direction(direction: BondDirection) -> BondDirection {
    match direction {
        BondDirection::BeginDash | BondDirection::BeginWedge | BondDirection::Unknown => direction,
        _ => BondDirection::None,
    }
}

pub(super) fn writer_parse_molfile_bond_cfg(bond: &Bond) -> Option<u32> {
    bond.prop("_MolFileBondCfg").and_then(|value| value.parse::<u32>().ok())
}

// BEGIN RDKIT CPP FUNCTION get_ringbond_cistrans_block
// RDKit✔️✔️: std::string get_ringbond_cistrans_block(
// RDKit✔️✔️:     const ROMol &mol, const std::vector<unsigned int> &atomOrder,
// RDKit✔️✔️:     const std::vector<unsigned int> &bondOrder) {
// RDKit✔️✔️:   if (!mol.getRingInfo()->isInitialized()) { return ""; }
// RDKit✔️✔️:   std::string c = "", t = "", ctu = "";
// RDKit✔️✔️:   for (unsigned int i = 0; i < bondOrder.size(); ++i) {
// RDKit✔️✔️:     auto idx = bondOrder[i];
// RDKit✔️✔️:     if (!rinfo->numBondRings(idx) ||
// RDKit✔️✔️:         rinfo->minBondRingSize(idx) < Chirality::minRingSizeForDoubleBondStereo) { continue; }
// RDKit✔️✔️:     if (bond->getBondType() != Bond::DOUBLE && bond->getBondType() != Bond::AROMATIC) { continue; }
// RDKit✔️✔️:     if (bstereo != Bond::STEREOANY && bstereo != Bond::STEREOCIS && bstereo != Bond::STEREOTRANS) { continue; }
// RDKit✔️✔️:     if (bstereo == Bond::STEREOANY) { ctu += label; } else {
// RDKit✔️✔️:       bool needSwap = false;
// RDKit✔️✔️:       // parity flips from atom-order reordering around stereo refs
// RDKit✔️✔️:       if (bstereo == Bond::STEREOCIS || needSwap) { c += label; } else { t += label; }
// RDKit✔️✔️:     }
// RDKit✔️✔️:   }
// RDKit✔️✔️:   return c + t + ctu;
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION get_ringbond_cistrans_block
pub(super) fn write_cx_ringbond_cistrans_block(
    molecule: &Molecule,
    atom_order: &[AtomId],
    bond_order: &[BondId],
) -> String {
    const MIN_RING_SIZE_FOR_DOUBLE_BOND_STEREO: usize = 8;
    let Some(ring_info) = molecule.derived_cache().rings.as_ref() else {
        return String::new();
    };
    if !ring_info.is_initialized() {
        return String::new();
    }

    let mut atom_order_positions: Vec<Option<usize>> = vec![None; molecule.atoms().len()];
    for (position, atom_id) in atom_order.iter().copied().enumerate() {
        if atom_id.index() < atom_order_positions.len() {
            atom_order_positions[atom_id.index()] = Some(position);
        }
    }

    let mut c_labels: Vec<String> = Vec::new();
    let mut t_labels: Vec<String> = Vec::new();
    let mut ctu_labels: Vec<String> = Vec::new();

    for (bond_output_idx, bond_id) in bond_order.iter().copied().enumerate() {
        if ring_info.num_bond_rings(bond_id) == 0
            || ring_info.min_bond_ring_size(bond_id) < MIN_RING_SIZE_FOR_DOUBLE_BOND_STEREO
        {
            continue;
        }
        let bond = &molecule.bonds()[bond_id.index()];
        if !matches!(bond.order(), BondOrder::Double | BondOrder::Aromatic) {
            continue;
        }
        if !matches!(bond.stereo(), BondStereo::Any | BondStereo::Cis | BondStereo::Trans) {
            continue;
        }

        let label = bond_output_idx.to_string();
        if bond.stereo() == BondStereo::Any {
            ctu_labels.push(label);
            continue;
        }

        let Some([stereo_begin, stereo_end]) = bond.stereo_atoms() else {
            continue;
        };

        let begin_atom = bond.begin();
        let end_atom = bond.end();
        let mut need_swap = false;

        if incident_bonds(molecule, begin_atom).len() > 2 {
            let Some(stereo_begin_order) = atom_order_positions
                .get(stereo_begin.index())
                .and_then(|position| *position)
            else {
                continue;
            };
            for neighbor_bond in incident_bonds(molecule, begin_atom) {
                if neighbor_bond == bond_id {
                    continue;
                }
                let Some(neighbor_atom) = bond_other_atom(&molecule.bonds()[neighbor_bond.index()], begin_atom) else {
                    continue;
                };
                if neighbor_atom == end_atom || neighbor_atom == stereo_begin {
                    continue;
                }
                if atom_order_positions
                    .get(neighbor_atom.index())
                    .and_then(|position| *position)
                    .is_some_and(|neighbor_order| neighbor_order < stereo_begin_order)
                {
                    need_swap = !need_swap;
                }
            }
        }

        if incident_bonds(molecule, end_atom).len() > 2 {
            let Some(stereo_end_order) = atom_order_positions
                .get(stereo_end.index())
                .and_then(|position| *position)
            else {
                continue;
            };
            for neighbor_bond in incident_bonds(molecule, end_atom) {
                if neighbor_bond == bond_id {
                    continue;
                }
                let Some(neighbor_atom) = bond_other_atom(&molecule.bonds()[neighbor_bond.index()], end_atom) else {
                    continue;
                };
                if neighbor_atom == begin_atom || neighbor_atom == stereo_end {
                    continue;
                }
                if atom_order_positions
                    .get(neighbor_atom.index())
                    .and_then(|position| *position)
                    .is_some_and(|neighbor_order| neighbor_order < stereo_end_order)
                {
                    need_swap = !need_swap;
                }
            }
        }

        if bond.stereo() == BondStereo::Cis || need_swap {
            c_labels.push(label);
        } else {
            t_labels.push(label);
        }
    }

    let c = if c_labels.is_empty() {
        String::new()
    } else {
        format!("c:{}", c_labels.join(","))
    };
    let t = if t_labels.is_empty() {
        String::new()
    } else {
        format!("t:{}", t_labels.join(","))
    };
    let ctu = if ctu_labels.is_empty() {
        String::new()
    } else {
        format!("ctu:{}", ctu_labels.join(","))
    };
    format!("{c}{t}{ctu}")
}

// BEGIN RDKIT CPP FUNCTION get_linknodes_block
// RDKit✔️✔️: std::string get_linknodes_block(const ROMol &mol,
// RDKit✔️✔️:                                 const std::vector<unsigned int> &atomOrder) {
// RDKit✔️✔️:   bool strict = false;
// RDKit✔️✔️:   auto linkNodes = MolEnumerator::utils::getMolLinkNodes(mol, strict);
// RDKit✔️✔️:   if (linkNodes.empty()) { return ""; }
// RDKit✔️✔️:   std::stringstream res;
// RDKit✔️✔️:   res << "LN:";
// RDKit✔️✔️:   for (const auto &ln : linkNodes) {
// RDKit✔️✔️:     unsigned int atomIdx = atomOrder[ln.bondAtoms[0].first];
// RDKit✔️✔️:     res << atomIdx << ":" << ln.minRep << "." << ln.maxRep;
// RDKit✔️✔️:     if (mol.getAtomWithIdx(ln.bondAtoms[0].first)->getDegree() > 2) {
// RDKit✔️✔️:       res << "." << atomOrder[ln.bondAtoms[0].second] << "."
// RDKit✔️✔️:           << atomOrder[ln.bondAtoms[1].second];
// RDKit✔️✔️:     }
// RDKit✔️✔️:     res << ",";
// RDKit✔️✔️:   }
// RDKit✔️✔️:   if (!resStr.empty() && resStr.back() == ',') { resStr.pop_back(); }
// RDKit✔️✔️:   return resStr;
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION get_linknodes_block
pub(super) fn write_cx_linknodes_block(molecule: &Molecule, atom_order: &[AtomId]) -> String {
    let Some(raw_link_nodes) = molecule.prop("_MolFileLinkNodes") else {
        return String::new();
    };

    let mut atom_order_positions: Vec<Option<usize>> = vec![None; molecule.atoms().len()];
    for (position, atom_id) in atom_order.iter().copied().enumerate() {
        if atom_id.index() < atom_order_positions.len() {
            atom_order_positions[atom_id.index()] = Some(position);
        }
    }

    let mut entries: Vec<String> = Vec::new();
    for record in raw_link_nodes.split('|').filter(|part| !part.trim().is_empty()) {
        let values: Vec<usize> = record
            .split_whitespace()
            .filter_map(|part| part.parse::<usize>().ok())
            .collect();
        if values.len() < 5 {
            continue;
        }
        let min_rep = values[0];
        let max_rep = values[1];
        let pair_count = values[2];
        let required = 3usize.saturating_add(pair_count.saturating_mul(2));
        if values.len() < required || pair_count == 0 {
            continue;
        }

        let Some(center_atom_one_based) = values.get(3).copied() else {
            continue;
        };
        let Some(center_atom_idx) = center_atom_one_based.checked_sub(1) else {
            continue;
        };
        let Some(center_output_idx) = atom_order_positions.get(center_atom_idx).and_then(|position| *position) else {
            continue;
        };

        let mut entry = format!("{center_output_idx}:{min_rep}.{max_rep}");
        let center_atom = AtomId::new(center_atom_idx);
        if incident_bonds(molecule, center_atom).len() > 2 && pair_count >= 2 {
            let Some(first_neighbor_one_based) = values.get(4).copied() else {
                continue;
            };
            let Some(second_neighbor_one_based) = values.get(6).copied() else {
                continue;
            };
            let Some(first_neighbor_idx) = first_neighbor_one_based.checked_sub(1) else {
                continue;
            };
            let Some(second_neighbor_idx) = second_neighbor_one_based.checked_sub(1) else {
                continue;
            };
            let Some(first_neighbor_out) = atom_order_positions
                .get(first_neighbor_idx)
                .and_then(|position| *position)
            else {
                continue;
            };
            let Some(second_neighbor_out) = atom_order_positions
                .get(second_neighbor_idx)
                .and_then(|position| *position)
            else {
                continue;
            };
            entry.push_str(&format!(".{first_neighbor_out}.{second_neighbor_out}"));
        }
        entries.push(entry);
    }

    if entries.is_empty() {
        String::new()
    } else {
        format!("LN:{}", entries.join(","))
    }
}

// BEGIN RDKIT CPP FUNCTION get_sgroup_polymer_block
// RDKit✔️✔️: std::string get_sgroup_polymer_block(
// RDKit✔️✔️:     const ROMol &mol, const std::vector<unsigned int> &atomOrder,
// RDKit✔️✔️:     const std::vector<unsigned int> &bondOrder) {
// RDKit✔️✔️:   // builds Sg:type:atoms:label:connect:headCrossings:tailCrossings:
// RDKit✔️✔️:   // using reverse typemap and output-order indices
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION get_sgroup_polymer_block
pub(super) fn write_cx_sgroup_polymer_block(
    molecule: &Molecule,
    atom_order: &[AtomId],
    bond_order: &[BondId],
) -> String {
    if molecule.substance_groups().is_empty() {
        return String::new();
    }

    let rev_atom_order = cx_atom_output_positions(atom_order, molecule.num_atoms());
    let rev_bond_order = cx_bond_output_positions(bond_order, molecule.num_bonds());

    let mut entries: Vec<String> = Vec::new();
    for sgroup in molecule.substance_groups() {
        let Some(type_code) = writer_polymer_sgroup_type_code(sgroup) else {
            continue;
        };
        if sgroup.atoms().is_empty() {
            continue;
        }

        let mut atom_parts: Vec<String> = Vec::new();
        for atom in sgroup.atoms() {
            let Some(out_idx) = rev_atom_order.get(atom.index()).and_then(|v| *v) else {
                continue;
            };
            atom_parts.push(out_idx.to_string());
        }
        if atom_parts.is_empty() {
            continue;
        }

        let label = sgroup.label().unwrap_or_default();
        let connect = sgroup
            .connection()
            .map(|value| match value {
                crate::sgroup::SGroupConnection::HeadToHead => "hh".to_string(),
                crate::sgroup::SGroupConnection::HeadToTail => "ht".to_string(),
                crate::sgroup::SGroupConnection::Either => "eu".to_string(),
                crate::sgroup::SGroupConnection::Unknown(text) => text.to_ascii_lowercase(),
            })
            .or_else(|| sgroup.props().get("CONNECT").map(|value| value.to_ascii_lowercase()))
            .unwrap_or_default();

        let head_crossings = writer_parse_sgroup_crossings(sgroup, "_headCrossings")
            .or_else(|| writer_parse_sgroup_crossings(sgroup, "XBHEAD"))
            .unwrap_or_default();
        let tail_crossings = writer_parse_sgroup_crossings(sgroup, "_tailCrossings")
            .or_else(|| writer_parse_sgroup_crossings(sgroup, "XBCORR"))
            .unwrap_or_default();

        let head_field = if head_crossings.len() > 1 {
            head_crossings
                .iter()
                .filter_map(|idx| rev_bond_order.get(*idx).and_then(|value| *value))
                .map(|idx| idx.to_string())
                .collect::<Vec<_>>()
                .join(",")
        } else {
            String::new()
        };
        let tail_field = if tail_crossings.len() > 2 {
            tail_crossings
                .iter()
                .enumerate()
                .filter_map(|(i, idx)| {
                    if i % 2 == 1 {
                        rev_bond_order
                            .get(*idx)
                            .and_then(|value| *value)
                            .map(|idx| idx.to_string())
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>()
                .join(",")
        } else {
            String::new()
        };

        entries.push(format!(
            "Sg:{type_code}:{}:{label}:{connect}:{head_field}:{tail_field}:",
            atom_parts.join(",")
        ));
    }

    entries.join(",")
}

pub(super) fn write_cx_polymer_sgroups(molecule: &Molecule, atom_order: &[AtomId], bond_order: &[BondId]) -> String {
    write_cx_sgroup_polymer_block(molecule, atom_order, bond_order)
}

pub(super) fn writer_is_data_sgroup(sgroup: &crate::SubstanceGroup) -> bool {
    matches!(sgroup.kind(), crate::sgroup::SubstanceGroupKind::Data)
        || sgroup.props().get("TYPE").is_some_and(|value| value == "DAT")
}

pub(super) fn writer_data_sgroup_field_name(sgroup: &crate::SubstanceGroup) -> String {
    sgroup
        .data()
        .and_then(|data| data.field_name.clone())
        .or_else(|| sgroup.props().get("FIELDNAME").cloned())
        .unwrap_or_default()
}

pub(super) fn writer_data_sgroup_values(sgroup: &crate::SubstanceGroup) -> Vec<String> {
    if !sgroup.data_fields().is_empty() {
        return sgroup.data_fields().to_vec();
    }
    if let Some(values) = sgroup
        .data()
        .map(|data| {
            data.values
                .iter()
                .filter(|value| !value.is_empty())
                .cloned()
                .collect::<Vec<_>>()
        })
        .filter(|values| !values.is_empty())
    {
        return values;
    }
    sgroup
        .props()
        .get("DATAFIELDS")
        .map(|value| vec![value.clone()])
        .unwrap_or_default()
}

pub(super) fn writer_data_sgroup_query_op(sgroup: &crate::SubstanceGroup) -> String {
    sgroup
        .data()
        .and_then(|data| data.query_op.clone())
        .or_else(|| sgroup.props().get("QUERYOP").cloned())
        .unwrap_or_default()
}

pub(super) fn writer_data_sgroup_field_info(sgroup: &crate::SubstanceGroup) -> String {
    sgroup
        .data()
        .and_then(|data| data.field_info.clone())
        .or_else(|| sgroup.props().get("FIELDINFO").cloned())
        .unwrap_or_default()
}

pub(super) fn writer_data_sgroup_field_tag(sgroup: &crate::SubstanceGroup) -> String {
    sgroup
        .data()
        .and_then(|data| data.units.clone())
        .or_else(|| sgroup.props().get("FIELDTAG").cloned())
        .unwrap_or_default()
}

pub(super) fn writer_polymer_sgroup_type_code(sgroup: &crate::SubstanceGroup) -> Option<String> {
    let code = match sgroup.kind() {
        crate::sgroup::SubstanceGroupKind::StructuralRepeatUnit => "n",
        crate::sgroup::SubstanceGroupKind::Monomer => "mon",
        crate::sgroup::SubstanceGroupKind::Mer => "mer",
        crate::sgroup::SubstanceGroupKind::Copolymer => {
            let subtype = sgroup
                .subtype()
                .map(|value| value.to_ascii_uppercase())
                .or_else(|| sgroup.props().get("SUBTYPE").map(|value| value.to_ascii_uppercase()));
            match subtype.as_deref() {
                Some("ALT") => "alt",
                Some("RAN") => "ran",
                Some("BLO") | Some("BLK") => "blk",
                _ => "co",
            }
        }
        crate::sgroup::SubstanceGroupKind::Crosslink => "xl",
        crate::sgroup::SubstanceGroupKind::Modification => "mod",
        crate::sgroup::SubstanceGroupKind::MixtureComponent => "mix",
        crate::sgroup::SubstanceGroupKind::Formulation => "f",
        crate::sgroup::SubstanceGroupKind::AnyPolymer => "any",
        crate::sgroup::SubstanceGroupKind::Graft => "grf",
        crate::sgroup::SubstanceGroupKind::Generic(value) if value.eq_ignore_ascii_case("GEN") => "gen",
        crate::sgroup::SubstanceGroupKind::Generic(value) if value.eq_ignore_ascii_case("COM") => "c",
        _ => return None,
    };
    Some(code.to_string())
}

pub(super) fn writer_parse_sgroup_crossings(sgroup: &crate::SubstanceGroup, key: &str) -> Option<Vec<usize>> {
    let raw = sgroup.props().get(key)?;
    let parsed: Vec<usize> = raw
        .split(',')
        .filter_map(|part| part.trim().parse::<usize>().ok())
        .collect();
    Some(parsed)
}

// BEGIN RDKIT CPP FUNCTION get_sgroup_hierarchy_block
// RDKit✔️✔️: std::string get_sgroup_hierarchy_block(const ROMol &mol) {
// RDKit✔️✔️:   // builds SgH:parentOutputIdx:childOutputIdx.childOutputIdx,...
// RDKit✔️✔️:   // using temporary _cxsmilesOutputIndex assigned by prior SGroup writers
// RDKit✔️✔️: }
// END RDKIT CPP FUNCTION get_sgroup_hierarchy_block
pub(super) fn write_cx_sgroup_hierarchy_block(
    molecule: &Molecule,
    atom_order: &[AtomId],
    bond_order: &[BondId],
    include_sgroups: bool,
    include_polymer: bool,
) -> String {
    let sgroups = molecule.substance_groups();
    if sgroups.is_empty() {
        return String::new();
    }
    let output_index_by_sgroup_id =
        writer_cx_hierarchy_output_indices(molecule, atom_order, bond_order, include_sgroups, include_polymer);

    if output_index_by_sgroup_id.is_empty() {
        return String::new();
    }

    let mut accum: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (fallback_idx, sgroup) in sgroups.iter().enumerate() {
        let child_id = writer_sgroup_index_value(sgroup, fallback_idx);
        let Some(child_output_idx) = output_index_by_sgroup_id.get(&child_id).copied() else {
            continue;
        };
        let Some(parent_id) = writer_sgroup_parent_value(sgroup) else {
            continue;
        };
        let Some(parent_output_idx) = output_index_by_sgroup_id.get(&parent_id).copied() else {
            continue;
        };
        accum.entry(parent_output_idx).or_default().push(child_output_idx);
    }

    if accum.is_empty() {
        return String::new();
    }

    let mut entries: Vec<String> = Vec::new();
    for (parent, children) in accum {
        if children.is_empty() {
            continue;
        }
        let child_text = children
            .iter()
            .map(|child| child.to_string())
            .collect::<Vec<_>>()
            .join(".");
        entries.push(format!("{parent}:{child_text}"));
    }

    if entries.is_empty() {
        String::new()
    } else {
        format!("SgH:{}", entries.join(","))
    }
}

pub(super) fn writer_cx_hierarchy_output_indices(
    molecule: &Molecule,
    atom_order: &[AtomId],
    bond_order: &[BondId],
    include_sgroups: bool,
    include_polymer: bool,
) -> BTreeMap<usize, usize> {
    let atom_set = atom_order.iter().copied().collect::<BTreeSet<_>>();
    let bond_set = bond_order.iter().copied().collect::<BTreeSet<_>>();
    let mut output_index_by_sgroup_id = BTreeMap::new();
    let mut next_output_index = 0usize;

    if include_sgroups {
        for (fallback_idx, sgroup) in molecule.substance_groups().iter().enumerate() {
            if !writer_is_data_sgroup(sgroup) || !sgroup.atoms().iter().any(|atom| atom_set.contains(atom)) {
                continue;
            }
            let sgroup_id = writer_sgroup_index_value(sgroup, fallback_idx);
            if let std::collections::btree_map::Entry::Vacant(entry) = output_index_by_sgroup_id.entry(sgroup_id) {
                entry.insert(next_output_index);
                next_output_index += 1;
            }
        }
        for (fallback_idx, sgroup) in molecule.substance_groups().iter().enumerate() {
            if writer_is_data_sgroup(sgroup)
                || writer_polymer_sgroup_type_code(sgroup).is_some()
                || (!sgroup.atoms().iter().any(|atom| atom_set.contains(atom))
                    && !sgroup.bonds().iter().any(|bond| bond_set.contains(bond)))
            {
                continue;
            }
            let sgroup_id = writer_sgroup_index_value(sgroup, fallback_idx);
            if let std::collections::btree_map::Entry::Vacant(entry) = output_index_by_sgroup_id.entry(sgroup_id) {
                entry.insert(next_output_index);
                next_output_index += 1;
            }
        }
    }

    if include_polymer {
        for (fallback_idx, sgroup) in molecule.substance_groups().iter().enumerate() {
            if writer_polymer_sgroup_type_code(sgroup).is_none()
                || !sgroup.atoms().iter().any(|atom| atom_set.contains(atom))
            {
                continue;
            }
            let sgroup_id = writer_sgroup_index_value(sgroup, fallback_idx);
            if let std::collections::btree_map::Entry::Vacant(entry) = output_index_by_sgroup_id.entry(sgroup_id) {
                entry.insert(next_output_index);
                next_output_index += 1;
            }
        }
    }

    output_index_by_sgroup_id
}

pub(super) fn writer_sgroup_index_value(sgroup: &crate::SubstanceGroup, _fallback_idx: usize) -> usize {
    sgroup
        .props()
        .get("index")
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or_else(|| sgroup.id().index())
}

pub(super) fn writer_sgroup_parent_value(sgroup: &crate::SubstanceGroup) -> Option<usize> {
    if let Some(parent) = sgroup.parent() {
        return Some(parent.index());
    }
    sgroup
        .props()
        .get("PARENT")
        .and_then(|value| value.parse::<usize>().ok())
}
