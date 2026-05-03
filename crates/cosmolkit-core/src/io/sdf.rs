use std::io::BufRead;

use crate::{
    AdjacencyList, Atom, Bond, BondDirection, BondOrder, BondStereo, ChiralTag,
    CoordinateDimension, Molecule,
};
use glam::{DVec2, DVec3};

/// One record extracted from an SDF stream.
#[derive(Debug, Clone, PartialEq)]
pub struct SdfRecord {
    pub title: String,
    pub program_line: Option<String>,
    pub comment_line: Option<String>,
    pub molecule: Molecule,
    pub data_fields: Vec<(String, String)>,
    pub raw_molblock: String,
}

/// Errors returned by SDF reading APIs.
#[derive(Debug, thiserror::Error)]
pub enum SdfReadError {
    #[error(transparent)]
    Io(#[from] std::io::Error),
    #[error("{0}")]
    Parse(String),
    #[error("SDF reader path is not implemented")]
    NotImplemented,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, Default)]
pub enum SdfCoordinateMode {
    #[default]
    Auto,
    Force2D,
    Force3D,
}

/// Streaming-capable SDF reader.
pub struct SdfReader<R> {
    reader: R,
    line_number: usize,
    at_end: bool,
    coordinate_mode: SdfCoordinateMode,
}

impl<R: BufRead> SdfReader<R> {
    /// Create a reader from any buffered input stream.
    #[must_use]
    pub fn new(reader: R) -> Self {
        Self::with_coordinate_mode(reader, SdfCoordinateMode::Auto)
    }

    #[must_use]
    pub fn with_coordinate_mode(reader: R, coordinate_mode: SdfCoordinateMode) -> Self {
        Self {
            reader,
            line_number: 0,
            at_end: false,
            coordinate_mode,
        }
    }

    /// Read next SDF record from stream.
    pub fn next_record(&mut self) -> Result<Option<SdfRecord>, SdfReadError> {
        let Some(lines) = self.read_next_raw_record()? else {
            return Ok(None);
        };

        // RDKit's forward supplier yields a null molecule for blank records
        // instead of manufacturing an empty molecule.
        if lines.iter().all(|line| strip_rdkit(line).is_empty()) {
            return Ok(None);
        }

        let title = lines.first().cloned().unwrap_or_default();
        let program_line = lines.get(1).cloned();
        let comment_line = lines.get(2).cloned();
        let source_dim = resolve_coordinate_dimension(
            self.coordinate_mode,
            infer_header_coordinate_dimension(lines.get(1).map(String::as_str)),
        );
        let (mut molecule, data_start) = parse_mol_data_stream(&lines, source_dim)?;
        let data_fields = parse_sdf_data_fields(&lines[data_start..])?;
        molecule.rebuild_adjacency();
        let raw_molblock = lines[..data_start].join("\n");

        Ok(Some(SdfRecord {
            title,
            program_line,
            comment_line,
            molecule,
            data_fields,
            raw_molblock,
        }))
    }

    fn read_next_raw_record(&mut self) -> Result<Option<Vec<String>>, SdfReadError> {
        if self.at_end {
            return Ok(None);
        }

        let mut lines = Vec::new();
        let mut line = String::new();
        loop {
            line.clear();
            let bytes_read = self.reader.read_line(&mut line)?;
            if bytes_read == 0 {
                self.at_end = true;
                return if lines.is_empty() {
                    Ok(None)
                } else {
                    Ok(Some(lines))
                };
            }
            self.line_number += 1;

            let record_line = trim_line_ending(&line).to_owned();
            if record_line.starts_with("$$$$") {
                return Ok(Some(lines));
            }
            lines.push(record_line);
        }
    }
}

pub fn read_sdf_record_from_str(s: &str) -> Result<SdfRecord, SdfReadError> {
    read_sdf_record_from_str_with_coordinate_mode(s, SdfCoordinateMode::Auto)
}

pub fn read_sdf_record_from_str_with_coordinate_mode(
    s: &str,
    coordinate_mode: SdfCoordinateMode,
) -> Result<SdfRecord, SdfReadError> {
    let mut reader =
        SdfReader::with_coordinate_mode(std::io::Cursor::new(s.as_bytes()), coordinate_mode);
    reader
        .next_record()?
        .ok_or_else(|| SdfReadError::Parse("SDF text did not contain a complete record".to_owned()))
}

pub fn read_sdf_records_from_str(s: &str) -> Result<Vec<SdfRecord>, SdfReadError> {
    read_sdf_records_from_str_with_coordinate_mode(s, SdfCoordinateMode::Auto)
}

pub fn read_sdf_records_from_str_with_coordinate_mode(
    s: &str,
    coordinate_mode: SdfCoordinateMode,
) -> Result<Vec<SdfRecord>, SdfReadError> {
    let mut reader =
        SdfReader::with_coordinate_mode(std::io::Cursor::new(s.as_bytes()), coordinate_mode);
    let mut records = Vec::new();
    while let Some(record) = reader.next_record()? {
        records.push(record);
    }
    Ok(records)
}

fn parse_mol_data_stream(
    lines: &[String],
    source_dim: Option<CoordinateDimension>,
) -> Result<(Molecule, usize), SdfReadError> {
    let counts = lines.get(3).ok_or_else(|| {
        SdfReadError::Parse(format!(
            "Counts line too short: '{}' on line4",
            lines.get(3).map_or("", String::as_str)
        ))
    })?;
    if counts.len() > 35 && counts.len() >= 39 && counts.as_bytes().get(34) == Some(&b'V') {
        match &counts[34..39] {
            "V2000" => parse_v2000_mol_data_stream(lines, source_dim),
            "V3000" => parse_v3000_mol_data_stream(lines, source_dim),
            version => Err(SdfReadError::Parse(format!(
                "Unsupported CTAB version: '{version}' at line 4"
            ))),
        }
    } else {
        parse_v2000_mol_data_stream(lines, source_dim)
    }
}

fn parse_v2000_mol_data_stream(
    lines: &[String],
    source_dim: Option<CoordinateDimension>,
) -> Result<(Molecule, usize), SdfReadError> {
    if lines.len() < 4 {
        return Err(SdfReadError::Parse(format!(
            "Counts line too short: '{}' on line4",
            lines.get(3).map_or("", String::as_str)
        )));
    }

    let counts = &lines[3];
    if counts.len() < 6 {
        return Err(SdfReadError::Parse(format!(
            "Counts line too short: '{counts}' on line4"
        )));
    }

    let n_atoms = parse_unsigned_field(counts, 0, 3, 4)?;
    let n_bonds = parse_unsigned_field(counts, 3, 6, 4)?;

    if counts.len() > 35 {
        if counts.len() < 39 || counts.as_bytes().get(34) != Some(&b'V') {
            return Err(SdfReadError::Parse(
                "CTAB version string invalid at line 4".to_owned(),
            ));
        }
        let version = &counts[34..39];
        if version != "V2000" {
            return Err(SdfReadError::Parse(format!(
                "Unsupported CTAB version: '{version}' at line 4"
            )));
        }
    }

    let expected_atom_end = 4 + n_atoms;
    let expected_bond_end = expected_atom_end + n_bonds;
    if lines.len() < expected_bond_end {
        return Err(SdfReadError::Parse(
            "EOF hit while reading atoms or bonds".to_owned(),
        ));
    }

    let mut molecule = Molecule::new();
    let mut coords_2d = Vec::with_capacity(n_atoms);
    let mut coords_3d = Vec::with_capacity(n_atoms);
    let mut has_nonzero_z = false;
    for (offset, line) in lines[4..expected_atom_end].iter().enumerate() {
        let (mut atom, coord) = parse_v2000_atom_line(line, offset + 5)?;
        atom.index = offset;
        molecule.atoms_mut().push(atom);
        has_nonzero_z |= coord.z.abs() > 1e-12;
        coords_2d.push(DVec2::new(coord.x, coord.y));
        coords_3d.push(coord);
    }
    let coord_dim = source_dim.unwrap_or(if has_nonzero_z {
        CoordinateDimension::ThreeD
    } else {
        CoordinateDimension::TwoD
    });
    if matches!(coord_dim, CoordinateDimension::ThreeD) {
        molecule.conformers_3d_mut().push(coords_3d);
        molecule.set_source_coordinate_dim(Some(CoordinateDimension::ThreeD));
    } else {
        molecule.set_coords_2d(Some(coords_2d));
        molecule.set_source_coordinate_dim(Some(CoordinateDimension::TwoD));
    }

    for (offset, line) in lines[expected_atom_end..expected_bond_end]
        .iter()
        .enumerate()
    {
        let mut bond = parse_v2000_bond_line(line, offset + expected_atom_end + 1)?;
        bond.index = offset;
        if bond.is_aromatic {
            molecule.atoms_mut()[bond.begin_atom].is_aromatic = true;
            molecule.atoms_mut()[bond.end_atom].is_aromatic = true;
        }
        molecule.bonds_mut().push(bond);
    }

    let mut cursor = expected_bond_end;
    while cursor < lines.len() {
        let line = &lines[cursor];
        if line.starts_with("M  END") {
            finalize_parsed_stereochemistry(&mut molecule);
            return Ok((molecule, cursor + 1));
        }
        if line.starts_with("M  CHG") {
            apply_v2000_charge_line(&mut molecule, line, cursor + 1)?;
        } else if line.starts_with("M  ISO") {
            apply_v2000_isotope_line(&mut molecule, line, cursor + 1)?;
        } else if line.starts_with("M  RAD") {
            apply_v2000_radical_line(&mut molecule, line, cursor + 1)?;
        } else if !strip_rdkit(line).is_empty() {
            return Err(SdfReadError::Parse(format!(
                "Unsupported V2000 property line before M  END: '{line}'"
            )));
        }
        cursor += 1;
    }

    Err(SdfReadError::Parse(format!(
        "Problems encountered parsing Mol data, M  END missing around line {}",
        lines.len()
    )))
}

fn parse_v2000_atom_line(line: &str, line_number: usize) -> Result<(Atom, DVec3), SdfReadError> {
    if line.len() < 32 {
        return Err(SdfReadError::Parse(format!(
            "Atom line too short: '{line}' on line {line_number}"
        )));
    }

    let x = parse_double_field(line, 0, 10, line_number)?;
    let y = parse_double_field(line, 10, 20, line_number)?;
    let z = parse_double_field(line, 20, 30, line_number)?;

    let symbol = field(line, 31, 34).unwrap_or("").trim();
    let (atomic_num, isotope) = atomic_number_from_molfile_symbol(symbol).ok_or_else(|| {
        SdfReadError::Parse(format!(
            "Unsupported atom symbol '{symbol}' on line {line_number}"
        ))
    })?;

    let mass_diff = if line.len() >= 36 && field(line, 34, 36) != Some(" 0") {
        parse_int_field(line, 34, 36, line_number)?
    } else {
        0
    };
    if mass_diff != 0 {
        return Err(SdfReadError::Parse(format!(
            "V2000 atom isotope mass-difference fields are not implemented on line {line_number}"
        )));
    }

    let charge_code = if line.len() >= 39 && field(line, 36, 39) != Some("  0") {
        parse_int_field(line, 36, 39, line_number)?
    } else {
        0
    };
    let formal_charge = if charge_code == 0 { 0 } else { 4 - charge_code };

    let h_count = if line.len() >= 45 && field(line, 42, 45) != Some("  0") {
        parse_int_field(line, 42, 45, line_number)?
    } else {
        0
    };
    if h_count != 0 {
        return Err(SdfReadError::Parse(format!(
            "V2000 atom query hydrogen count fields are not implemented on line {line_number}"
        )));
    }

    let atom_map_num = if line.len() >= 63 && field(line, 60, 63) != Some("  0") {
        Some(parse_unsigned_field(line, 60, 63, line_number)? as u32)
    } else {
        None
    };
    let total_valence = if line.len() >= 51 && field(line, 48, 51) != Some("  0") {
        Some(parse_int_field(line, 48, 51, line_number)?)
    } else {
        None
    };
    let mut props = std::collections::BTreeMap::new();
    if let Some(total_valence) = total_valence {
        props.insert("_MolFileTotValence".to_owned(), total_valence.to_string());
    }

    Ok((
        Atom {
            index: 0,
            atomic_num,
            is_aromatic: false,
            formal_charge: formal_charge as i8,
            explicit_hydrogens: 0,
            no_implicit: false,
            num_radical_electrons: 0,
            chiral_tag: ChiralTag::Unspecified,
            isotope,
            atom_map_num,
            props,
            rdkit_cip_rank: None,
        },
        DVec3::new(x, y, z),
    ))
}

fn parse_v2000_bond_line(line: &str, line_number: usize) -> Result<Bond, SdfReadError> {
    if line.len() < 9 {
        return Err(SdfReadError::Parse(format!(
            "Bond line too short: '{line}' on line {line_number}"
        )));
    }

    let begin_atom = parse_unsigned_field(line, 0, 3, line_number)?
        .checked_sub(1)
        .ok_or_else(|| SdfReadError::Parse(format!("Bad bond atom index on line {line_number}")))?;
    let end_atom = parse_unsigned_field(line, 3, 6, line_number)?
        .checked_sub(1)
        .ok_or_else(|| SdfReadError::Parse(format!("Bad bond atom index on line {line_number}")))?;
    let bond_type = parse_unsigned_field(line, 6, 9, line_number)?;

    let order = match bond_type {
        0 | 8 => BondOrder::Null,
        1 => BondOrder::Single,
        2 => BondOrder::Double,
        3 => BondOrder::Triple,
        4 => BondOrder::Aromatic,
        9 => BondOrder::Dative,
        5..=7 => {
            return Err(SdfReadError::Parse(format!(
                "V2000 query bond type {bond_type} is not implemented on line {line_number}"
            )));
        }
        _ => BondOrder::Null,
    };

    let mut direction = BondDirection::None;
    let mut stereo = BondStereo::None;
    if line.len() >= 12 && field(line, 9, 12) != Some("  0") {
        match parse_unsigned_field(line, 9, 12, line_number)? {
            0 => {}
            1 => direction = BondDirection::EndUpRight,
            6 => direction = BondDirection::EndDownRight,
            3 => stereo = BondStereo::Any,
            4 => {
                return Err(SdfReadError::Parse(format!(
                    "V2000 unknown single-bond stereo is not implemented on line {line_number}"
                )));
            }
            value => {
                return Err(SdfReadError::Parse(format!(
                    "Unrecognized V2000 bond stereo value {value} on line {line_number}"
                )));
            }
        }
    }

    Ok(Bond {
        index: 0,
        begin_atom,
        end_atom,
        order,
        is_aromatic: matches!(order, BondOrder::Aromatic),
        direction,
        stereo,
        stereo_atoms: Vec::new(),
    })
}

fn apply_v2000_charge_line(
    molecule: &mut Molecule,
    line: &str,
    line_number: usize,
) -> Result<(), SdfReadError> {
    let tokens: Vec<_> = line.split_whitespace().collect();
    if tokens.len() < 3 || tokens[0] != "M" || tokens[1] != "CHG" {
        return Err(SdfReadError::Parse(format!(
            "Bad charge line: '{line}' on line {line_number}"
        )));
    }
    let count = tokens[2].parse::<usize>().map_err(|_| {
        SdfReadError::Parse(format!(
            "Cannot convert '{}' to int on line {line_number}",
            tokens[2]
        ))
    })?;
    if tokens.len() < 3 + count * 2 {
        return Err(SdfReadError::Parse(format!(
            "Bad charge line: '{line}' on line {line_number}"
        )));
    }
    for pair in tokens[3..3 + count * 2].chunks_exact(2) {
        let atom_idx = pair[0].parse::<usize>().map_err(|_| {
            SdfReadError::Parse(format!(
                "Cannot convert '{}' to int on line {line_number}",
                pair[0]
            ))
        })?;
        let charge = pair[1].parse::<i8>().map_err(|_| {
            SdfReadError::Parse(format!(
                "Cannot convert '{}' to int on line {line_number}",
                pair[1]
            ))
        })?;
        let zero_based = atom_idx.checked_sub(1).ok_or_else(|| {
            SdfReadError::Parse(format!("Bad charge atom index on line {line_number}"))
        })?;
        let atom = molecule.atoms_mut().get_mut(zero_based).ok_or_else(|| {
            SdfReadError::Parse(format!("Bad charge atom index on line {line_number}"))
        })?;
        atom.formal_charge = charge;
    }
    Ok(())
}

fn apply_v2000_isotope_line(
    molecule: &mut Molecule,
    line: &str,
    line_number: usize,
) -> Result<(), SdfReadError> {
    for (atom_idx, value) in parse_v2000_atom_int_pairs(line, line_number, "ISO")? {
        if value < 0 {
            return Err(SdfReadError::Parse(format!(
                "Negative isotope value on line {line_number}"
            )));
        }
        let atom = molecule
            .atoms_mut()
            .get_mut(atom_idx)
            .ok_or_else(|| SdfReadError::Parse(format!("Bad atom index on line {line_number}")))?;
        atom.isotope = Some(value as u16);
    }
    Ok(())
}

fn apply_v2000_radical_line(
    molecule: &mut Molecule,
    line: &str,
    line_number: usize,
) -> Result<(), SdfReadError> {
    for (atom_idx, value) in parse_v2000_atom_int_pairs(line, line_number, "RAD")? {
        let radicals = match value {
            0 => 0,
            1 => 2,
            2 => 1,
            3 => 2,
            _ => {
                return Err(SdfReadError::Parse(format!(
                    "Unrecognized radical value {value} on line {line_number}"
                )));
            }
        };
        let atom = molecule
            .atoms_mut()
            .get_mut(atom_idx)
            .ok_or_else(|| SdfReadError::Parse(format!("Bad atom index on line {line_number}")))?;
        atom.num_radical_electrons = radicals;
    }
    Ok(())
}

fn parse_v2000_atom_int_pairs(
    line: &str,
    line_number: usize,
    tag: &str,
) -> Result<Vec<(usize, i32)>, SdfReadError> {
    let expected = match tag {
        "ISO" => "M  ISO",
        "RAD" => "M  RAD",
        _ => unreachable!(),
    };
    if line.len() < 9 || field(line, 0, 6) != Some(expected) {
        return Err(SdfReadError::Parse(format!(
            "Bad {tag} line: '{line}' on line {line_number}"
        )));
    }
    let count = parse_usize_token(field(line, 6, 9).unwrap_or(""), line_number)?;
    let mut cursor = 9;
    let mut updates = Vec::with_capacity(count);
    for _ in 0..count {
        let atom_idx =
            parse_usize_token(field(line, cursor, cursor + 4).unwrap_or(""), line_number)?
                .checked_sub(1)
                .ok_or_else(|| {
                    SdfReadError::Parse(format!("Bad atom index on line {line_number}"))
                })?;
        cursor += 4;
        let value = parse_i32_token(field(line, cursor, cursor + 4).unwrap_or(""), line_number)?;
        cursor += 4;
        updates.push((atom_idx, value));
    }
    Ok(updates)
}

fn parse_v3000_mol_data_stream(
    lines: &[String],
    source_dim: Option<CoordinateDimension>,
) -> Result<(Molecule, usize), SdfReadError> {
    let mut cursor = 4;
    expect_v3000_line(lines, cursor, "BEGIN CTAB")?;
    cursor += 1;

    let counts = v3000_content(
        lines
            .get(cursor)
            .ok_or_else(|| SdfReadError::Parse("V3000 COUNTS line not found".to_owned()))?,
    )
    .ok_or_else(|| SdfReadError::Parse("V3000 COUNTS line not found".to_owned()))?;
    let count_tokens: Vec<_> = counts.split_whitespace().collect();
    if count_tokens.len() < 3 || count_tokens[0] != "COUNTS" {
        return Err(SdfReadError::Parse(
            "V3000 COUNTS line not found".to_owned(),
        ));
    }
    let n_atoms = parse_usize_token(count_tokens[1], cursor + 1)?;
    let n_bonds = parse_usize_token(count_tokens[2], cursor + 1)?;
    cursor += 1;

    expect_v3000_line(lines, cursor, "BEGIN ATOM")?;
    cursor += 1;
    let mut molecule = Molecule::new();
    let mut coords_2d = Vec::with_capacity(n_atoms);
    let mut coords_3d = Vec::with_capacity(n_atoms);
    let mut has_nonzero_z = false;
    for _ in 0..n_atoms {
        let (mut atom, coord) = parse_v3000_atom_line(
            v3000_content(lines.get(cursor).ok_or_else(|| {
                SdfReadError::Parse("EOF hit while reading V3000 atoms".to_owned())
            })?)
            .ok_or_else(|| SdfReadError::Parse("Bad V3000 atom line".to_owned()))?,
            cursor + 1,
        )?;
        atom.index = molecule.atoms().len();
        molecule.atoms_mut().push(atom);
        has_nonzero_z |= coord.z.abs() > 1e-12;
        coords_2d.push(DVec2::new(coord.x, coord.y));
        coords_3d.push(coord);
        cursor += 1;
    }
    let coord_dim = source_dim.unwrap_or(if has_nonzero_z {
        CoordinateDimension::ThreeD
    } else {
        CoordinateDimension::TwoD
    });
    if matches!(coord_dim, CoordinateDimension::ThreeD) {
        molecule.conformers_3d_mut().push(coords_3d);
        molecule.set_source_coordinate_dim(Some(CoordinateDimension::ThreeD));
    } else {
        molecule.set_coords_2d(Some(coords_2d));
        molecule.set_source_coordinate_dim(Some(CoordinateDimension::TwoD));
    }
    expect_v3000_line(lines, cursor, "END ATOM")?;
    cursor += 1;

    let has_bond_block = lines
        .get(cursor)
        .and_then(|line| v3000_content(line))
        .is_some_and(|content| content == "BEGIN BOND");
    if n_bonds > 0 || has_bond_block {
        expect_v3000_line(lines, cursor, "BEGIN BOND")?;
        cursor += 1;
        for _ in 0..n_bonds {
            let mut bond = parse_v3000_bond_line(
                v3000_content(lines.get(cursor).ok_or_else(|| {
                    SdfReadError::Parse("EOF hit while reading V3000 bonds".to_owned())
                })?)
                .ok_or_else(|| SdfReadError::Parse("Bad V3000 bond line".to_owned()))?,
                cursor + 1,
            )?;
            bond.index = molecule.bonds().len();
            if bond.is_aromatic {
                molecule.atoms_mut()[bond.begin_atom].is_aromatic = true;
                molecule.atoms_mut()[bond.end_atom].is_aromatic = true;
            }
            molecule.bonds_mut().push(bond);
            cursor += 1;
        }
        expect_v3000_line(lines, cursor, "END BOND")?;
        cursor += 1;
    }

    expect_v3000_line(lines, cursor, "END CTAB")?;
    cursor += 1;
    if lines
        .get(cursor)
        .is_some_and(|line| line.starts_with("M  END"))
    {
        finalize_parsed_stereochemistry(&mut molecule);
        Ok((molecule, cursor + 1))
    } else {
        Err(SdfReadError::Parse(format!(
            "Problems encountered parsing Mol data, M  END missing around line {}",
            cursor + 1
        )))
    }
}

fn finalize_parsed_stereochemistry(molecule: &mut Molecule) {
    let explicit_valences = molfile_explicit_valences(molecule);
    apply_molfile_total_valence(molecule, &explicit_valences);
    cache_molfile_explicit_valence_for_aromaticity(molecule, &explicit_valences);
    crate::smiles::sanitize_molfile_aromaticity(molecule);
    apply_aromatic_n_h_from_molfile_valence(molecule);
    let has_bond_stereo_markers = molecule
        .bonds()
        .iter()
        .any(|bond| !matches!(bond.direction, BondDirection::None));
    if matches!(
        molecule.source_coordinate_dim(),
        Some(CoordinateDimension::ThreeD)
    ) {
        crate::stereo::assign_chiral_types_from_3d_rdkit_subset(molecule);
    } else if has_bond_stereo_markers {
        crate::stereo::assign_chiral_types_from_bond_dirs_rdkit_subset(molecule);
    }
    clear_single_bond_dir_flags(molecule);
    set_double_bond_neighbor_directions_from_coords_rdkit_subset(molecule);
    crate::smiles::assign_double_bond_stereo_from_directions(molecule);
    crate::smiles::cleanup_nonstereo_double_bond_dirs(molecule);
    crate::stereo::cache_rdkit_legacy_cip_ranks(molecule);
}

fn molfile_explicit_valences(molecule: &Molecule) -> Vec<i32> {
    let mut explicit_valences = vec![0; molecule.atoms().len()];
    for bond in molecule.bonds() {
        let contribution = match bond.order {
            BondOrder::Single | BondOrder::Dative => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Quadruple => 4,
            BondOrder::Aromatic => 1,
            BondOrder::Null => 0,
        };
        explicit_valences[bond.begin_atom] += contribution;
        explicit_valences[bond.end_atom] += contribution;
    }
    explicit_valences
}

fn cache_molfile_explicit_valence_for_aromaticity(
    molecule: &mut Molecule,
    explicit_valences: &[i32],
) {
    for atom in molecule.atoms_mut() {
        atom.props.insert(
            "_MolFileExplicitValence".to_owned(),
            explicit_valences[atom.index].to_string(),
        );
    }
}

fn apply_aromatic_n_h_from_molfile_valence(molecule: &mut Molecule) {
    for atom in molecule.atoms_mut() {
        let explicit_valence = atom
            .props
            .remove("_MolFileExplicitValence")
            .and_then(|value| value.parse::<i32>().ok());
        if atom.atomic_num == 7
            && atom.formal_charge == 0
            && atom.is_aromatic
            && atom.explicit_hydrogens == 0
            && explicit_valence == Some(2)
        {
            atom.explicit_hydrogens = 1;
            atom.no_implicit = true;
        }
    }
}

fn apply_molfile_total_valence(molecule: &mut Molecule, explicit_valences: &[i32]) {
    for atom in molecule.atoms_mut() {
        let Some(value) = atom
            .props
            .remove("_MolFileTotValence")
            .and_then(|value| value.parse::<i32>().ok())
        else {
            continue;
        };
        if value == 0 {
            continue;
        }
        atom.no_implicit = true;
        if value == 15 || value == -1 || explicit_valences[atom.index] > value {
            atom.explicit_hydrogens = 0;
        } else {
            atom.explicit_hydrogens = (value - explicit_valences[atom.index]) as u8;
        }
    }
}

fn clear_single_bond_dir_flags(molecule: &mut Molecule) {
    for bond in molecule.bonds_mut() {
        if matches!(bond.order, BondOrder::Single) {
            bond.direction = BondDirection::None;
        }
    }
}

fn set_double_bond_neighbor_directions_from_coords_rdkit_subset(molecule: &mut Molecule) {
    let coords = if let Some(coords_2d) = molecule.coords_2d() {
        coords_2d
            .iter()
            .map(|coord| DVec3::new(coord.x, coord.y, 0.0))
            .collect::<Vec<_>>()
    } else if let Some(coords_3d) = molecule.coords_3d() {
        coords_3d.to_vec()
    } else {
        return;
    };
    if coords.len() != molecule.atoms().len() {
        return;
    }
    let adjacency = AdjacencyList::from_topology(molecule.atoms().len(), molecule.bonds());
    let mut cycle_cache = BondCycleCache::new(molecule.bonds().len());

    let mut single_bond_counts = vec![0usize; molecule.bonds().len()];
    let mut double_bond_neighbors = vec![Vec::<usize>::new(); molecule.bonds().len()];
    let mut single_bond_neighbors = vec![Vec::<usize>::new(); molecule.bonds().len()];
    let mut needs_dir = vec![false; molecule.bonds().len()];
    let mut bonds_in_play = Vec::<usize>::new();

    for bond in molecule.bonds() {
        if !is_double_bond_candidate_for_stereo(molecule, &adjacency, &mut cycle_cache, bond.index)
        {
            continue;
        }
        for atom_idx in [bond.begin_atom, bond.end_atom] {
            for neighbor in adjacency.neighbors_of(atom_idx) {
                let neighbor_idx = neighbor.bond_index;
                let neighbor = &molecule.bonds()[neighbor_idx];
                if matches!(neighbor.order, BondOrder::Single | BondOrder::Aromatic) {
                    single_bond_counts[neighbor_idx] += 1;
                    needs_dir[bond.index] = true;
                    if matches!(
                        neighbor.direction,
                        BondDirection::None
                            | BondDirection::EndDownRight
                            | BondDirection::EndUpRight
                    ) {
                        needs_dir[neighbor_idx] = true;
                        double_bond_neighbors[bond.index].push(neighbor_idx);
                        if !single_bond_neighbors[neighbor_idx].contains(&bond.index) {
                            single_bond_neighbors[neighbor_idx].push(bond.index);
                        }
                    }
                }
            }
        }
        bonds_in_play.push(bond.index);
    }

    bonds_in_play.sort_by_key(|bond_idx| {
        let mut count = double_bond_neighbors[*bond_idx]
            .iter()
            .map(|neighbor_idx| single_bond_counts[*neighbor_idx])
            .sum::<usize>();
        if cycle_cache
            .min_cycle_size_for_bond(molecule, &adjacency, *bond_idx)
            .is_none()
        {
            count *= 10;
        }
        count
    });
    for bond_idx in bonds_in_play.into_iter().rev() {
        update_double_bond_neighbors(
            molecule,
            &coords,
            &adjacency,
            bond_idx,
            &mut needs_dir,
            &single_bond_counts,
            &single_bond_neighbors,
        );
    }
}

fn update_double_bond_neighbors(
    molecule: &mut Molecule,
    coords: &[DVec3],
    adjacency: &AdjacencyList,
    double_bond_idx: usize,
    needs_dir: &mut [bool],
    single_bond_counts: &[usize],
    single_bond_neighbors: &[Vec<usize>],
) {
    if !needs_dir[double_bond_idx] {
        return;
    }
    needs_dir[double_bond_idx] = false;
    let double_bond = molecule.bonds()[double_bond_idx].clone();
    let Some((Some(mut bond1_idx), mut obond1_idx)) = controlling_bond_from_atom(
        molecule,
        needs_dir,
        single_bond_counts,
        double_bond_idx,
        double_bond.begin_atom,
        adjacency,
    ) else {
        return;
    };
    let Some((Some(mut bond2_idx), mut obond2_idx)) = controlling_bond_from_atom(
        molecule,
        needs_dir,
        single_bond_counts,
        double_bond_idx,
        double_bond.end_atom,
        adjacency,
    ) else {
        return;
    };

    let mut bond1_other = other_atom_index(&molecule.bonds()[bond1_idx], double_bond.begin_atom);
    let mut bond2_other = other_atom_index(&molecule.bonds()[bond2_idx], double_bond.end_atom);
    let begin_p = coords[double_bond.begin_atom];
    let end_p = coords[double_bond.end_atom];
    let mut bond1_p = coords[bond1_other];
    let mut bond2_p = coords[bond2_other];
    let mut linear = false;
    let mut p1 = bond1_p - begin_p;
    let p2 = end_p - begin_p;
    if is_linear_arrangement(p1, p2) {
        if let Some(obond1) = obond1_idx {
            obond1_idx = Some(bond1_idx);
            bond1_idx = obond1;
            bond1_other = other_atom_index(&molecule.bonds()[bond1_idx], double_bond.begin_atom);
            bond1_p = coords[bond1_other];
            p1 = bond1_p - begin_p;
            if is_linear_arrangement(p1, p2) {
                linear = true;
            }
        } else {
            linear = true;
        }
    }
    if !linear {
        p1 = bond2_p - end_p;
        let p2 = begin_p - end_p;
        if is_linear_arrangement(p1, p2) {
            if let Some(obond2) = obond2_idx {
                obond2_idx = Some(bond2_idx);
                bond2_idx = obond2;
                bond2_other = other_atom_index(&molecule.bonds()[bond2_idx], double_bond.end_atom);
                bond2_p = coords[bond2_other];
                p1 = bond2_p - begin_p;
                if is_linear_arrangement(p1, p2) {
                    linear = true;
                }
            } else {
                linear = true;
            }
        }
    }
    if linear {
        molecule.bonds_mut()[double_bond_idx].stereo = BondStereo::Any;
        return;
    }

    let same_torsion_dir =
        compute_dihedral_angle(bond1_p, begin_p, end_p, bond2_p) >= std::f64::consts::FRAC_PI_2;
    let mut reverse_bond_dir = same_torsion_dir;
    let atom1 = double_bond.begin_atom;
    let atom2 = double_bond.end_atom;
    let mut followup_bonds = Vec::<usize>::new();
    if needs_dir[bond1_idx] {
        for bidx in &single_bond_neighbors[bond1_idx] {
            if needs_dir[*bidx] {
                followup_bonds.push(*bidx);
            }
        }
    }
    if needs_dir[bond2_idx] {
        for bidx in &single_bond_neighbors[bond2_idx] {
            if needs_dir[*bidx] {
                followup_bonds.push(*bidx);
            }
        }
    }
    if !needs_dir[bond1_idx] {
        if needs_dir[bond2_idx] {
            if molecule.bonds()[bond1_idx].begin_atom != atom1 {
                reverse_bond_dir = !reverse_bond_dir;
            }
            let dir = molecule.bonds()[bond1_idx].direction;
            molecule.bonds_mut()[bond2_idx].direction = bond_dir_relative_to_atom(
                &molecule.bonds()[bond2_idx],
                atom2,
                dir,
                reverse_bond_dir,
            );
        }
    } else if !needs_dir[bond2_idx] {
        if molecule.bonds()[bond2_idx].begin_atom != atom2 {
            reverse_bond_dir = !reverse_bond_dir;
        }
        let dir = molecule.bonds()[bond2_idx].direction;
        molecule.bonds_mut()[bond1_idx].direction =
            bond_dir_relative_to_atom(&molecule.bonds()[bond1_idx], atom1, dir, reverse_bond_dir);
    } else {
        molecule.bonds_mut()[bond1_idx].direction = bond_dir_relative_to_atom(
            &molecule.bonds()[bond1_idx],
            atom1,
            BondDirection::EndDownRight,
            false,
        );
        molecule.bonds_mut()[bond2_idx].direction = bond_dir_relative_to_atom(
            &molecule.bonds()[bond2_idx],
            atom2,
            BondDirection::EndDownRight,
            reverse_bond_dir,
        );
    }
    needs_dir[bond1_idx] = false;
    needs_dir[bond2_idx] = false;
    if let Some(obond1_idx) = obond1_idx
        && needs_dir[obond1_idx]
    {
        let dir = molecule.bonds()[bond1_idx].direction;
        let reverse = molecule.bonds()[bond1_idx].begin_atom == atom1;
        molecule.bonds_mut()[obond1_idx].direction =
            bond_dir_relative_to_atom(&molecule.bonds()[obond1_idx], atom1, dir, reverse);
        needs_dir[obond1_idx] = false;
    }
    if let Some(obond2_idx) = obond2_idx
        && needs_dir[obond2_idx]
    {
        let dir = molecule.bonds()[bond2_idx].direction;
        let reverse = molecule.bonds()[bond2_idx].begin_atom == atom2;
        molecule.bonds_mut()[obond2_idx].direction =
            bond_dir_relative_to_atom(&molecule.bonds()[obond2_idx], atom2, dir, reverse);
        needs_dir[obond2_idx] = false;
    }
    for followup in followup_bonds {
        update_double_bond_neighbors(
            molecule,
            coords,
            adjacency,
            followup,
            needs_dir,
            single_bond_counts,
            single_bond_neighbors,
        );
    }
}

fn controlling_bond_from_atom(
    molecule: &Molecule,
    needs_dir: &[bool],
    single_bond_counts: &[usize],
    double_bond_idx: usize,
    atom_idx: usize,
    adjacency: &AdjacencyList,
) -> Option<(Option<usize>, Option<usize>)> {
    let mut bond = None;
    let mut obond = None;
    for neighbor in adjacency.neighbors_of(atom_idx) {
        let bond_idx = neighbor.bond_index;
        if bond_idx == double_bond_idx {
            continue;
        }
        let candidate = &molecule.bonds()[bond_idx];
        if matches!(candidate.order, BondOrder::Single | BondOrder::Aromatic)
            && matches!(
                candidate.direction,
                BondDirection::None | BondDirection::EndDownRight | BondDirection::EndUpRight
            )
        {
            if let Some(current) = bond {
                if needs_dir[bond_idx] {
                    if single_bond_counts[bond_idx] > single_bond_counts[current] {
                        obond = bond;
                        bond = Some(bond_idx);
                    } else {
                        obond = Some(bond_idx);
                    }
                } else {
                    obond = bond;
                    bond = Some(bond_idx);
                }
            } else {
                bond = Some(bond_idx);
            }
        }
    }
    Some((bond, obond))
}

fn other_atom_index(bond: &Bond, atom_idx: usize) -> usize {
    if bond.begin_atom == atom_idx {
        bond.end_atom
    } else {
        bond.begin_atom
    }
}

fn opposite_bond_direction(direction: BondDirection) -> BondDirection {
    match direction {
        BondDirection::EndUpRight => BondDirection::EndDownRight,
        BondDirection::EndDownRight => BondDirection::EndUpRight,
        BondDirection::None => BondDirection::None,
    }
}

fn bond_dir_relative_to_atom(
    bond: &Bond,
    atom_idx: usize,
    direction: BondDirection,
    mut reverse: bool,
) -> BondDirection {
    if bond.begin_atom != atom_idx {
        reverse = !reverse;
    }
    if reverse {
        opposite_bond_direction(direction)
    } else {
        direction
    }
}

struct BondCycleCache {
    min_cycle_sizes: Vec<Option<Option<usize>>>,
}

impl BondCycleCache {
    fn new(bond_count: usize) -> Self {
        Self {
            min_cycle_sizes: vec![None; bond_count],
        }
    }

    fn min_cycle_size_for_bond(
        &mut self,
        molecule: &Molecule,
        adjacency: &AdjacencyList,
        bond_idx: usize,
    ) -> Option<usize> {
        if let Some(cached) = self.min_cycle_sizes[bond_idx] {
            return cached;
        }
        let result = min_cycle_size_for_bond(molecule, adjacency, bond_idx);
        self.min_cycle_sizes[bond_idx] = Some(result);
        result
    }
}

fn min_cycle_size_for_bond(
    molecule: &Molecule,
    adjacency: &AdjacencyList,
    bond_idx: usize,
) -> Option<usize> {
    let bond = molecule.bonds().get(bond_idx)?;
    let mut visited = vec![false; molecule.atoms().len()];
    let mut queue = std::collections::VecDeque::new();
    visited[bond.begin_atom] = true;
    queue.push_back((bond.begin_atom, 0usize));
    while let Some((atom_idx, distance)) = queue.pop_front() {
        for neighbor in adjacency.neighbors_of(atom_idx) {
            if neighbor.bond_index == bond_idx {
                continue;
            }
            let next = neighbor.atom_index;
            if next == bond.end_atom {
                return Some(distance + 2);
            }
            if !visited[next] {
                visited[next] = true;
                queue.push_back((next, distance + 1));
            }
        }
    }
    None
}

fn is_double_bond_candidate_for_stereo(
    molecule: &Molecule,
    adjacency: &AdjacencyList,
    cycle_cache: &mut BondCycleCache,
    bond_idx: usize,
) -> bool {
    let bond = &molecule.bonds()[bond_idx];
    matches!(bond.order, BondOrder::Double)
        && !matches!(bond.stereo, BondStereo::Any)
        && cycle_cache
            .min_cycle_size_for_bond(molecule, adjacency, bond_idx)
            .is_none_or(|size| size >= 8)
        && adjacency.neighbors_of(bond.begin_atom).len() > 1
        && adjacency.neighbors_of(bond.end_atom).len() > 1
}

fn is_linear_arrangement(v1: DVec3, v2: DVec3) -> bool {
    let length_sq = v1.length_squared() * v2.length_squared();
    if length_sq < 1.0e-6 {
        return true;
    }
    let cos178 = -0.999388_f64;
    v1.dot(v2) < cos178 * length_sq.sqrt()
}

fn compute_dihedral_angle(pt1: DVec3, pt2: DVec3, pt3: DVec3, pt4: DVec3) -> f64 {
    let beg_end_vec = pt3 - pt2;
    let beg_nbr_vec = pt1 - pt2;
    let crs1 = beg_nbr_vec.cross(beg_end_vec);
    let end_nbr_vec = pt4 - pt3;
    let crs2 = end_nbr_vec.cross(beg_end_vec);
    let denom = crs1.length() * crs2.length();
    if denom == 0.0 {
        return 0.0;
    }
    (crs1.dot(crs2) / denom).clamp(-1.0, 1.0).acos()
}

fn parse_v3000_atom_line(line: &str, line_number: usize) -> Result<(Atom, DVec3), SdfReadError> {
    let tokens: Vec<_> = line.split_whitespace().collect();
    if tokens.len() < 6 {
        return Err(SdfReadError::Parse(format!(
            "Bad atom line : '{line}' on line {line_number}"
        )));
    }
    let symbol = tokens[1];
    let (atomic_num, mut isotope) = atomic_number_from_molfile_symbol(symbol).ok_or_else(|| {
        SdfReadError::Parse(format!(
            "Unsupported atom symbol '{symbol}' on line {line_number}"
        ))
    })?;
    let x = parse_f64_token(tokens[2], line_number)?;
    let y = parse_f64_token(tokens[3], line_number)?;
    let z = parse_f64_token(tokens[4], line_number)?;
    let atom_map_num = match parse_i32_token(tokens[5], line_number)? {
        value if value > 0 => Some(value as u32),
        _ => None,
    };

    let mut formal_charge = 0;
    let mut radicals = 0;
    let mut props = std::collections::BTreeMap::new();
    for token in tokens.iter().skip(6) {
        let Some((prop, value)) = token.split_once('=') else {
            return Err(SdfReadError::Parse(format!(
                "Invalid atom property: '{token}' for atom on line {line_number}"
            )));
        };
        match prop {
            "CHG" => formal_charge = parse_i32_token(value, line_number)? as i8,
            "MASS" => {
                let value = parse_i32_token(value, line_number)?;
                if value < 0 {
                    return Err(SdfReadError::Parse(format!(
                        "Bad value for MASS :{value} for atom on line {line_number}"
                    )));
                }
                isotope = Some(value as u16);
            }
            "RAD" => {
                radicals = match parse_i32_token(value, line_number)? {
                    0 => 0,
                    1 => 2,
                    2 => 1,
                    3 => 2,
                    other => {
                        return Err(SdfReadError::Parse(format!(
                            "Unrecognized RAD value {other} for atom on line {line_number}"
                        )));
                    }
                };
            }
            "VAL" => {
                let total_valence = parse_i32_token(value, line_number)?;
                if total_valence != 0 {
                    props.insert("_MolFileTotValence".to_owned(), total_valence.to_string());
                }
            }
            _ => {
                props.insert(prop.to_owned(), value.to_owned());
            }
        }
    }

    Ok((
        Atom {
            index: 0,
            atomic_num,
            is_aromatic: false,
            formal_charge,
            explicit_hydrogens: 0,
            no_implicit: false,
            num_radical_electrons: radicals,
            chiral_tag: ChiralTag::Unspecified,
            isotope,
            atom_map_num,
            props,
            rdkit_cip_rank: None,
        },
        DVec3::new(x, y, z),
    ))
}

fn parse_v3000_bond_line(line: &str, line_number: usize) -> Result<Bond, SdfReadError> {
    let tokens: Vec<_> = line.split_whitespace().collect();
    if tokens.len() < 4 {
        return Err(SdfReadError::Parse(format!(
            "bond line {line_number} is too short"
        )));
    }
    let bond_type = parse_usize_token(tokens[1], line_number)?;
    let begin_atom = parse_usize_token(tokens[2], line_number)?
        .checked_sub(1)
        .ok_or_else(|| SdfReadError::Parse(format!("Bad bond atom index on line {line_number}")))?;
    let end_atom = parse_usize_token(tokens[3], line_number)?
        .checked_sub(1)
        .ok_or_else(|| SdfReadError::Parse(format!("Bad bond atom index on line {line_number}")))?;
    let order = match bond_type {
        0 | 8 => BondOrder::Null,
        1 => BondOrder::Single,
        2 => BondOrder::Double,
        3 => BondOrder::Triple,
        4 => BondOrder::Aromatic,
        9 => BondOrder::Dative,
        5..=7 | 10 => {
            return Err(SdfReadError::Parse(format!(
                "V3000 bond type {bond_type} is not implemented on line {line_number}"
            )));
        }
        _ => BondOrder::Null,
    };
    let mut direction = BondDirection::None;
    let mut stereo = BondStereo::None;
    for token in tokens.iter().skip(4) {
        let Some((prop, value)) = token.split_once('=') else {
            return Err(SdfReadError::Parse(format!(
                "bad bond property '{token}' on line {line_number}"
            )));
        };
        match prop {
            "CFG" => match parse_usize_token(value, line_number)? {
                0 => {}
                1 => direction = BondDirection::EndUpRight,
                2 => {
                    if bond_type == 2 {
                        stereo = BondStereo::Any;
                    } else {
                        return Err(SdfReadError::Parse(format!(
                            "V3000 unknown single-bond CFG is not implemented on line {line_number}"
                        )));
                    }
                }
                3 => direction = BondDirection::EndDownRight,
                cfg => {
                    return Err(SdfReadError::Parse(format!(
                        "bad bond CFG {cfg}' on line {line_number}"
                    )));
                }
            },
            _ => {
                return Err(SdfReadError::Parse(format!(
                    "Unsupported V3000 bond property '{prop}' on line {line_number}"
                )));
            }
        }
    }
    Ok(Bond {
        index: 0,
        begin_atom,
        end_atom,
        order,
        is_aromatic: matches!(order, BondOrder::Aromatic),
        direction,
        stereo,
        stereo_atoms: Vec::new(),
    })
}

fn expect_v3000_line(lines: &[String], index: usize, expected: &str) -> Result<(), SdfReadError> {
    let line = lines
        .get(index)
        .ok_or_else(|| SdfReadError::Parse(format!("{expected} line not found")))?;
    let content = v3000_content(line)
        .ok_or_else(|| SdfReadError::Parse(format!("{expected} line not found")))?;
    if content == expected {
        Ok(())
    } else {
        Err(SdfReadError::Parse(format!("{expected} line not found")))
    }
}

fn v3000_content(line: &str) -> Option<&str> {
    line.strip_prefix("M  V30 ").map(strip_rdkit)
}

fn infer_header_coordinate_dimension(line: Option<&str>) -> Option<CoordinateDimension> {
    let line = line?;
    let trimmed = line.trim();
    if trimmed.ends_with("3D") || trimmed.contains(" 3D") {
        Some(CoordinateDimension::ThreeD)
    } else if trimmed.ends_with("2D") || trimmed.contains(" 2D") {
        Some(CoordinateDimension::TwoD)
    } else {
        None
    }
}

fn resolve_coordinate_dimension(
    coordinate_mode: SdfCoordinateMode,
    inferred: Option<CoordinateDimension>,
) -> Option<CoordinateDimension> {
    match coordinate_mode {
        SdfCoordinateMode::Auto => inferred,
        SdfCoordinateMode::Force2D => Some(CoordinateDimension::TwoD),
        SdfCoordinateMode::Force3D => Some(CoordinateDimension::ThreeD),
    }
}

fn parse_sdf_data_fields(lines: &[String]) -> Result<Vec<(String, String)>, SdfReadError> {
    let mut fields = Vec::new();
    let mut cursor = 0;

    while cursor < lines.len() {
        let temp = strip_rdkit(&lines[cursor]);
        if temp.is_empty() {
            cursor += 1;
            continue;
        }
        if !temp.starts_with('>') {
            return Err(SdfReadError::Parse(
                "Problems encountered parsing data fields".to_owned(),
            ));
        }

        let header = &temp[1..];
        let Some(label_start) = header.find('<') else {
            cursor = skip_until_blank(lines, cursor + 1);
            continue;
        };
        let Some(label_end) = header.rfind('>') else {
            cursor = skip_until_blank(lines, cursor + 1);
            continue;
        };
        if label_end == label_start + 1 {
            cursor = skip_until_blank(lines, cursor + 1);
            continue;
        }

        let label = header[label_start + 1..label_end].to_owned();
        cursor += 1;
        let mut value = String::new();
        let mut n_lines = 0;
        while cursor < lines.len() {
            let raw = &lines[cursor];
            let stripped = strip_rdkit(raw);
            if stripped.is_empty() && !raw.starts_with(' ') && !raw.starts_with('\t') {
                break;
            }
            if n_lines > 0 {
                value.push('\n');
            }
            value.push_str(raw.strip_suffix('\r').unwrap_or(raw));
            n_lines += 1;
            cursor += 1;
        }
        fields.push((label, value));
    }

    Ok(fields)
}

fn skip_until_blank(lines: &[String], mut cursor: usize) -> usize {
    while cursor < lines.len() && !strip_rdkit(&lines[cursor]).is_empty() {
        cursor += 1;
    }
    cursor
}

fn parse_unsigned_field(
    line: &str,
    start: usize,
    end: usize,
    line_number: usize,
) -> Result<usize, SdfReadError> {
    let token = field(line, start, end).unwrap_or("");
    token.trim().parse::<usize>().map_err(|_| {
        SdfReadError::Parse(format!(
            "Cannot convert '{token}' to unsigned int on line {line_number}"
        ))
    })
}

fn parse_int_field(
    line: &str,
    start: usize,
    end: usize,
    line_number: usize,
) -> Result<i32, SdfReadError> {
    let token = field(line, start, end).unwrap_or("");
    token.trim().parse::<i32>().map_err(|_| {
        SdfReadError::Parse(format!(
            "Cannot convert '{token}' to int on line {line_number}"
        ))
    })
}

fn parse_double_field(
    line: &str,
    start: usize,
    end: usize,
    line_number: usize,
) -> Result<f64, SdfReadError> {
    let token = field(line, start, end).unwrap_or("");
    token.trim().parse::<f64>().map_err(|_| {
        SdfReadError::Parse(format!("Cannot process coordinates on line {line_number}"))
    })
}

fn parse_usize_token(token: &str, line_number: usize) -> Result<usize, SdfReadError> {
    token.trim().parse::<usize>().map_err(|_| {
        SdfReadError::Parse(format!(
            "Cannot convert '{token}' to unsigned int on line {line_number}"
        ))
    })
}

fn parse_i32_token(token: &str, line_number: usize) -> Result<i32, SdfReadError> {
    token.trim().parse::<i32>().map_err(|_| {
        SdfReadError::Parse(format!(
            "Cannot convert '{token}' to int on line {line_number}"
        ))
    })
}

fn parse_f64_token(token: &str, line_number: usize) -> Result<f64, SdfReadError> {
    token.trim().parse::<f64>().map_err(|_| {
        SdfReadError::Parse(format!("Cannot process coordinates on line {line_number}"))
    })
}

fn field(line: &str, start: usize, end: usize) -> Option<&str> {
    line.get(start..end)
}

fn trim_line_ending(line: &str) -> &str {
    line.strip_suffix('\n')
        .unwrap_or(line)
        .strip_suffix('\r')
        .unwrap_or(line.strip_suffix('\n').unwrap_or(line))
}

fn strip_rdkit(line: &str) -> &str {
    line.trim_matches([' ', '\t', '\r', '\n'])
}

fn atomic_number_from_molfile_symbol(symbol: &str) -> Option<(u8, Option<u16>)> {
    match symbol {
        "*" | "A" | "Q" | "L" | "LP" | "R" | "R#" => Some((0, None)),
        "D" => Some((1, Some(2))),
        "T" => Some((1, Some(3))),
        _ => crate::periodic_table::atomic_number(symbol).map(|atomic_num| (atomic_num, None)),
    }
}
