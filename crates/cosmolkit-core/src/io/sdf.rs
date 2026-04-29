use std::io::BufRead;

use crate::{Atom, Bond, BondDirection, BondOrder, BondStereo, ChiralTag, Molecule};
use glam::DVec2;

/// One record extracted from an SDF stream.
#[derive(Debug, Clone, PartialEq)]
pub struct SdfRecord {
    pub molecule: Molecule,
    pub data_fields: Vec<(String, String)>,
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

/// Streaming-capable SDF reader.
pub struct SdfReader<R> {
    reader: R,
    line_number: usize,
    at_end: bool,
}

impl<R: BufRead> SdfReader<R> {
    /// Create a reader from any buffered input stream.
    #[must_use]
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            line_number: 0,
            at_end: false,
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

        let (mut molecule, data_start) = parse_mol_data_stream(&lines)?;
        let data_fields = parse_sdf_data_fields(&lines[data_start..])?;
        molecule.rebuild_adjacency();

        Ok(Some(SdfRecord {
            molecule,
            data_fields,
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

fn parse_mol_data_stream(lines: &[String]) -> Result<(Molecule, usize), SdfReadError> {
    let counts = lines.get(3).ok_or_else(|| {
        SdfReadError::Parse(format!(
            "Counts line too short: '{}' on line4",
            lines.get(3).map_or("", String::as_str)
        ))
    })?;
    if counts.len() > 35 && counts.len() >= 39 && counts.as_bytes().get(34) == Some(&b'V') {
        match &counts[34..39] {
            "V2000" => parse_v2000_mol_data_stream(lines),
            "V3000" => parse_v3000_mol_data_stream(lines),
            version => Err(SdfReadError::Parse(format!(
                "Unsupported CTAB version: '{version}' at line 4"
            ))),
        }
    } else {
        parse_v2000_mol_data_stream(lines)
    }
}

fn parse_v2000_mol_data_stream(lines: &[String]) -> Result<(Molecule, usize), SdfReadError> {
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
    let mut coords = Vec::with_capacity(n_atoms);
    for (offset, line) in lines[4..expected_atom_end].iter().enumerate() {
        let (atom, coord) = parse_v2000_atom_line(line, offset + 5)?;
        molecule.add_atom(atom);
        coords.push(coord);
    }
    molecule.coords_2d = Some(coords);

    for (offset, line) in lines[expected_atom_end..expected_bond_end]
        .iter()
        .enumerate()
    {
        let bond = parse_v2000_bond_line(line, offset + expected_atom_end + 1)?;
        if matches!(bond.order, BondOrder::Aromatic) {
            molecule.atoms[bond.begin_atom].is_aromatic = true;
            molecule.atoms[bond.end_atom].is_aromatic = true;
        }
        molecule.add_bond(bond);
    }

    let mut cursor = expected_bond_end;
    while cursor < lines.len() {
        let line = &lines[cursor];
        if line.starts_with("M  END") {
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

fn parse_v2000_atom_line(line: &str, line_number: usize) -> Result<(Atom, DVec2), SdfReadError> {
    if line.len() < 32 {
        return Err(SdfReadError::Parse(format!(
            "Atom line too short: '{line}' on line {line_number}"
        )));
    }

    let x = parse_double_field(line, 0, 10, line_number)?;
    let y = parse_double_field(line, 10, 20, line_number)?;
    parse_double_field(line, 20, 30, line_number)?;

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
        },
        DVec2::new(x, y),
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
        let atom = molecule.atoms.get_mut(zero_based).ok_or_else(|| {
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
            .atoms
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
            .atoms
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

fn parse_v3000_mol_data_stream(lines: &[String]) -> Result<(Molecule, usize), SdfReadError> {
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
    let mut coords = Vec::with_capacity(n_atoms);
    for _ in 0..n_atoms {
        let (atom, coord) = parse_v3000_atom_line(
            v3000_content(lines.get(cursor).ok_or_else(|| {
                SdfReadError::Parse("EOF hit while reading V3000 atoms".to_owned())
            })?)
            .ok_or_else(|| SdfReadError::Parse("Bad V3000 atom line".to_owned()))?,
            cursor + 1,
        )?;
        molecule.add_atom(atom);
        coords.push(coord);
        cursor += 1;
    }
    molecule.coords_2d = Some(coords);
    expect_v3000_line(lines, cursor, "END ATOM")?;
    cursor += 1;

    expect_v3000_line(lines, cursor, "BEGIN BOND")?;
    cursor += 1;
    for _ in 0..n_bonds {
        let bond = parse_v3000_bond_line(
            v3000_content(lines.get(cursor).ok_or_else(|| {
                SdfReadError::Parse("EOF hit while reading V3000 bonds".to_owned())
            })?)
            .ok_or_else(|| SdfReadError::Parse("Bad V3000 bond line".to_owned()))?,
            cursor + 1,
        )?;
        if matches!(bond.order, BondOrder::Aromatic) {
            molecule.atoms[bond.begin_atom].is_aromatic = true;
            molecule.atoms[bond.end_atom].is_aromatic = true;
        }
        molecule.add_bond(bond);
        cursor += 1;
    }
    expect_v3000_line(lines, cursor, "END BOND")?;
    cursor += 1;

    expect_v3000_line(lines, cursor, "END CTAB")?;
    cursor += 1;
    if lines
        .get(cursor)
        .is_some_and(|line| line.starts_with("M  END"))
    {
        Ok((molecule, cursor + 1))
    } else {
        Err(SdfReadError::Parse(format!(
            "Problems encountered parsing Mol data, M  END missing around line {}",
            cursor + 1
        )))
    }
}

fn parse_v3000_atom_line(line: &str, line_number: usize) -> Result<(Atom, DVec2), SdfReadError> {
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
    parse_f64_token(tokens[4], line_number)?;

    let mut formal_charge = 0;
    let mut radicals = 0;
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
            "CFG" => {
                return Err(SdfReadError::Parse(format!(
                    "V3000 atom CFG is not implemented on line {line_number}"
                )));
            }
            _ => {
                return Err(SdfReadError::Parse(format!(
                    "Unsupported V3000 atom property '{prop}' on line {line_number}"
                )));
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
        },
        DVec2::new(x, y),
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
        "H" => Some((1, None)),
        "B" => Some((5, None)),
        "C" => Some((6, None)),
        "N" => Some((7, None)),
        "O" => Some((8, None)),
        "F" => Some((9, None)),
        "Na" => Some((11, None)),
        "Si" => Some((14, None)),
        "P" => Some((15, None)),
        "S" => Some((16, None)),
        "Cl" => Some((17, None)),
        "Cu" => Some((29, None)),
        "Se" => Some((34, None)),
        "Br" => Some((35, None)),
        "Rh" => Some((45, None)),
        "I" => Some((53, None)),
        _ => None,
    }
}
