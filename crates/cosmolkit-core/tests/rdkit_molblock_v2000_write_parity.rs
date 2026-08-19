use std::fs::File;
use std::io::{BufRead, BufReader};

use cosmolkit_core::Molecule;
use cosmolkit_core::io::molblock::mol_to_v2000_block;
use cosmolkit_test_support::{golden_path, smiles_path};
use serde::Deserialize;

fn body(block: &str) -> String {
    let lines: Vec<_> = block.lines().collect();
    lines[3..].join("\n")
}

fn normalize_signed_zero(text: &str) -> String {
    text.replace("-0.0000", " 0.0000")
}

fn atom_symbol_from_v2000_line(line: &str) -> String {
    line.get(31..34).unwrap_or_default().trim().to_owned()
}

fn atom_symbols_equivalent(ours: &str, expected: &str) -> bool {
    ours == expected || matches!((ours, expected), ("*", "R") | ("R", "*"))
}

fn parse_bond_line(line: &str) -> (usize, usize, usize) {
    let a = line[0..3].trim().parse::<usize>().expect("bond a index");
    let b = line[3..6].trim().parse::<usize>().expect("bond b index");
    let order = line[6..9].trim().parse::<usize>().expect("bond order");
    let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
    (lo, hi, order)
}

#[derive(Debug)]
struct ParsedBody {
    atoms: Vec<String>,
    coords: Vec<(f64, f64)>,
    bonds: Vec<(usize, usize, usize)>,
}

fn parse_v2000_body(block_body: &str) -> Option<ParsedBody> {
    let lines: Vec<_> = block_body.lines().collect();
    let counts = *lines.first()?;
    let atom_count = counts.get(0..3)?.trim().parse::<usize>().ok()?;
    let bond_count = counts.get(3..6)?.trim().parse::<usize>().ok()?;
    if lines.len() < 1 + atom_count + bond_count {
        return None;
    }

    let mut atoms = Vec::with_capacity(atom_count);
    let mut coords = Vec::with_capacity(atom_count);
    for line in &lines[1..1 + atom_count] {
        atoms.push(atom_symbol_from_v2000_line(line));
        coords.push((
            line[0..10].trim().parse::<f64>().ok()?,
            line[10..20].trim().parse::<f64>().ok()?,
        ));
    }
    let bonds = lines[1 + atom_count..1 + atom_count + bond_count]
        .iter()
        .map(|line| parse_bond_line(line))
        .collect();
    Some(ParsedBody {
        atoms,
        coords,
        bonds,
    })
}

fn parse_v3000_body(block_body: &str) -> Option<ParsedBody> {
    let mut in_atom = false;
    let mut in_bond = false;
    let mut atoms_raw: Vec<(usize, String, f64, f64)> = Vec::new();
    let mut bonds = Vec::new();

    for line in block_body.lines() {
        let trimmed = line.trim();
        match trimmed {
            "M  V30 BEGIN ATOM" => {
                in_atom = true;
                continue;
            }
            "M  V30 END ATOM" => {
                in_atom = false;
                continue;
            }
            "M  V30 BEGIN BOND" => {
                in_bond = true;
                continue;
            }
            "M  V30 END BOND" => {
                in_bond = false;
                continue;
            }
            _ => {}
        }
        if in_atom && trimmed.starts_with("M  V30 ") {
            let tokens: Vec<_> = trimmed.split_whitespace().collect();
            if tokens.len() < 6 {
                return None;
            }
            atoms_raw.push((
                tokens[2].parse::<usize>().ok()?,
                tokens[3].to_owned(),
                tokens[4].parse::<f64>().ok()?,
                tokens[5].parse::<f64>().ok()?,
            ));
        } else if in_bond && trimmed.starts_with("M  V30 ") {
            let tokens: Vec<_> = trimmed.split_whitespace().collect();
            if tokens.len() < 6 {
                return None;
            }
            let order = tokens[3].parse::<usize>().ok()?;
            let a = tokens[4].parse::<usize>().ok()?;
            let b = tokens[5].parse::<usize>().ok()?;
            let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
            bonds.push((lo, hi, order));
        }
    }

    if atoms_raw.is_empty() {
        return None;
    }
    atoms_raw.sort_by_key(|(index, _, _, _)| *index);
    let mut atoms = Vec::with_capacity(atoms_raw.len());
    let mut coords = Vec::with_capacity(atoms_raw.len());
    for (_, symbol, x, y) in atoms_raw {
        atoms.push(symbol);
        coords.push((x, y));
    }
    Some(ParsedBody {
        atoms,
        coords,
        bonds,
    })
}

fn parse_body_for_compare(block_body: &str) -> Option<ParsedBody> {
    if block_body
        .lines()
        .next()
        .unwrap_or_default()
        .contains("V3000")
    {
        parse_v3000_body(block_body)
    } else {
        parse_v2000_body(block_body)
    }
}

fn canonical_bonds(bonds: &[(usize, usize, usize)]) -> Vec<(usize, usize, usize)> {
    let mut bonds = bonds.to_vec();
    bonds.sort_unstable();
    bonds
}

fn bodies_match(ours_body: &str, expected_body: &str) -> bool {
    let Some(ours) = parse_body_for_compare(ours_body) else {
        return false;
    };
    let Some(expected) = parse_body_for_compare(expected_body) else {
        return false;
    };
    ours.atoms.len() == expected.atoms.len()
        && ours.coords.len() == expected.coords.len()
        && ours.bonds.len() == expected.bonds.len()
        && ours
            .atoms
            .iter()
            .zip(&expected.atoms)
            .all(|(ours, expected)| atom_symbols_equivalent(ours, expected))
        && ours
            .coords
            .iter()
            .zip(&expected.coords)
            .all(|((ox, oy), (ex, ey))| (ox - ex).abs() <= 1e-3 && (oy - ey).abs() <= 1e-3)
        && canonical_bonds(&ours.bonds) == canonical_bonds(&expected.bonds)
}

fn compare_against_expected(
    ours_body: &str,
    expected_body: &str,
    smiles: &str,
    row: usize,
    variant: &str,
) {
    let ours = normalize_signed_zero(ours_body);
    let expected = normalize_signed_zero(expected_body);
    assert!(
        bodies_match(&ours, &expected),
        "molblock mismatch at row {row} ({smiles}) against {variant}\nours:\n{ours}\nexpected:\n{expected}"
    );
}

#[derive(Debug, Deserialize)]
struct GoldenRecord {
    smiles: String,
    parse_ok: bool,
    parse_error: Option<String>,
    v2000_ok: bool,
    v2000_body: Option<String>,
    v2000_error: Option<String>,
    v3000_ok: bool,
    v3000_body: Option<String>,
    v3000_error: Option<String>,
}

fn load_smiles() -> Vec<String> {
    let path = smiles_path();
    std::fs::read_to_string(&path)
        .unwrap_or_else(|error| panic!("failed to read {}: {error}", path.display()))
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(ToOwned::to_owned)
        .collect()
}

fn load_golden() -> Vec<GoldenRecord> {
    let path = golden_path("molblock_v2000_minimal.jsonl");
    BufReader::new(
        File::open(&path)
            .unwrap_or_else(|error| panic!("failed to open {}: {error}", path.display())),
    )
    .lines()
    .enumerate()
    .map(|(index, line)| {
        let line = line.unwrap_or_else(|error| {
            panic!(
                "failed to read {} line {}: {error}",
                path.display(),
                index + 1
            )
        });
        serde_json::from_str(&line).unwrap_or_else(|error| {
            panic!(
                "failed to parse {} line {}: {error}",
                path.display(),
                index + 1
            )
        })
    })
    .collect()
}

#[test]
fn molblock_v2000_golden_has_one_record_per_smiles() {
    let smiles = load_smiles();
    let golden = load_golden();
    assert_eq!(golden.len(), smiles.len());
    for (row, (record, input)) in golden.iter().zip(&smiles).enumerate() {
        assert_eq!(
            &record.smiles,
            input,
            "golden smiles mismatch at row {}",
            row + 1
        );
        assert_eq!(
            record.parse_ok,
            record.parse_error.is_none(),
            "parse fields at row {}",
            row + 1
        );
        assert_eq!(
            record.v2000_ok,
            record.v2000_body.is_some(),
            "V2000 body at row {}",
            row + 1
        );
        assert_eq!(
            record.v2000_ok,
            record.v2000_error.is_none(),
            "V2000 error at row {}",
            row + 1
        );
        assert_eq!(
            record.v3000_ok,
            record.v3000_body.is_some(),
            "V3000 body at row {}",
            row + 1
        );
        assert_eq!(
            record.v3000_ok,
            record.v3000_error.is_none(),
            "V3000 error at row {}",
            row + 1
        );
    }
}

#[test]
fn molblock_v2000_body_matches_rdkit_coordinates_and_topology() {
    for (index, record) in load_golden().iter().enumerate() {
        let row = index + 1;
        let parsed = Molecule::from_smiles(&record.smiles);
        if !record.parse_ok {
            assert!(
                parsed.is_err(),
                "parse should fail at row {row} ({})",
                record.smiles
            );
            continue;
        }
        let molecule = parsed
            .unwrap_or_else(|error| panic!("parse failed at row {row}: {error}"))
            .with_2d_coordinates()
            .unwrap_or_else(|error| panic!("2D coordinates failed at row {row}: {error}"));
        let block = mol_to_v2000_block(&molecule);
        if !record.v2000_ok && !record.v3000_ok {
            assert!(
                block.is_err(),
                "write should fail at row {row} ({}); RDKit V2000 error: {:?}; RDKit V3000 error: {:?}",
                record.smiles,
                record.v2000_error,
                record.v3000_error
            );
            continue;
        }
        let block = block.unwrap_or_else(|error| panic!("write failed at row {row}: {error}"));
        let actual_body = body(&block);
        if let Some(expected) = &record.v2000_body {
            compare_against_expected(&actual_body, expected, &record.smiles, row, "v2000");
        }
        if let Some(expected) = &record.v3000_body {
            compare_against_expected(&actual_body, expected, &record.smiles, row, "v3000");
        }
    }
}
