use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::process::Command;

use super::mol_to_v2000_block;
use crate::Molecule;
use serde::Deserialize;

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

#[derive(Debug, Deserialize)]
struct KekulizeGoldenRecord {
    smiles: String,
    parse_ok: bool,
    parse_error: Option<String>,
    kekulize_ok: bool,
    kekulize_error: Option<String>,
    v2000_ok: bool,
    v2000_body: Option<String>,
    v2000_error: Option<String>,
    v3000_ok: bool,
    v3000_body: Option<String>,
    v3000_error: Option<String>,
}

#[derive(Debug, thiserror::Error)]
enum TestDataError {
    #[error(transparent)]
    Io(#[from] std::io::Error),
    #[error("invalid jsonl at line {line_no}: {source}")]
    Json {
        line_no: usize,
        #[source]
        source: serde_json::Error,
    },
}

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("crates/")
        .parent()
        .expect("repo root")
        .to_path_buf()
}

fn body(block: &str) -> String {
    let lines: Vec<_> = block.lines().collect();
    lines[3..].join("\n")
}

fn normalize_signed_zero(text: &str) -> String {
    text.replace("-0.0000", " 0.0000")
}

fn atom_symbol_from_v2000_line(line: &str) -> String {
    if line.len() >= 34 {
        return line[31..34].trim().to_owned();
    }
    String::new()
}

fn atom_symbols_equivalent(ours: &str, expected: &str) -> bool {
    if ours == expected {
        return true;
    }
    // RDKit writes mapped dummies as "R" atoms in V2000 while our minimal
    // subset keeps "*" placeholders.
    (ours == "*" && expected == "R") || (ours == "R" && expected == "*")
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
    for i in 0..atom_count {
        let line = lines[1 + i];
        atoms.push(atom_symbol_from_v2000_line(line));
        let x = line[0..10].trim().parse::<f64>().ok()?;
        let y = line[10..20].trim().parse::<f64>().ok()?;
        coords.push((x, y));
    }

    let mut bonds = Vec::with_capacity(bond_count);
    for i in 0..bond_count {
        bonds.push(parse_bond_line(lines[1 + atom_count + i]));
    }
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
    let mut bonds: Vec<(usize, usize, usize)> = Vec::new();

    for line in block_body.lines() {
        let t = line.trim();
        if t == "M  V30 BEGIN ATOM" {
            in_atom = true;
            continue;
        }
        if t == "M  V30 END ATOM" {
            in_atom = false;
            continue;
        }
        if t == "M  V30 BEGIN BOND" {
            in_bond = true;
            continue;
        }
        if t == "M  V30 END BOND" {
            in_bond = false;
            continue;
        }

        if in_atom && t.starts_with("M  V30 ") {
            let toks: Vec<_> = t.split_whitespace().collect();
            if toks.len() < 5 {
                return None;
            }
            let idx = toks[2].parse::<usize>().ok()?;
            let symbol = toks[3].to_owned();
            let x = toks[4].parse::<f64>().ok()?;
            let y = toks[5].parse::<f64>().ok()?;
            atoms_raw.push((idx, symbol, x, y));
        } else if in_bond && t.starts_with("M  V30 ") {
            let toks: Vec<_> = t.split_whitespace().collect();
            if toks.len() < 6 {
                return None;
            }
            let order = toks[3].parse::<usize>().ok()?;
            let a = toks[4].parse::<usize>().ok()?;
            let b = toks[5].parse::<usize>().ok()?;
            let (lo, hi) = if a <= b { (a, b) } else { (b, a) };
            bonds.push((lo, hi, order));
        }
    }

    if atoms_raw.is_empty() {
        return None;
    }
    atoms_raw.sort_by_key(|(idx, _, _, _)| *idx);
    let mut atoms = Vec::with_capacity(atoms_raw.len());
    let mut coords = Vec::with_capacity(atoms_raw.len());
    for (_, s, x, y) in atoms_raw {
        atoms.push(s);
        coords.push((x, y));
    }
    Some(ParsedBody {
        atoms,
        coords,
        bonds,
    })
}

fn parse_body_for_compare(block_body: &str) -> Option<ParsedBody> {
    let first = block_body.lines().next().unwrap_or_default();
    if first.contains("V3000") {
        parse_v3000_body(block_body)
    } else {
        parse_v2000_body(block_body)
    }
}

fn canonical_bonds_for_compare(bonds: &[(usize, usize, usize)]) -> Vec<(usize, usize, usize)> {
    let mut out = bonds.to_vec();
    out.sort_unstable();
    out
}

fn coords_match_strict(ours: &[(f64, f64)], expected: &[(f64, f64)]) -> bool {
    if ours.len() != expected.len() {
        return false;
    }
    let coord_tol = 1e-3f64;
    for i in 0..ours.len() {
        let (ox, oy) = ours[i];
        let (ex, ey) = expected[i];
        if (ox - ex).abs() > coord_tol || (oy - ey).abs() > coord_tol {
            return false;
        }
    }
    true
}

fn align_rigid_2d(source: &[(f64, f64)], target: &[(f64, f64)]) -> Option<(f64, f64, f64, f64)> {
    if source.len() != target.len() || source.is_empty() {
        return None;
    }

    let n = source.len() as f64;
    let src_cx = source.iter().map(|(x, _)| x).sum::<f64>() / n;
    let src_cy = source.iter().map(|(_, y)| y).sum::<f64>() / n;
    let tgt_cx = target.iter().map(|(x, _)| x).sum::<f64>() / n;
    let tgt_cy = target.iter().map(|(_, y)| y).sum::<f64>() / n;

    let mut a = 0.0f64;
    let mut b = 0.0f64;
    for i in 0..source.len() {
        let sx = source[i].0 - src_cx;
        let sy = source[i].1 - src_cy;
        let tx = target[i].0 - tgt_cx;
        let ty = target[i].1 - tgt_cy;
        a += sx * tx + sy * ty;
        b += sx * ty - sy * tx;
    }

    let norm = (a * a + b * b).sqrt();
    if norm <= 1e-12 {
        return None;
    }
    let cos_t = a / norm;
    let sin_t = b / norm;
    let tx = tgt_cx - (cos_t * src_cx - sin_t * src_cy);
    let ty = tgt_cy - (sin_t * src_cx + cos_t * src_cy);
    Some((cos_t, sin_t, tx, ty))
}

fn apply_rigid_2d(p: (f64, f64), tf: (f64, f64, f64, f64)) -> (f64, f64) {
    let (cos_t, sin_t, tx, ty) = tf;
    let x = cos_t * p.0 - sin_t * p.1 + tx;
    let y = sin_t * p.0 + cos_t * p.1 + ty;
    (x, y)
}

fn rmsd_for_indices(
    ours: &[(f64, f64)],
    expected: &[(f64, f64)],
    indices_1based: &[usize],
    tf: (f64, f64, f64, f64),
) -> f64 {
    let mut sse = 0.0f64;
    for &idx in indices_1based {
        let i = idx - 1;
        let po = apply_rigid_2d(ours[i], tf);
        let pe = expected[i];
        let dx = po.0 - pe.0;
        let dy = po.1 - pe.1;
        sse += dx * dx + dy * dy;
    }
    (sse / (indices_1based.len() as f64)).sqrt()
}

fn non_coordinate_sections_match(ours_body: &str, expected_body: &str) -> bool {
    let Some(ours) = parse_body_for_compare(ours_body) else {
        return false;
    };
    let Some(expected) = parse_body_for_compare(expected_body) else {
        return false;
    };

    if ours.atoms.len() != expected.atoms.len()
        || ours.coords.len() != expected.coords.len()
        || ours.bonds.len() != expected.bonds.len()
    {
        return false;
    }
    let coords_ok = coords_match_strict(&ours.coords, &expected.coords);

    for i in 0..ours.atoms.len() {
        let ours_symbol = &ours.atoms[i];
        let expected_symbol = &expected.atoms[i];
        if !atom_symbols_equivalent(ours_symbol, expected_symbol) {
            return false;
        }
    }

    canonical_bonds_for_compare(&ours.bonds) == canonical_bonds_for_compare(&expected.bonds)
        && coords_ok
}

fn compare_against_expected(
    ours_body: &str,
    expected_body: &str,
    smiles: &str,
    row_idx_1based: usize,
    variant: &str,
) {
    let ours_norm = normalize_signed_zero(ours_body);
    let expected_norm = normalize_signed_zero(expected_body);
    let mismatch_detail = {
        let ours = parse_body_for_compare(&ours_norm);
        let expected = parse_body_for_compare(&expected_norm);
        match (ours, expected) {
            (Some(ours), Some(expected)) => {
                let mut detail = String::new();
                if ours.atoms.len() == expected.atoms.len() {
                    if let Some((idx, (o, e))) = ours
                        .coords
                        .iter()
                        .zip(expected.coords.iter())
                        .enumerate()
                        .find(|(_, (o, e))| {
                            (o.0 - e.0).abs() > 1e-3f64 || (o.1 - e.1).abs() > 1e-3f64
                        })
                    {
                        detail.push_str(&format!(
                                "first coordinate mismatch at atom {}: ours=({:.4},{:.4}) expected=({:.4},{:.4})\n",
                                idx + 1,
                                o.0,
                                o.1,
                                e.0,
                                e.1
                            ));
                    }
                }
                detail.push_str("ours body:\n");
                detail.push_str(&ours_norm);
                detail.push_str("\nexpected body:\n");
                detail.push_str(&expected_norm);
                detail
            }
            _ => format!(
                "failed to parse ours/expected body\nours:\n{ours_norm}\nexpected:\n{expected_norm}"
            ),
        }
    };
    assert!(
        non_coordinate_sections_match(&ours_norm, &expected_norm),
        "molblock mismatch (including coordinates) at row {} ({}) against {}\n{}",
        row_idx_1based,
        smiles,
        variant,
        mismatch_detail
    );
}

fn bonds_and_atoms_match_strict(ours_body: &str, expected_body: &str) -> bool {
    let Some(ours) = parse_body_for_compare(ours_body) else {
        return false;
    };
    let Some(expected) = parse_body_for_compare(expected_body) else {
        return false;
    };

    if ours.atoms.len() != expected.atoms.len() || ours.bonds.len() != expected.bonds.len() {
        return false;
    }
    for i in 0..ours.atoms.len() {
        if !atom_symbols_equivalent(&ours.atoms[i], &expected.atoms[i]) {
            return false;
        }
    }
    let mut ours_bonds = ours.bonds;
    let mut expected_bonds = expected.bonds;
    ours_bonds.sort_unstable();
    expected_bonds.sort_unstable();
    ours_bonds == expected_bonds
}

fn load_smiles() -> Result<Vec<String>, TestDataError> {
    let path = repo_root().join("tests/smiles.smi");
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut rows = Vec::new();
    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        rows.push(trimmed.to_owned());
    }
    Ok(rows)
}

fn load_golden() -> Result<Vec<GoldenRecord>, TestDataError> {
    let path = repo_root().join("tests/golden/molblock_v2000_minimal.jsonl");
    ensure_golden_exists(&path);
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut rows = Vec::new();
    for (idx, line) in reader.lines().enumerate() {
        let line_no = idx + 1;
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let record = serde_json::from_str::<GoldenRecord>(&line)
            .map_err(|source| TestDataError::Json { line_no, source })?;
        rows.push(record);
    }
    Ok(rows)
}

fn load_kekulize_golden() -> Result<Vec<KekulizeGoldenRecord>, TestDataError> {
    let path = repo_root().join("tests/golden/molblock_v2000_kekulized.jsonl");
    ensure_kekulize_golden_exists(&path);
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut rows = Vec::new();
    for (idx, line) in reader.lines().enumerate() {
        let line_no = idx + 1;
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let record = serde_json::from_str::<KekulizeGoldenRecord>(&line)
            .map_err(|source| TestDataError::Json { line_no, source })?;
        rows.push(record);
    }
    Ok(rows)
}

fn ensure_golden_exists(golden_path: &Path) {
    if golden_path.exists() {
        return;
    }
    let repo = repo_root();
    let script = repo.join("tests/scripts/gen_rdkit_v2000_minimal_golden.py");
    let input = repo.join("tests/smiles.smi");

    let candidates = [
        std::env::var("COSMOLKIT_PYTHON").ok(),
        Some(repo.join(".venv/bin/python").display().to_string()),
        Some(String::from("python3")),
    ];

    let mut last_error = String::new();
    for candidate in candidates.iter().flatten() {
        let output = Command::new(candidate)
            .arg(&script)
            .arg("--input")
            .arg(&input)
            .arg("--output")
            .arg(golden_path)
            .output();
        match output {
            Ok(out) if out.status.success() => return,
            Ok(out) => {
                last_error = format!(
                    "python={} exit={} stderr={}",
                    candidate,
                    out.status,
                    String::from_utf8_lossy(&out.stderr)
                );
            }
            Err(err) => {
                last_error = format!("python={} spawn error={}", candidate, err);
            }
        }
    }

    panic!(
        "golden file missing and auto-generation failed.\n\
             expected: {}\n\
             tried COSMOLKIT_PYTHON, .venv/bin/python, python3.\n\
             last error: {}\n\
             please run:\n\
             uv sync --group dev && .venv/bin/python tests/scripts/gen_rdkit_v2000_minimal_golden.py --input tests/smiles.smi --output tests/golden/molblock_v2000_minimal.jsonl",
        golden_path.display(),
        last_error
    );
}

fn ensure_kekulize_golden_exists(golden_path: &Path) {
    if golden_path.exists() {
        return;
    }
    let repo = repo_root();
    let script = repo.join("tests/scripts/gen_rdkit_kekulize_molblock_golden.py");
    let input = repo.join("tests/smiles.smi");

    let candidates = [
        std::env::var("COSMOLKIT_PYTHON").ok(),
        Some(repo.join(".venv/bin/python").display().to_string()),
        Some(String::from("python3")),
    ];

    let mut last_error = String::new();
    for candidate in candidates.iter().flatten() {
        let output = Command::new(candidate)
            .arg(&script)
            .arg("--input")
            .arg(&input)
            .arg("--output")
            .arg(golden_path)
            .output();
        match output {
            Ok(out) if out.status.success() => return,
            Ok(out) => {
                last_error = format!(
                    "python={} exit={} stderr={}",
                    candidate,
                    out.status,
                    String::from_utf8_lossy(&out.stderr)
                );
            }
            Err(err) => {
                last_error = format!("python={} spawn error={}", candidate, err);
            }
        }
    }

    panic!(
        "kekulize golden file missing and auto-generation failed.\n\
             expected: {}\n\
             tried COSMOLKIT_PYTHON, .venv/bin/python, python3.\n\
             last error: {}\n\
             please run:\n\
             uv sync --group dev && .venv/bin/python tests/scripts/gen_rdkit_kekulize_molblock_golden.py --input tests/smiles.smi --output tests/golden/molblock_v2000_kekulized.jsonl",
        golden_path.display(),
        last_error
    );
}

#[test]
fn molblock_v2000_golden_has_one_record_per_smiles() {
    let smiles = load_smiles().expect("read tests/smiles.smi");
    let golden = load_golden().expect("read tests/golden/molblock_v2000_minimal.jsonl");
    assert_eq!(
        golden.len(),
        smiles.len(),
        "golden rows must match input smiles rows"
    );

    for (idx, (record, input_smiles)) in golden.iter().zip(smiles.iter()).enumerate() {
        assert_eq!(
            record.smiles,
            *input_smiles,
            "smiles mismatch at row {}",
            idx + 1
        );

        if record.parse_ok {
            assert!(
                record.parse_error.is_none(),
                "parse_ok=true should not carry parse_error at row {}",
                idx + 1
            );
        } else {
            assert!(
                record.parse_error.is_some(),
                "parse_ok=false should carry parse_error at row {}",
                idx + 1
            );
        }

        if record.v2000_ok {
            assert!(
                record.v2000_body.is_some(),
                "v2000_ok=true requires v2000_body at row {}",
                idx + 1
            );
            assert!(
                record.v2000_error.is_none(),
                "v2000_ok=true should not carry v2000_error at row {}",
                idx + 1
            );
        } else {
            assert!(
                record.v2000_body.is_none(),
                "v2000_ok=false should not carry v2000_body at row {}",
                idx + 1
            );
            assert!(
                record.v2000_error.is_some(),
                "v2000_ok=false should carry v2000_error at row {}",
                idx + 1
            );
        }

        if record.v3000_ok {
            assert!(
                record.v3000_body.is_some(),
                "v3000_ok=true requires v3000_body at row {}",
                idx + 1
            );
            assert!(
                record.v3000_error.is_none(),
                "v3000_ok=true should not carry v3000_error at row {}",
                idx + 1
            );
        } else {
            assert!(
                record.v3000_body.is_none(),
                "v3000_ok=false should not carry v3000_body at row {}",
                idx + 1
            );
            assert!(
                record.v3000_error.is_some(),
                "v3000_ok=false should carry v3000_error at row {}",
                idx + 1
            );
        }
    }
}

#[test]
fn molblock_v2000_body_matches_rdkit_coordinates_and_topology() {
    let golden = load_golden().expect("read tests/golden/molblock_v2000_minimal.jsonl");
    for (idx, record) in golden.iter().enumerate() {
        let parsed = Molecule::from_smiles(&record.smiles);
        if record.parse_ok {
            assert!(
                parsed.is_ok(),
                "parse should succeed at row {} ({})",
                idx + 1,
                record.smiles
            );
        } else {
            assert!(
                parsed.is_err(),
                "parse should fail at row {} ({})",
                idx + 1,
                record.smiles
            );
            continue;
        }
        let mut mol = parsed.expect("parse checked above");
        mol.compute_2d_coords().unwrap_or_else(|e| {
            panic!("2D coordinate generation failed at row {}: {}", idx + 1, e)
        });
        let ours = mol_to_v2000_block(&mol)
            .unwrap_or_else(|e| panic!("write failed at row {}: {}", idx + 1, e));
        let ours_body = body(&ours);

        if record.v2000_ok {
            let expected_v2000 = record
                .v2000_body
                .as_ref()
                .expect("v2000_ok=true requires v2000_body");
            compare_against_expected(&ours_body, expected_v2000, &record.smiles, idx + 1, "v2000");
        }
        if record.v3000_ok {
            let expected_v3000 = record
                .v3000_body
                .as_ref()
                .expect("v3000_ok=true requires v3000_body");
            compare_against_expected(&ours_body, expected_v3000, &record.smiles, idx + 1, "v3000");
        }
    }
}

#[test]
fn molblock_kekulized_golden_has_one_record_per_smiles() {
    let smiles = load_smiles().expect("read tests/smiles.smi");
    let golden = load_kekulize_golden().expect("read tests/golden/molblock_v2000_kekulized.jsonl");
    assert_eq!(
        golden.len(),
        smiles.len(),
        "kekulize golden rows must match input smiles rows"
    );
    for (idx, (record, input_smiles)) in golden.iter().zip(smiles.iter()).enumerate() {
        assert_eq!(
            record.smiles,
            *input_smiles,
            "kekulize smiles mismatch at row {}",
            idx + 1
        );
        if record.parse_ok {
            assert!(
                record.parse_error.is_none(),
                "kekulize parse_ok=true should not carry parse_error at row {}",
                idx + 1
            );
        } else {
            assert!(
                record.parse_error.is_some(),
                "kekulize parse_ok=false should carry parse_error at row {}",
                idx + 1
            );
        }
        if record.kekulize_ok {
            assert!(
                record.kekulize_error.is_none(),
                "kekulize_ok=true should not carry kekulize_error at row {}",
                idx + 1
            );
            if record.v2000_ok {
                assert!(
                    record.v2000_body.is_some() && record.v2000_error.is_none(),
                    "kekulize v2000_ok=true requires body and no error at row {}",
                    idx + 1
                );
            } else {
                assert!(
                    record.v2000_body.is_none() && record.v2000_error.is_some(),
                    "kekulize v2000_ok=false requires error and no body at row {}",
                    idx + 1
                );
            }
            if record.v3000_ok {
                assert!(
                    record.v3000_body.is_some() && record.v3000_error.is_none(),
                    "kekulize v3000_ok=true requires body and no error at row {}",
                    idx + 1
                );
            } else {
                assert!(
                    record.v3000_body.is_none() && record.v3000_error.is_some(),
                    "kekulize v3000_ok=false requires error and no body at row {}",
                    idx + 1
                );
            }
        } else {
            assert!(
                record.kekulize_error.is_some(),
                "kekulize_ok=false should carry kekulize_error at row {}",
                idx + 1
            );
            assert!(
                !record.v2000_ok && !record.v3000_ok,
                "kekulize_ok=false should not have successful molblock variants at row {}",
                idx + 1
            );
        }
    }
}

#[test]
fn molblock_kekulized_topology_matches_rdkit_golden() {
    let golden = load_kekulize_golden().expect("read tests/golden/molblock_v2000_kekulized.jsonl");
    for (idx, record) in golden.iter().enumerate() {
        let parsed = Molecule::from_smiles(&record.smiles);
        if record.parse_ok {
            assert!(
                parsed.is_ok(),
                "parse should succeed at row {} ({})",
                idx + 1,
                record.smiles
            );
        } else {
            assert!(
                parsed.is_err(),
                "parse should fail at row {} ({})",
                idx + 1,
                record.smiles
            );
            continue;
        }
        let mut mol = parsed.expect("parse checked above");

        let ours_kek = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            crate::kekulize::kekulize_in_place(&mut mol, true)
        }));
        let kek_ok = matches!(ours_kek, Ok(Ok(())));

        if record.kekulize_ok {
            assert!(
                kek_ok,
                "kekulize should succeed at row {} ({})",
                idx + 1,
                record.smiles
            );
        } else {
            assert!(
                !kek_ok,
                "kekulize should fail at row {} ({})",
                idx + 1,
                record.smiles
            );
            continue;
        }

        mol.compute_2d_coords().unwrap_or_else(|e| {
            panic!("2D coordinate generation failed at row {}: {}", idx + 1, e)
        });
        let ours = mol_to_v2000_block(&mol)
            .unwrap_or_else(|e| panic!("write failed at row {}: {}", idx + 1, e));
        let ours_body = body(&ours);

        if record.v2000_ok {
            let expected_v2000 = record
                .v2000_body
                .as_ref()
                .expect("v2000_ok=true requires v2000_body");
            let ours_norm = normalize_signed_zero(&ours_body);
            let expected_norm = normalize_signed_zero(expected_v2000);
            assert!(
                bonds_and_atoms_match_strict(&ours_norm, &expected_norm),
                "kekulized bond block mismatch at row {} ({}) against v2000\nours:\n{}\nexpected:\n{}",
                idx + 1,
                record.smiles,
                ours_norm,
                expected_norm
            );
        }
        if record.v3000_ok {
            let expected_v3000 = record
                .v3000_body
                .as_ref()
                .expect("v3000_ok=true requires v3000_body");
            let ours_norm = normalize_signed_zero(&ours_body);
            let expected_norm = normalize_signed_zero(expected_v3000);
            assert!(
                bonds_and_atoms_match_strict(&ours_norm, &expected_norm),
                "kekulized bond block mismatch at row {} ({}) against v3000\nours:\n{}\nexpected:\n{}",
                idx + 1,
                record.smiles,
                ours_norm,
                expected_norm
            );
        }
    }
}

#[test]
fn long_biaryl_ester_anchor_aligned_ring_rmsd_diagnostic() {
    let golden = load_golden().expect("read tests/golden/molblock_v2000_minimal.jsonl");
    let diagnostic_smiles = "CCCCCCCCOc1ccc(C(=O)COC(=O)c2ccc(CCCCCCC)cc2)cc1";
    let record = golden
        .iter()
        .find(|record| record.smiles == diagnostic_smiles)
        .expect("diagnostic molecule should exist in golden");
    assert!(
        record.v2000_ok,
        "diagnostic molecule should have v2000 golden body"
    );

    let mut mol = Molecule::from_smiles(&record.smiles).expect("parse diagnostic molecule");
    mol.compute_2d_coords()
        .expect("compute diagnostic 2D coordinates");
    let ours_block = mol_to_v2000_block(&mol).expect("write diagnostic molecule");
    let ours = parse_v2000_body(&body(&ours_block)).expect("parse ours diagnostic body");
    let expected = parse_v2000_body(
        record
            .v2000_body
            .as_ref()
            .expect("diagnostic v2000_ok=true requires body"),
    )
    .expect("parse rdkit diagnostic body");

    // O=C-CO-C=O six-atom anchor: [O, C, C, O, C, O]
    let anchor = [17usize, 16, 18, 19, 20, 21];
    let ring_a = [10usize, 11, 12, 13, 33, 34];
    let ring_b = [20usize, 21, 22, 23, 31, 32];

    let src_anchor: Vec<(f64, f64)> = anchor.iter().map(|&i| ours.coords[i - 1]).collect();
    let tgt_anchor: Vec<(f64, f64)> = anchor.iter().map(|&i| expected.coords[i - 1]).collect();
    let tf = align_rigid_2d(&src_anchor, &tgt_anchor).expect("anchor rigid align");
    let anchor_rmsd = rmsd_for_indices(&ours.coords, &expected.coords, &anchor, tf);
    let ring_a_rmsd = rmsd_for_indices(&ours.coords, &expected.coords, &ring_a, tf);
    let ring_b_rmsd = rmsd_for_indices(&ours.coords, &expected.coords, &ring_b, tf);

    assert!(
        anchor_rmsd <= 1e-4 && ring_a_rmsd <= 1e-4 && ring_b_rmsd <= 1e-4,
        "long biaryl ester anchor-fit RMSD: anchor={:.6}, ring1={:.6}, ring2={:.6}",
        anchor_rmsd,
        ring_a_rmsd,
        ring_b_rmsd
    );
}

#[test]
fn long_biaryl_ester_coordinate_diagnostic() {
    let golden = load_golden().expect("read tests/golden/molblock_v2000_minimal.jsonl");
    let diagnostic_smiles = "CCCCCCCCOc1ccc(C(=O)COC(=O)c2ccc(CCCCCCC)cc2)cc1";
    let record = golden
        .iter()
        .find(|record| record.smiles == diagnostic_smiles)
        .expect("diagnostic molecule should exist in golden");
    let mut mol = Molecule::from_smiles(&record.smiles).expect("parse diagnostic molecule");
    mol.compute_2d_coords()
        .expect("compute diagnostic 2D coordinates");
    let ours_block = mol_to_v2000_block(&mol).expect("write diagnostic molecule");
    let ours = parse_v2000_body(&body(&ours_block)).expect("parse ours diagnostic body");
    let expected = parse_v2000_body(
        record
            .v2000_body
            .as_ref()
            .expect("diagnostic v2000_ok=true requires body"),
    )
    .expect("parse rdkit diagnostic body");

    assert_eq!(ours.atoms, expected.atoms, "diagnostic atom symbols differ");
    let mut ours_bonds = ours.bonds.clone();
    let mut expected_bonds = expected.bonds.clone();
    ours_bonds.sort_unstable();
    expected_bonds.sort_unstable();
    assert_eq!(ours_bonds, expected_bonds, "diagnostic bonds differ");

    let mut sse = 0.0f64;
    let mut sse_rot = 0.0f64;
    let mut max = (0usize, 0.0f64, (0.0f64, 0.0f64), (0.0f64, 0.0f64));
    let mut max_rot = (0usize, 0.0f64, (0.0f64, 0.0f64), (0.0f64, 0.0f64));
    for i in 0..ours.coords.len() {
        let dx = ours.coords[i].0 - expected.coords[i].0;
        let dy = ours.coords[i].1 - expected.coords[i].1;
        let d = (dx * dx + dy * dy).sqrt();
        sse += d * d;
        if d > max.1 {
            max = (i + 1, d, ours.coords[i], expected.coords[i]);
        }
        let rot = (-ours.coords[i].0, -ours.coords[i].1);
        let dxr = rot.0 - expected.coords[i].0;
        let dyr = rot.1 - expected.coords[i].1;
        let dr = (dxr * dxr + dyr * dyr).sqrt();
        sse_rot += dr * dr;
        if dr > max_rot.1 {
            max_rot = (i + 1, dr, rot, expected.coords[i]);
        }
    }
    let rmsd = (sse / ours.coords.len() as f64).sqrt();
    assert!(
        rmsd <= 1e-4,
        "long biaryl ester direct RMSD={:.6}; rot180 RMSD={:.6}; direct max atom {} d={:.6} ours=({:.4},{:.4}) rdkit=({:.4},{:.4}); rot max atom {} d={:.6} ours=({:.4},{:.4}) rdkit=({:.4},{:.4})",
        rmsd,
        (sse_rot / ours.coords.len() as f64).sqrt(),
        max.0,
        max.1,
        max.2.0,
        max.2.1,
        max.3.0,
        max.3.1,
        max_rot.0,
        max_rot.1,
        max_rot.2.0,
        max_rot.2.1,
        max_rot.3.0,
        max_rot.3.1
    );
}
