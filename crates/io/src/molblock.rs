use core::fmt;

use cosmolkit_chem_core::{BondOrder, Molecule};

#[derive(Debug)]
pub enum MolWriteError {
    UnsupportedSubset(&'static str),
}

impl fmt::Display for MolWriteError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnsupportedSubset(msg) => f.write_str(msg),
        }
    }
}

impl std::error::Error for MolWriteError {}

pub fn mol_to_v2000_block_minimal(mol: &Molecule) -> Result<String, MolWriteError> {
    let coords = compute_2d_coords_minimal(mol)?;

    let mut out = String::new();
    out.push('\n');
    out.push_str("     COSMolKit      2D\n");
    out.push('\n');
    out.push_str(&format!(
        "{:>3}{:>3}  0  0  0  0  0  0  0  0999 V2000\n",
        mol.atoms.len(),
        mol.bonds.len()
    ));

    for (idx, atom) in mol.atoms.iter().enumerate() {
        let (x, y) = coords[idx];
        out.push_str(&format!(
            "{:>10.4}{:>10.4}{:>10.4} {:<3} 0  0  0  0  0  0  0  0  0  0  0  0\n",
            x,
            y,
            0.0,
            atom_symbol(atom.atomic_num)
        ));
    }

    for bond in &mol.bonds {
        out.push_str(&format!(
            "{:>3}{:>3}{:>3}  0\n",
            bond.begin_atom + 1,
            bond.end_atom + 1,
            bond_type_code(bond.order)
        ));
    }
    out.push_str("M  END\n");
    Ok(out)
}

pub fn mol_to_sdf_record_minimal(mol: &Molecule) -> Result<String, MolWriteError> {
    let mut block = mol_to_v2000_block_minimal(mol)?;
    block.push_str("$$$$\n");
    Ok(block)
}

fn compute_2d_coords_minimal(mol: &Molecule) -> Result<Vec<(f64, f64)>, MolWriteError> {
    let n = mol.atoms.len();
    if n == 0 {
        return Err(MolWriteError::UnsupportedSubset(
            "empty molecule is unsupported",
        ));
    }
    if n > 3 {
        return Err(MolWriteError::UnsupportedSubset(
            "minimal deterministic coordinate subset currently supports up to 3 atoms",
        ));
    }
    if n == 1 {
        return Ok(vec![(0.0, 0.0)]);
    }

    // This subset intentionally handles only connected linear paths.
    if mol.bonds.len() != n - 1 {
        return Err(MolWriteError::UnsupportedSubset(
            "minimal deterministic coordinate subset requires a connected linear path",
        ));
    }

    let mut degree = vec![0usize; n];
    for b in &mol.bonds {
        degree[b.begin_atom] += 1;
        degree[b.end_atom] += 1;
    }
    if degree.iter().any(|d| *d > 2) {
        return Err(MolWriteError::UnsupportedSubset(
            "branches/rings are outside the minimal deterministic subset",
        ));
    }

    if n == 2 {
        return Ok(vec![(-0.7500, 0.0), (0.7500, -0.0)]);
    }
    Ok(vec![(-1.2990, -0.2500), (0.0, 0.5000), (1.2990, -0.2500)])
}

fn atom_symbol(atomic_num: u8) -> &'static str {
    match atomic_num {
        1 => "H",
        5 => "B",
        6 => "C",
        7 => "N",
        8 => "O",
        9 => "F",
        15 => "P",
        16 => "S",
        17 => "Cl",
        35 => "Br",
        53 => "I",
        _ => "*",
    }
}

fn bond_type_code(order: BondOrder) -> usize {
    match order {
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Aromatic => 4,
        BondOrder::Dative => 9,
        BondOrder::Null => 1,
        BondOrder::Quadruple => 0,
    }
}

#[cfg(test)]
mod tests {
    use core::fmt;
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    use std::path::{Path, PathBuf};
    use std::process::Command;

    use super::mol_to_v2000_block_minimal;
    use cosmolkit_chem_core::Molecule;
    use serde::Deserialize;

    #[derive(Debug, Deserialize)]
    struct GoldenRecord {
        smiles: String,
        rdkit_ok: bool,
        body: Option<String>,
        error: Option<String>,
    }

    #[derive(Debug)]
    enum TestDataError {
        Io(std::io::Error),
        Json {
            line_no: usize,
            source: serde_json::Error,
        },
    }

    impl fmt::Display for TestDataError {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            match self {
                Self::Io(err) => write!(f, "{err}"),
                Self::Json { line_no, source } => {
                    write!(f, "invalid jsonl at line {line_no}: {source}")
                }
            }
        }
    }

    impl std::error::Error for TestDataError {}

    impl From<std::io::Error> for TestDataError {
        fn from(value: std::io::Error) -> Self {
            Self::Io(value)
        }
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

    fn load_smiles() -> Result<Vec<String>, TestDataError> {
        let path = repo_root().join("tests/molblock_minimal_smiles.txt");
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

    fn ensure_golden_exists(golden_path: &Path) {
        if golden_path.exists() {
            return;
        }
        let repo = repo_root();
        let script = repo.join("tests/scripts/gen_rdkit_v2000_minimal_golden.py");
        let input = repo.join("tests/molblock_minimal_smiles.txt");

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
             uv sync --group dev && .venv/bin/python tests/scripts/gen_rdkit_v2000_minimal_golden.py --input tests/molblock_minimal_smiles.txt --output tests/golden/molblock_v2000_minimal.jsonl",
            golden_path.display(),
            last_error
        );
    }

    #[test]
    fn molblock_minimal_golden_has_entry_for_each_smiles() {
        let smiles = load_smiles().expect("read tests/molblock_minimal_smiles.txt");
        let golden =
            load_golden().expect("read tests/golden/molblock_v2000_minimal.jsonl");
        assert_eq!(
            golden.len(),
            smiles.len(),
            "golden rows must match input smiles rows"
        );

        for (idx, (record, input_smiles)) in golden.iter().zip(smiles.iter()).enumerate() {
            assert_eq!(
                record.smiles, *input_smiles,
                "smiles mismatch at row {}",
                idx + 1
            );
            if record.rdkit_ok {
                assert!(
                    record.body.is_some(),
                    "rdkit_ok=true requires body at row {}",
                    idx + 1
                );
                assert!(
                    record.error.is_none(),
                    "rdkit_ok=true should not carry error at row {}",
                    idx + 1
                );
            } else {
                assert!(
                    record.body.is_none(),
                    "rdkit_ok=false should not carry body at row {}",
                    idx + 1
                );
                assert!(
                    record.error.is_some(),
                    "rdkit_ok=false should carry error at row {}",
                    idx + 1
                );
            }
        }
    }

    #[test]
    fn molblock_minimal_matches_rdkit_golden() {
        let golden =
            load_golden().expect("read tests/golden/molblock_v2000_minimal.jsonl");
        for (idx, record) in golden.iter().enumerate() {
            let mol = Molecule::from_smiles(&record.smiles)
                .unwrap_or_else(|e| panic!("parse failed at row {}: {}", idx + 1, e));
            let ours = mol_to_v2000_block_minimal(&mol)
                .unwrap_or_else(|e| panic!("write failed at row {}: {}", idx + 1, e));
            let ours_body = body(&ours);
            let expected = record
                .body
                .as_ref()
                .expect("rdkit_ok=true requires body")
                .to_owned();
            assert_eq!(
                normalize_signed_zero(&ours_body),
                normalize_signed_zero(&expected),
                "molblock body mismatch at row {} ({})",
                idx + 1,
                record.smiles
            );
        }
    }
}
