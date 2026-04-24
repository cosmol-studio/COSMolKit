use core::fmt;
use std::collections::VecDeque;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::process::Command;

use cosmolkit_chem_core::{
    BondOrder, Molecule, ValenceModel, add_hydrogens_in_place, assign_radicals_rdkit_2025,
    assign_valence,
};
use serde::Deserialize;

#[derive(Debug, Deserialize)]
struct AtomFeature {
    atomic_num: u8,
    chirality: String,
    degree: usize,
    formal_charge: i8,
    num_hs: usize,
    num_radical_electrons: usize,
    hybridization: String,
    is_aromatic: bool,
    is_in_ring: bool,
    explicit_valence: i32,
    implicit_hs: i32,
    total_valence: i32,
}

#[derive(Debug, Deserialize)]
struct BondFeature {
    begin_atom: usize,
    end_atom: usize,
    bond_type: String,
    stereo: String,
    is_conjugated: bool,
}

#[derive(Debug, Deserialize)]
struct FeatureSet {
    atom_features: Vec<AtomFeature>,
    bond_features: Vec<BondFeature>,
}

#[derive(Debug, Deserialize)]
struct GoldenRecord {
    smiles: String,
    rdkit_ok: bool,
    direct: Option<FeatureSet>,
    with_hs: Option<FeatureSet>,
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

fn load_smiles() -> Result<Vec<String>, TestDataError> {
    let path = repo_root().join("tests/smiles.txt");
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
    let path = repo_root().join("tests/golden/atomic_nums.jsonl");
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

fn ensure_golden_exists(golden_path: &PathBuf) {
    if golden_path.exists() {
        return;
    }

    let repo = repo_root();
    let script = repo.join("tests/scripts/gen_rdkit_atomic_nums.py");
    let smiles = repo.join("tests/smiles.txt");

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
            .arg(&smiles)
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
         uv sync --group dev && .venv/bin/python tests/scripts/gen_rdkit_atomic_nums.py --input tests/smiles.txt --output tests/golden/atomic_nums.jsonl",
        golden_path.display(),
        last_error
    );
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct OursAtomFeature {
    atomic_num: u8,
    chirality: String,
    degree: usize,
    formal_charge: i8,
    num_hs: usize,
    num_radical_electrons: usize,
    hybridization: String,
    is_aromatic: bool,
    is_in_ring: bool,
    explicit_valence: i32,
    implicit_hs: i32,
    total_valence: i32,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct OursBondFeature {
    begin_atom: usize,
    end_atom: usize,
    bond_type: String,
    stereo: String,
    is_conjugated: bool,
}

fn bond_type_name(order: BondOrder) -> String {
    match order {
        BondOrder::Single => "SINGLE".to_string(),
        BondOrder::Double => "DOUBLE".to_string(),
        BondOrder::Triple => "TRIPLE".to_string(),
        BondOrder::Quadruple => "QUADRUPLE".to_string(),
        BondOrder::Aromatic => "AROMATIC".to_string(),
        BondOrder::Dative => "DATIVE".to_string(),
        BondOrder::Null => "UNSPECIFIED".to_string(),
    }
}

fn compute_ring_flags(mol: &Molecule) -> Vec<bool> {
    let n = mol.atoms.len();
    let mut adj = vec![Vec::<usize>::new(); n];
    for b in &mol.bonds {
        adj[b.begin_atom].push(b.end_atom);
        adj[b.end_atom].push(b.begin_atom);
    }
    let mut in_ring = vec![false; n];

    for a in 0..n {
        if adj[a].len() < 2 {
            continue;
        }
        let neigh = &adj[a];
        let mut found = false;
        for i in 0..neigh.len() {
            for j in (i + 1)..neigh.len() {
                let src = neigh[i];
                let dst = neigh[j];
                let mut q = VecDeque::new();
                let mut seen = vec![false; n];
                seen[a] = true;
                seen[src] = true;
                q.push_back(src);
                while let Some(v) = q.pop_front() {
                    if v == dst {
                        found = true;
                        break;
                    }
                    for &nb in &adj[v] {
                        if !seen[nb] {
                            seen[nb] = true;
                            q.push_back(nb);
                        }
                    }
                }
                if found {
                    break;
                }
            }
            if found {
                break;
            }
        }
        in_ring[a] = found;
    }
    in_ring
}

fn extract_ours_features(mol: &Molecule) -> (Vec<OursAtomFeature>, Vec<OursBondFeature>) {
    let ring_flags = compute_ring_flags(mol);
    let assignment = assign_valence(mol, ValenceModel::RdkitLike).unwrap_or_else(|e| {
        panic!("assign_valence failed in test extraction: {:?}", e);
    });
    let radicals =
        assign_radicals_rdkit_2025(mol, &assignment.explicit_valence).unwrap_or_else(|e| {
            panic!(
                "assign_radicals_rdkit_2025 failed in test extraction: {:?}",
                e
            );
        });
    let mut atom_degree = vec![0usize; mol.atoms.len()];
    let mut atom_has_multi = vec![false; mol.atoms.len()];

    for b in &mol.bonds {
        atom_degree[b.begin_atom] += 1;
        atom_degree[b.end_atom] += 1;
        if matches!(
            b.order,
            BondOrder::Double | BondOrder::Triple | BondOrder::Quadruple | BondOrder::Aromatic
        ) {
            atom_has_multi[b.begin_atom] = true;
            atom_has_multi[b.end_atom] = true;
        }
    }

    let mut atoms = Vec::with_capacity(mol.atoms.len());
    for (i, a) in mol.atoms.iter().enumerate() {
        let implicit_hs = assignment.implicit_hydrogens[i] as i32;
        let explicit_valence = assignment.explicit_valence[i] as i32;
        let num_hs = (a.explicit_hydrogens as i32 + implicit_hs).max(0) as usize;
        let degree = atom_degree[i] + num_hs;
        let total_valence = explicit_valence + implicit_hs;
        let hybridization = if a.is_aromatic {
            "SP2"
        } else if explicit_valence >= 3 && atom_has_multi[i] {
            "SP"
        } else if atom_has_multi[i] {
            "SP2"
        } else {
            "SP3"
        };
        atoms.push(OursAtomFeature {
            atomic_num: a.atomic_num,
            chirality: "CHI_UNSPECIFIED".to_string(),
            degree,
            formal_charge: a.formal_charge,
            num_hs,
            num_radical_electrons: radicals[i] as usize,
            hybridization: hybridization.to_string(),
            is_aromatic: a.is_aromatic,
            is_in_ring: ring_flags[i],
            explicit_valence,
            implicit_hs,
            total_valence,
        });
    }

    let mut bonds = Vec::with_capacity(mol.bonds.len());
    for b in &mol.bonds {
        bonds.push(OursBondFeature {
            begin_atom: b.begin_atom,
            end_atom: b.end_atom,
            bond_type: bond_type_name(b.order),
            stereo: "STEREONONE".to_string(),
            is_conjugated: matches!(b.order, BondOrder::Aromatic | BondOrder::Double),
        });
    }
    (atoms, bonds)
}

fn compare_features(
    ours_atoms: &[OursAtomFeature],
    ours_bonds: &[OursBondFeature],
    expected: &FeatureSet,
    row_idx: usize,
    smiles: &str,
    mode: &str,
) {
    assert_eq!(
        ours_atoms.len(),
        expected.atom_features.len(),
        "atom count mismatch at row {} ({}) [{}]",
        row_idx,
        smiles,
        mode
    );
    assert_eq!(
        ours_bonds.len(),
        expected.bond_features.len(),
        "bond count mismatch at row {} ({}) [{}]",
        row_idx,
        smiles,
        mode
    );

    for i in 0..ours_atoms.len() {
        let a = &ours_atoms[i];
        let e = &expected.atom_features[i];
        assert_eq!(
            a.atomic_num, e.atomic_num,
            "atomic_num mismatch row {} atom {} ({}) [{}]",
            row_idx, i, smiles, mode
        );
        assert_eq!(
            a.degree, e.degree,
            "degree mismatch row {} atom {} ({}) [{}]",
            row_idx, i, smiles, mode
        );
        assert_eq!(
            a.formal_charge, e.formal_charge,
            "formal_charge mismatch row {} atom {} ({}) [{}]",
            row_idx, i, smiles, mode
        );
        assert_eq!(
            a.num_hs, e.num_hs,
            "num_hs mismatch row {} atom {} ({}) [{}]",
            row_idx, i, smiles, mode
        );
        assert_eq!(
            a.num_radical_electrons, e.num_radical_electrons,
            "num_radical_electrons mismatch row {} atom {} ({}) [{}]",
            row_idx, i, smiles, mode
        );
        assert_eq!(
            a.is_aromatic, e.is_aromatic,
            "is_aromatic mismatch row {} atom {} ({}) [{}]",
            row_idx, i, smiles, mode
        );
        assert_eq!(
            a.is_in_ring, e.is_in_ring,
            "is_in_ring mismatch row {} atom {} ({}) [{}]",
            row_idx, i, smiles, mode
        );
        assert_eq!(
            a.explicit_valence, e.explicit_valence,
            "explicit_valence mismatch row {} atom {} ({}) [{}]",
            row_idx, i, smiles, mode
        );
        assert_eq!(
            a.implicit_hs, e.implicit_hs,
            "implicit_hs mismatch row {} atom {} ({}) [{}]",
            row_idx, i, smiles, mode
        );
        assert_eq!(
            a.total_valence, e.total_valence,
            "total_valence mismatch row {} atom {} ({}) [{}]",
            row_idx, i, smiles, mode
        );

        // currently unsupported-by-model fields are still part of the golden:
        let _ = (&a.chirality, &e.chirality);
        let _ = (&a.hybridization, &e.hybridization);
    }

    let normalize = |mut b: OursBondFeature| -> OursBondFeature {
        if b.bond_type != "DATIVE" && b.begin_atom > b.end_atom {
            std::mem::swap(&mut b.begin_atom, &mut b.end_atom);
        }
        b
    };
    let mut ours_bonds_sorted: Vec<OursBondFeature> =
        ours_bonds.iter().cloned().map(normalize).collect();
    ours_bonds_sorted.sort_by(|l, r| {
        (
            l.begin_atom,
            l.end_atom,
            &l.bond_type,
            &l.stereo,
            l.is_conjugated,
        )
            .cmp(&(
                r.begin_atom,
                r.end_atom,
                &r.bond_type,
                &r.stereo,
                r.is_conjugated,
            ))
    });
    let mut expected_bonds_sorted: Vec<OursBondFeature> = expected
        .bond_features
        .iter()
        .map(|b| OursBondFeature {
            begin_atom: b.begin_atom,
            end_atom: b.end_atom,
            bond_type: b.bond_type.clone(),
            stereo: b.stereo.clone(),
            is_conjugated: b.is_conjugated,
        })
        .map(normalize)
        .collect();
    expected_bonds_sorted.sort_by(|l, r| {
        (
            l.begin_atom,
            l.end_atom,
            &l.bond_type,
            &l.stereo,
            l.is_conjugated,
        )
            .cmp(&(
                r.begin_atom,
                r.end_atom,
                &r.bond_type,
                &r.stereo,
                r.is_conjugated,
            ))
    });

    for i in 0..ours_bonds_sorted.len() {
        let b = &ours_bonds_sorted[i];
        let e = &expected_bonds_sorted[i];
        assert_eq!(
            b.begin_atom, e.begin_atom,
            "bond begin mismatch row {} bond {} ({}) [{}]",
            row_idx, i, smiles, mode
        );
        assert_eq!(
            b.end_atom, e.end_atom,
            "bond end mismatch row {} bond {} ({}) [{}]",
            row_idx, i, smiles, mode
        );
        assert_eq!(
            b.bond_type, e.bond_type,
            "bond type mismatch row {} bond {} ({}) [{}]",
            row_idx, i, smiles, mode
        );

        // currently unsupported-by-model field:
        let _ = (&b.stereo, &e.stereo, b.is_conjugated, e.is_conjugated);
    }
}

#[test]
fn golden_has_entry_for_each_smiles_input() {
    let smiles = load_smiles().expect("should read tests/smiles.txt");
    let golden = load_golden().expect("should read tests/golden/atomic_nums.jsonl");
    assert_eq!(
        golden.len(),
        smiles.len(),
        "golden rows must match input smiles rows"
    );

    for (idx, (record, raw_smiles)) in golden.iter().zip(smiles.iter()).enumerate() {
        assert_eq!(
            record.smiles,
            *raw_smiles,
            "smiles mismatch at row {}",
            idx + 1
        );

        if record.rdkit_ok {
            assert!(
                record.direct.is_some() && record.with_hs.is_some(),
                "rdkit_ok=true requires both direct and with_hs at row {}",
                idx + 1
            );
            assert!(
                record.error.is_none(),
                "rdkit_ok=true should not carry error at row {}",
                idx + 1
            );
        } else {
            assert!(
                record.direct.is_none() && record.with_hs.is_none(),
                "rdkit_ok=false should not carry feature sets at row {}",
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
fn graph_features_match_rdkit_golden_dual_track() {
    let golden = load_golden().expect("should read tests/golden/atomic_nums.jsonl");

    for (idx, record) in golden.iter().enumerate() {
        let ours = Molecule::from_smiles(&record.smiles);
        match (record.rdkit_ok, ours) {
            (true, Ok(mol_direct)) => {
                let direct_expected = record.direct.as_ref().expect("direct missing");
                let (ours_atoms_direct, ours_bonds_direct) = extract_ours_features(&mol_direct);
                compare_features(
                    &ours_atoms_direct,
                    &ours_bonds_direct,
                    direct_expected,
                    idx + 1,
                    &record.smiles,
                    "direct",
                );

                let with_h_expected = record.with_hs.as_ref().expect("with_hs missing");
                let mut mol_with_h = mol_direct.clone();
                add_hydrogens_in_place(&mut mol_with_h).unwrap_or_else(|e| {
                    panic!(
                        "add_hydrogens failed at row {} ({}): {:?}",
                        idx + 1,
                        record.smiles,
                        e
                    )
                });
                let (ours_atoms_h, ours_bonds_h) = extract_ours_features(&mol_with_h);
                compare_features(
                    &ours_atoms_h,
                    &ours_bonds_h,
                    with_h_expected,
                    idx + 1,
                    &record.smiles,
                    "with_hs",
                );
            }
            (false, Err(_)) => {}
            (true, Err(err)) => {
                panic!(
                    "unexpected parse error at row {} ({}): {}",
                    idx + 1,
                    record.smiles,
                    err
                );
            }
            (false, Ok(_)) => {
                panic!(
                    "expected parse failure at row {} ({})",
                    idx + 1,
                    record.smiles
                );
            }
        }
    }
}
