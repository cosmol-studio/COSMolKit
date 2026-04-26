use core::fmt;
use std::collections::VecDeque;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;
use std::process::Command;

use cosmolkit_core::{
    BondOrder, BondStereo, ChiralTag as OursChiralTag, Molecule, ValenceModel,
    add_hydrogens_in_place, assign_radicals_rdkit_2025, assign_valence, rdkit_valence_list,
};
use serde::{Deserialize, Deserializer};

#[derive(Debug, Clone, PartialEq, Eq)]
enum ChiralTag {
    ChiUnspecified,
    ChiTetrahedralCcw,
    ChiTetrahedralCw,
    Unknown(String),
}

impl<'de> Deserialize<'de> for ChiralTag {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let value = String::deserialize(deserializer)?;
        let tag = match value.as_str() {
            "CHI_UNSPECIFIED" => Self::ChiUnspecified,
            "CHI_TETRAHEDRAL_CCW" => Self::ChiTetrahedralCcw,
            "CHI_TETRAHEDRAL_CW" => Self::ChiTetrahedralCw,
            _ => Self::Unknown(value),
        };
        Ok(tag)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum CipCode {
    R,
    S,
    LowerR,
    LowerS,
    Unknown(String),
}

impl<'de> Deserialize<'de> for CipCode {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let value = String::deserialize(deserializer)?;
        let code = match value.as_str() {
            "R" => Self::R,
            "S" => Self::S,
            "r" => Self::LowerR,
            "s" => Self::LowerS,
            _ => Self::Unknown(value),
        };
        Ok(code)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum Hybridization {
    S,
    Sp,
    Sp2,
    Sp3,
    Sp2D,
    Sp3D,
    Sp3D2,
    Other,
    Unspecified,
    Unknown(String),
}

impl<'de> Deserialize<'de> for Hybridization {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let value = String::deserialize(deserializer)?;
        let hybridization = match value.as_str() {
            "S" => Self::S,
            "SP" => Self::Sp,
            "SP2" => Self::Sp2,
            "SP3" => Self::Sp3,
            "SP2D" => Self::Sp2D,
            "SP3D" => Self::Sp3D,
            "SP3D2" => Self::Sp3D2,
            "OTHER" => Self::Other,
            "UNSPECIFIED" => Self::Unspecified,
            _ => Self::Unknown(value),
        };
        Ok(hybridization)
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
enum RdkitBondType {
    Unspecified,
    Single,
    Double,
    Triple,
    Quadruple,
    Quintuple,
    Hextuple,
    OneAndAHalf,
    TwoAndAHalf,
    ThreeAndAHalf,
    FourAndAHalf,
    FiveAndAHalf,
    Aromatic,
    Ionic,
    Hydrogen,
    ThreeCenter,
    DativeOne,
    Dative,
    DativeL,
    DativeR,
    Other,
    Zero,
    Unknown(String),
}

impl<'de> Deserialize<'de> for RdkitBondType {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        // Mirrors rdkit::bondTypeToString() naming in Code/GraphMol/Bond.cpp.
        let value = String::deserialize(deserializer)?;
        let bond_type = match value.as_str() {
            "UNSPECIFIED" => Self::Unspecified,
            "SINGLE" => Self::Single,
            "DOUBLE" => Self::Double,
            "TRIPLE" => Self::Triple,
            "QUADRUPLE" => Self::Quadruple,
            "QUINTUPLE" => Self::Quintuple,
            "HEXTUPLE" => Self::Hextuple,
            "ONEANDAHALF" => Self::OneAndAHalf,
            "TWOANDAHALF" => Self::TwoAndAHalf,
            "THREEANDAHALF" => Self::ThreeAndAHalf,
            "FOURANDAHALF" => Self::FourAndAHalf,
            "FIVEANDAHALF" => Self::FiveAndAHalf,
            "AROMATIC" => Self::Aromatic,
            "IONIC" => Self::Ionic,
            "HYDROGEN" => Self::Hydrogen,
            "THREECENTER" => Self::ThreeCenter,
            "DATIVEONE" => Self::DativeOne,
            "DATIVE" => Self::Dative,
            "DATIVEL" => Self::DativeL,
            "DATIVER" => Self::DativeR,
            "OTHER" => Self::Other,
            "ZERO" => Self::Zero,
            _ => Self::Unknown(value),
        };
        Ok(bond_type)
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
enum RdkitBondStereo {
    StereoNone,
    StereoAny,
    StereoZ,
    StereoE,
    StereoCis,
    StereoTrans,
    StereoAtropCw,
    StereoAtropCcw,
    Unknown(String),
}

impl<'de> Deserialize<'de> for RdkitBondStereo {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        // Mirrors rdkit::bondStereoToString() naming in Code/GraphMol/Bond.cpp.
        let value = String::deserialize(deserializer)?;
        let stereo = match value.as_str() {
            "STEREONONE" => Self::StereoNone,
            "STEREOANY" => Self::StereoAny,
            "STEREOZ" => Self::StereoZ,
            "STEREOE" => Self::StereoE,
            "STEREOCIS" => Self::StereoCis,
            "STEREOTRANS" => Self::StereoTrans,
            "STEREOATROPCW" => Self::StereoAtropCw,
            "STEREOATROPCCW" => Self::StereoAtropCcw,
            _ => Self::Unknown(value),
        };
        Ok(stereo)
    }
}

#[derive(Debug, Deserialize)]
struct AtomFeature {
    atomic_num: u8,
    chirality: ChiralTag,
    cip_code: Option<CipCode>,
    cip_rank: Option<i64>,
    chirality_possible: bool,
    degree: usize,
    formal_charge: i8,
    num_hs: usize,
    num_radical_electrons: usize,
    hybridization: Hybridization,
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
    bond_type: RdkitBondType,
    stereo: RdkitBondStereo,
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
    let path = repo_root().join("tests/golden/graph_features.jsonl");
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
    let script = repo.join("tests/scripts/gen_rdkit_graph_features.py");
    let smiles = repo.join("tests/smiles.smi");

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
         uv sync --group dev && .venv/bin/python tests/scripts/gen_rdkit_graph_features.py --input tests/smiles.smi --output tests/golden/graph_features.jsonl",
        golden_path.display(),
        last_error
    );
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct OursAtomFeature {
    atomic_num: u8,
    chirality: ChiralTag,
    degree: usize,
    formal_charge: i8,
    num_hs: usize,
    num_radical_electrons: usize,
    hybridization: Hybridization,
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
    bond_type: RdkitBondType,
    stereo: RdkitBondStereo,
    is_conjugated: bool,
}

fn bond_type_name(order: BondOrder) -> RdkitBondType {
    match order {
        BondOrder::Single => RdkitBondType::Single,
        BondOrder::Double => RdkitBondType::Double,
        BondOrder::Triple => RdkitBondType::Triple,
        BondOrder::Quadruple => RdkitBondType::Quadruple,
        BondOrder::Aromatic => RdkitBondType::Aromatic,
        BondOrder::Dative => RdkitBondType::Dative,
        BondOrder::Null => RdkitBondType::Unspecified,
    }
}

fn ours_chiral_tag_name(tag: OursChiralTag) -> ChiralTag {
    match tag {
        OursChiralTag::Unspecified => ChiralTag::ChiUnspecified,
        OursChiralTag::TetrahedralCw => ChiralTag::ChiTetrahedralCw,
        OursChiralTag::TetrahedralCcw => ChiralTag::ChiTetrahedralCcw,
    }
}

fn ours_bond_stereo_name(stereo: BondStereo) -> RdkitBondStereo {
    match stereo {
        BondStereo::None => RdkitBondStereo::StereoNone,
        BondStereo::Any => RdkitBondStereo::StereoAny,
        BondStereo::Cis => RdkitBondStereo::StereoZ,
        BondStereo::Trans => RdkitBondStereo::StereoE,
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

fn rdkit_default_valence(atomic_num: u8) -> Option<i32> {
    let vals = rdkit_valence_list(atomic_num)?;
    vals.iter().copied().find(|v| *v >= 0)
}

fn rdkit_n_outer_electrons(atomic_num: u8) -> Option<i32> {
    match atomic_num {
        0 => Some(0),
        1 => Some(1),
        2 => Some(2),
        3 => Some(1),
        4 => Some(2),
        5 => Some(3),
        6 => Some(4),
        7 => Some(5),
        8 => Some(6),
        9 => Some(7),
        10 => Some(8),
        11 => Some(1),
        12 => Some(2),
        13 => Some(3),
        14 => Some(4),
        15 => Some(5),
        16 => Some(6),
        17 => Some(7),
        18 => Some(8),
        19 => Some(1),
        20 => Some(2),
        21 => Some(3),
        22 => Some(4),
        23 => Some(5),
        24 => Some(6),
        25 => Some(7),
        26 => Some(8),
        27 => Some(9),
        28 => Some(10),
        29 => Some(11),
        30 => Some(2),
        31 => Some(3),
        32 => Some(4),
        33 => Some(5),
        34 => Some(6),
        35 => Some(7),
        36 => Some(8),
        37 => Some(1),
        38 => Some(2),
        39 => Some(3),
        40 => Some(4),
        41 => Some(5),
        42 => Some(6),
        43 => Some(7),
        44 => Some(8),
        45 => Some(9),
        46 => Some(10),
        47 => Some(11),
        48 => Some(2),
        49 => Some(3),
        50 => Some(4),
        51 => Some(5),
        52 => Some(6),
        53 => Some(7),
        54 => Some(8),
        55 => Some(1),
        56 => Some(2),
        57 => Some(3),
        58 => Some(4),
        59 => Some(3),
        60 => Some(4),
        61 => Some(5),
        62 => Some(6),
        63 => Some(7),
        64 => Some(8),
        65 => Some(9),
        66 => Some(10),
        67 => Some(11),
        68 => Some(12),
        69 => Some(13),
        70 => Some(14),
        71 => Some(15),
        72 => Some(4),
        73 => Some(5),
        74 => Some(6),
        75 => Some(7),
        76 => Some(8),
        77 => Some(9),
        78 => Some(10),
        79 => Some(11),
        80 => Some(2),
        81 => Some(3),
        82 => Some(4),
        83 => Some(5),
        84 => Some(6),
        85 => Some(7),
        86 => Some(8),
        87 => Some(1),
        88 => Some(2),
        89 => Some(3),
        90 => Some(4),
        91 => Some(3),
        92 => Some(4),
        93 => Some(5),
        94 => Some(6),
        95 => Some(7),
        96 => Some(8),
        97 => Some(9),
        98 => Some(10),
        99 => Some(11),
        100 => Some(12),
        101 => Some(13),
        102 => Some(14),
        103 => Some(15),
        104..=118 => Some(2),
        _ => None,
    }
}

fn bond_valence_contrib_for_atom(b: &cosmolkit_core::Bond, atom_index: usize) -> f64 {
    if b.begin_atom != atom_index && b.end_atom != atom_index {
        return 0.0;
    }
    match b.order {
        BondOrder::Null => 0.0,
        BondOrder::Single => 1.0,
        BondOrder::Double => 2.0,
        BondOrder::Triple => 3.0,
        BondOrder::Quadruple => 4.0,
        BondOrder::Aromatic => 1.5,
        BondOrder::Dative => {
            if b.end_atom == atom_index {
                1.0
            } else {
                0.0
            }
        }
    }
}

fn count_atom_electrons_rdkit(
    mol: &Molecule,
    assignment: &cosmolkit_core::ValenceAssignment,
    atom_degree: &[usize],
    atom_index: usize,
) -> i32 {
    let atom = &mol.atoms[atom_index];
    let Some(dv) = rdkit_default_valence(atom.atomic_num) else {
        return -1;
    };
    if dv <= 1 {
        return -1;
    }
    let mut degree = atom_degree[atom_index] as i32
        + atom.explicit_hydrogens as i32
        + assignment.implicit_hydrogens[atom_index] as i32;
    for b in &mol.bonds {
        if (b.begin_atom == atom_index || b.end_atom == atom_index)
            && bond_valence_contrib_for_atom(b, atom_index) == 0.0
        {
            degree -= 1;
        }
    }
    if degree > 3 {
        return -1;
    }
    let Some(nouter) = rdkit_n_outer_electrons(atom.atomic_num) else {
        return -1;
    };
    let nlp = (nouter - dv - atom.formal_charge as i32).max(0);
    let radicals = atom.num_radical_electrons as i32;
    let mut res = (dv - degree) + nlp - radicals;
    if res > 1 {
        let n_unsaturations =
            assignment.explicit_valence[atom_index] as i32 - atom_degree[atom_index] as i32;
        if n_unsaturations > 1 {
            res = 1;
        }
    }
    res
}

fn is_atom_conjug_cand(
    mol: &Molecule,
    assignment: &cosmolkit_core::ValenceAssignment,
    atom_degree: &[usize],
    atom_index: usize,
) -> bool {
    let at = &mol.atoms[atom_index];
    if let Some(vals) = rdkit_valence_list(at.atomic_num)
        && at.formal_charge == 0
        && !vals.is_empty()
        && vals[0] >= 0
    {
        let total_valence = assignment.explicit_valence[atom_index] as i32
            + assignment.implicit_hydrogens[atom_index] as i32;
        if total_valence > vals[0] {
            return false;
        }
    }
    let nouter = rdkit_n_outer_electrons(at.atomic_num).unwrap_or(0);
    let row_ok = at.atomic_num <= 10
        || (nouter != 5 && nouter != 6)
        || (nouter == 6 && atom_degree[atom_index] < 2);
    row_ok && count_atom_electrons_rdkit(mol, assignment, atom_degree, atom_index) > 0
}

fn compute_conjugated_bonds(
    mol: &Molecule,
    assignment: &cosmolkit_core::ValenceAssignment,
    atom_degree: &[usize],
) -> Vec<bool> {
    let mut conjugated = vec![false; mol.bonds.len()];
    for (bi, b) in mol.bonds.iter().enumerate() {
        conjugated[bi] = matches!(b.order, BondOrder::Aromatic);
    }
    for at in 0..mol.atoms.len() {
        if !is_atom_conjug_cand(mol, assignment, atom_degree, at) {
            continue;
        }
        let sbo = atom_degree[at]
            + mol.atoms[at].explicit_hydrogens as usize
            + assignment.implicit_hydrogens[at] as usize;
        if !(2..=3).contains(&sbo) {
            continue;
        }
        let bnds: Vec<usize> = mol
            .bonds
            .iter()
            .enumerate()
            .filter_map(|(bi, b)| {
                if b.begin_atom == at || b.end_atom == at {
                    Some(bi)
                } else {
                    None
                }
            })
            .collect();
        for &b1 in &bnds {
            let bond1 = &mol.bonds[b1];
            if bond_valence_contrib_for_atom(bond1, at) < 1.5 {
                continue;
            }
            let o1 = if bond1.begin_atom == at {
                bond1.end_atom
            } else {
                bond1.begin_atom
            };
            if !is_atom_conjug_cand(mol, assignment, atom_degree, o1) {
                continue;
            }
            for &b2 in &bnds {
                if b1 == b2 {
                    continue;
                }
                let bond2 = &mol.bonds[b2];
                let o2 = if bond2.begin_atom == at {
                    bond2.end_atom
                } else {
                    bond2.begin_atom
                };
                let sbo2 = atom_degree[o2]
                    + mol.atoms[o2].explicit_hydrogens as usize
                    + assignment.implicit_hydrogens[o2] as usize;
                if sbo2 > 3 {
                    continue;
                }
                if is_atom_conjug_cand(mol, assignment, atom_degree, o2) {
                    conjugated[b1] = true;
                    conjugated[b2] = true;
                }
            }
        }
    }
    conjugated
}

fn compute_hybridization(
    mol: &Molecule,
    assignment: &cosmolkit_core::ValenceAssignment,
    atom_degree: &[usize],
    atom_has_conjugated_bond: &[bool],
    atom_index: usize,
) -> Hybridization {
    let atom = &mol.atoms[atom_index];
    if atom.atomic_num == 1 {
        // RDKit parity corpus behavior: hydrogens introduced by AddHs() stay
        // UNSPECIFIED, while bracket/input hydrogens (including isotopes) are S.
        if !atom.no_implicit && atom.isotope.is_none() {
            return Hybridization::Unspecified;
        }
    }
    if atom.atomic_num == 0 {
        return Hybridization::Unspecified;
    }
    // RDKit shortcut for tetrahedral stereocenters.
    if !matches!(atom.chiral_tag, OursChiralTag::Unspecified)
        && atom_degree[atom_index]
            + atom.explicit_hydrogens as usize
            + assignment.implicit_hydrogens[atom_index] as usize
            == 4
    {
        return Hybridization::Sp3;
    }
    let mut deg = atom_degree[atom_index] as i32
        + atom.explicit_hydrogens as i32
        + assignment.implicit_hydrogens[atom_index] as i32;
    for b in &mol.bonds {
        if (b.begin_atom == atom_index || b.end_atom == atom_index)
            && (matches!(b.order, BondOrder::Dative) && b.end_atom != atom_index)
        {
            deg -= 1;
        }
    }
    if atom.atomic_num <= 1 {
        return match deg {
            0 | 1 => Hybridization::S,
            2 => Hybridization::Sp,
            3 => Hybridization::Sp2,
            4 => Hybridization::Sp3,
            5 => Hybridization::Sp3D,
            6 => Hybridization::Sp3D2,
            _ => Hybridization::Unspecified,
        };
    }
    let nouter = rdkit_n_outer_electrons(atom.atomic_num).unwrap_or(0);
    let total_valence = assignment.explicit_valence[atom_index] as i32
        + assignment.implicit_hydrogens[atom_index] as i32;
    let num_free = nouter - (total_valence + atom.formal_charge as i32);
    let norbs = if total_valence + nouter - (atom.formal_charge as i32) < 8 {
        let radicals = atom.num_radical_electrons as i32;
        let lone_pairs = (num_free - radicals) / 2;
        deg + lone_pairs + radicals
    } else {
        let lone_pairs = num_free / 2;
        deg + lone_pairs
    };
    match norbs {
        0 | 1 => Hybridization::S,
        2 => Hybridization::Sp,
        3 => Hybridization::Sp2,
        4 => {
            let total_degree = atom_degree[atom_index]
                + atom.explicit_hydrogens as usize
                + assignment.implicit_hydrogens[atom_index] as usize;
            if total_degree > 3 || !atom_has_conjugated_bond[atom_index] {
                Hybridization::Sp3
            } else {
                Hybridization::Sp2
            }
        }
        5 => Hybridization::Sp3D,
        6 => Hybridization::Sp3D2,
        _ => Hybridization::Unspecified,
    }
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

    let conjugated_bonds = compute_conjugated_bonds(mol, &assignment, &atom_degree);
    let mut atom_has_conjugated_bond = vec![false; mol.atoms.len()];
    for (bi, b) in mol.bonds.iter().enumerate() {
        if conjugated_bonds[bi] {
            atom_has_conjugated_bond[b.begin_atom] = true;
            atom_has_conjugated_bond[b.end_atom] = true;
        }
    }

    let mut atoms = Vec::with_capacity(mol.atoms.len());
    for (i, a) in mol.atoms.iter().enumerate() {
        let implicit_hs = assignment.implicit_hydrogens[i] as i32;
        let explicit_valence = assignment.explicit_valence[i] as i32;
        let num_hs = (a.explicit_hydrogens as i32 + implicit_hs).max(0) as usize;
        let degree = atom_degree[i] + num_hs;
        let total_valence = explicit_valence + implicit_hs;
        let _ = atom_has_multi[i];
        let hybridization =
            compute_hybridization(mol, &assignment, &atom_degree, &atom_has_conjugated_bond, i);
        atoms.push(OursAtomFeature {
            atomic_num: a.atomic_num,
            chirality: ours_chiral_tag_name(a.chiral_tag),
            degree,
            formal_charge: a.formal_charge,
            num_hs,
            num_radical_electrons: radicals[i] as usize,
            hybridization,
            is_aromatic: a.is_aromatic,
            is_in_ring: ring_flags[i],
            explicit_valence,
            implicit_hs,
            total_valence,
        });
    }

    let mut bonds = Vec::with_capacity(mol.bonds.len());
    for (bi, b) in mol.bonds.iter().enumerate() {
        bonds.push(OursBondFeature {
            begin_atom: b.begin_atom,
            end_atom: b.end_atom,
            bond_type: bond_type_name(b.order),
            stereo: ours_bond_stereo_name(b.stereo),
            is_conjugated: conjugated_bonds[bi],
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

        assert_eq!(
            a.chirality, e.chirality,
            "chirality mismatch row {} atom {} ({}) [{}]",
            row_idx, i, smiles, mode
        );
        assert_eq!(
            a.hybridization, e.hybridization,
            "hybridization mismatch row {} atom {} ({}) [{}]",
            row_idx, i, smiles, mode
        );
    }

    let normalize = |mut b: OursBondFeature| -> OursBondFeature {
        if b.bond_type != RdkitBondType::Dative && b.begin_atom > b.end_atom {
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

        assert_eq!(
            b.stereo, e.stereo,
            "bond stereo mismatch row {} bond {} ({}) [{}]",
            row_idx, i, smiles, mode
        );
        assert_eq!(
            b.is_conjugated, e.is_conjugated,
            "bond conjugation mismatch row {} bond {} ({}) [{}]",
            row_idx, i, smiles, mode
        );
    }
}

#[test]
fn graph_feature_golden_has_one_record_per_smiles() {
    let smiles = load_smiles().expect("should read tests/smiles.smi");
    let golden = load_golden().expect("should read tests/golden/graph_features.jsonl");
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
fn graph_feature_golden_records_cip_for_chiral_atoms() {
    let golden = load_golden().expect("should read tests/golden/graph_features.jsonl");
    let mut chiral_atoms = 0usize;
    let mut cip_labeled_chiral_atoms = 0usize;

    for (idx, record) in golden.iter().enumerate() {
        if !record.rdkit_ok {
            continue;
        }
        let direct = record.direct.as_ref().expect("direct missing");
        for (atom_idx, atom) in direct.atom_features.iter().enumerate() {
            if matches!(atom.chirality, ChiralTag::ChiUnspecified) {
                continue;
            }
            chiral_atoms += 1;
            assert!(
                atom.cip_rank.is_some(),
                "chiral atom must carry RDKit _CIPRank at row {} atom {} ({})",
                idx + 1,
                atom_idx,
                record.smiles
            );
            assert!(
                matches!(
                    atom.cip_code,
                    Some(CipCode::R | CipCode::S | CipCode::LowerR | CipCode::LowerS)
                ),
                "chiral atom must carry RDKit _CIPCode at row {} atom {} ({})",
                idx + 1,
                atom_idx,
                record.smiles
            );
            assert!(
                !atom.chirality_possible,
                "chiral atom should not carry RDKit _ChiralityPossible in this corpus at row {} atom {} ({})",
                idx + 1,
                atom_idx,
                record.smiles
            );
            cip_labeled_chiral_atoms += 1;
        }
    }

    assert!(
        chiral_atoms >= 8,
        "SMILES corpus should include multiple RDKit-recognized chiral atoms"
    );
    assert_eq!(cip_labeled_chiral_atoms, chiral_atoms);
}

#[test]
fn graph_features_match_rdkit_golden_for_direct_and_explicit_hydrogen_molecules() {
    let golden = load_golden().expect("should read tests/golden/graph_features.jsonl");

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
