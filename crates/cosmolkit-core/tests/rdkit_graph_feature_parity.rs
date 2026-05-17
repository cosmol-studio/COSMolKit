use std::collections::VecDeque;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

use cosmolkit_core::{
    BatchErrorMode, BatchRecord, BondDirection, BondOrder, BondStereo, ChiralTag as OursChiralTag,
    Molecule, MoleculeBatch, ValenceModel, assign_valence, rdkit_valence_list,
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
    isotope: Option<u16>,
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
    addhs_removehs: Option<FeatureSet>,
    possible_stereo: Option<FeatureSet>,
    chiral_centers: Option<ChiralCentersRecord>,
    error: Option<String>,
}

#[derive(Debug, Deserialize)]
struct ChiralCenter {
    atom_idx: usize,
    label: String,
}

#[derive(Debug, Deserialize)]
struct ChiralCentersRecord {
    include_unassigned_false: Vec<ChiralCenter>,
    include_unassigned_true: Vec<ChiralCenter>,
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
    assert!(
        golden_path.exists(),
        "missing RDKit graph feature golden: {}. Generate it before running tests:\n\
         uv sync --group dev && .venv/bin/python tests/scripts/gen_rdkit_graph_features.py --input tests/smiles.smi --output tests/golden/graph_features.jsonl",
        golden_path.display()
    );
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct OursAtomFeature {
    atomic_num: u8,
    isotope: Option<u16>,
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
struct OursPossibleStereoAtomFeature {
    base: OursAtomFeature,
    cip_code: Option<CipCode>,
    cip_rank: Option<i64>,
    chirality_possible: bool,
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
        BondOrder::Unspecified | BondOrder::Null => RdkitBondType::Unspecified,
        BondOrder::Single => RdkitBondType::Single,
        BondOrder::Double => RdkitBondType::Double,
        BondOrder::Triple => RdkitBondType::Triple,
        BondOrder::Quadruple => RdkitBondType::Quadruple,
        BondOrder::Quintuple => RdkitBondType::Quintuple,
        BondOrder::Hextuple => RdkitBondType::Hextuple,
        BondOrder::OneAndHalf => RdkitBondType::OneAndAHalf,
        BondOrder::TwoAndHalf => RdkitBondType::TwoAndAHalf,
        BondOrder::ThreeAndHalf => RdkitBondType::ThreeAndAHalf,
        BondOrder::FourAndHalf => RdkitBondType::FourAndAHalf,
        BondOrder::FiveAndHalf => RdkitBondType::FiveAndAHalf,
        BondOrder::Aromatic => RdkitBondType::Aromatic,
        BondOrder::Ionic => RdkitBondType::Ionic,
        BondOrder::ThreeCenter => RdkitBondType::ThreeCenter,
        BondOrder::DativeOne => RdkitBondType::DativeOne,
        BondOrder::Dative => RdkitBondType::Dative,
        BondOrder::DativeLeft => RdkitBondType::DativeL,
        BondOrder::DativeRight => RdkitBondType::DativeR,
        BondOrder::Hydrogen => RdkitBondType::Hydrogen,
        BondOrder::Other => RdkitBondType::Other,
        BondOrder::Zero => RdkitBondType::Zero,
    }
}

fn ours_chiral_tag_name(tag: OursChiralTag) -> ChiralTag {
    match tag {
        OursChiralTag::Unspecified => ChiralTag::ChiUnspecified,
        OursChiralTag::TetrahedralCw => ChiralTag::ChiTetrahedralCw,
        OursChiralTag::TetrahedralCcw => ChiralTag::ChiTetrahedralCcw,
        other => ChiralTag::Unknown(other.rdkit_name().to_owned()),
    }
}

fn ours_bond_stereo_name(stereo: BondStereo) -> RdkitBondStereo {
    match stereo {
        BondStereo::None => RdkitBondStereo::StereoNone,
        BondStereo::Any => RdkitBondStereo::StereoAny,
        BondStereo::Z | BondStereo::Cis => RdkitBondStereo::StereoZ,
        BondStereo::E | BondStereo::Trans => RdkitBondStereo::StereoE,
        BondStereo::AtropCw => RdkitBondStereo::StereoAtropCw,
        BondStereo::AtropCcw => RdkitBondStereo::StereoAtropCcw,
    }
}

fn ours_cip_code_name(code: Option<&str>) -> Option<CipCode> {
    match code {
        Some("R") => Some(CipCode::R),
        Some("S") => Some(CipCode::S),
        Some("r") => Some(CipCode::LowerR),
        Some("s") => Some(CipCode::LowerS),
        Some(other) => Some(CipCode::Unknown(other.to_owned())),
        None => None,
    }
}

fn compute_ring_flags(mol: &Molecule) -> Vec<bool> {
    let n = mol.atoms().len();
    let mut adj = vec![Vec::<usize>::new(); n];
    for b in mol.bonds() {
        if matches!(b.order(), BondOrder::Null | BondOrder::Dative) {
            continue;
        }
        adj[b.begin().index()].push(b.end().index());
        adj[b.end().index()].push(b.begin().index());
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
    let vals = rdkit_valence_list(atomic_num).ok()??;
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
    if b.begin().index() != atom_index && b.end().index() != atom_index {
        return 0.0;
    }
    match b.order() {
        BondOrder::Null | BondOrder::Unspecified | BondOrder::Zero => 0.0,
        BondOrder::Single => 1.0,
        BondOrder::Double => 2.0,
        BondOrder::Triple => 3.0,
        BondOrder::Quadruple => 4.0,
        BondOrder::Quintuple => 5.0,
        BondOrder::Hextuple => 6.0,
        BondOrder::OneAndHalf | BondOrder::Aromatic => 1.5,
        BondOrder::TwoAndHalf => 2.5,
        BondOrder::ThreeAndHalf => 3.5,
        BondOrder::FourAndHalf => 4.5,
        BondOrder::FiveAndHalf => 5.5,
        BondOrder::Ionic | BondOrder::Hydrogen | BondOrder::ThreeCenter | BondOrder::Other => 0.0,
        BondOrder::DativeOne | BondOrder::DativeLeft | BondOrder::DativeRight => 1.0,
        BondOrder::Dative => {
            if b.end().index() == atom_index {
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
    let atom = &mol.atoms()[atom_index];
    let Some(dv) = rdkit_default_valence(atom.atomic_number()) else {
        return -1;
    };
    if dv <= 1 {
        return -1;
    }
    let mut degree = atom_degree[atom_index] as i32
        + atom.explicit_hydrogens() as i32
        + assignment.implicit_hydrogens[atom_index] as i32;
    for b in mol.bonds() {
        if (b.begin().index() == atom_index || b.end().index() == atom_index)
            && bond_valence_contrib_for_atom(b, atom_index) == 0.0
        {
            degree -= 1;
        }
    }
    if degree > 3 {
        return -1;
    }
    let Some(nouter) = rdkit_n_outer_electrons(atom.atomic_number()) else {
        return -1;
    };
    let nlp = (nouter - dv - atom.formal_charge() as i32).max(0);
    let radicals = atom.radical_electrons() as i32;
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
    let at = &mol.atoms()[atom_index];
    if let Ok(Some(vals)) = rdkit_valence_list(at.atomic_number())
        && at.formal_charge() == 0
        && !vals.is_empty()
        && vals[0] >= 0
    {
        let total_valence = assignment.explicit_valence[atom_index] as i32
            + assignment.implicit_hydrogens[atom_index] as i32;
        if total_valence > vals[0] {
            return false;
        }
    }
    let nouter = rdkit_n_outer_electrons(at.atomic_number()).unwrap_or(0);
    let total_degree = atom_degree[atom_index]
        + at.explicit_hydrogens() as usize
        + assignment.implicit_hydrogens[atom_index] as usize;
    let row_ok = at.atomic_number() <= 10
        || (nouter != 5 && nouter != 6)
        || (nouter == 6 && total_degree < 2);
    row_ok && count_atom_electrons_rdkit(mol, assignment, atom_degree, atom_index) > 0
}

fn compute_conjugated_bonds(
    mol: &Molecule,
    assignment: &cosmolkit_core::ValenceAssignment,
    atom_degree: &[usize],
) -> Vec<bool> {
    let mut conjugated = vec![false; mol.bonds().len()];
    for (bi, b) in mol.bonds().iter().enumerate() {
        conjugated[bi] = b.is_aromatic();
    }
    for at in 0..mol.atoms().len() {
        if !is_atom_conjug_cand(mol, assignment, atom_degree, at) {
            continue;
        }
        let sbo = atom_degree[at]
            + mol.atoms()[at].explicit_hydrogens() as usize
            + assignment.implicit_hydrogens[at] as usize;
        if !(2..=3).contains(&sbo) {
            continue;
        }
        let bnds: Vec<usize> = mol
            .bonds()
            .iter()
            .enumerate()
            .filter_map(|(bi, b)| {
                if b.begin().index() == at || b.end().index() == at {
                    Some(bi)
                } else {
                    None
                }
            })
            .collect();
        for &b1 in &bnds {
            let bond1 = &mol.bonds()[b1];
            if bond_valence_contrib_for_atom(bond1, at) < 1.5 {
                continue;
            }
            let o1 = if bond1.begin().index() == at {
                bond1.end().index()
            } else {
                bond1.begin().index()
            };
            if !is_atom_conjug_cand(mol, assignment, atom_degree, o1) {
                continue;
            }
            for &b2 in &bnds {
                if b1 == b2 {
                    continue;
                }
                let bond2 = &mol.bonds()[b2];
                let o2 = if bond2.begin().index() == at {
                    bond2.end().index()
                } else {
                    bond2.begin().index()
                };
                let sbo2 = atom_degree[o2]
                    + mol.atoms()[o2].explicit_hydrogens() as usize
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
    let atom = &mol.atoms()[atom_index];
    if atom.atomic_number() == 1 {
        // RDKit parity corpus behavior: hydrogens introduced by AddHs() stay
        // UNSPECIFIED, while bracket/input hydrogens (including isotopes) are S.
        if !atom.no_implicit() && atom.isotope().is_none() {
            return Hybridization::Unspecified;
        }
    }
    if atom.atomic_number() == 0 {
        return Hybridization::Unspecified;
    }
    // RDKit shortcut for tetrahedral stereocenters.
    if !matches!(atom.chiral_tag(), OursChiralTag::Unspecified)
        && atom_degree[atom_index]
            + atom.explicit_hydrogens() as usize
            + assignment.implicit_hydrogens[atom_index] as usize
            == 4
    {
        return Hybridization::Sp3;
    }
    let mut deg = atom_degree[atom_index] as i32
        + atom.explicit_hydrogens() as i32
        + assignment.implicit_hydrogens[atom_index] as i32;
    for b in mol.bonds() {
        if (b.begin().index() == atom_index || b.end().index() == atom_index)
            && (matches!(b.order(), BondOrder::Dative) && b.end().index() != atom_index)
        {
            deg -= 1;
        }
    }
    if atom.atomic_number() <= 1 {
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
    let nouter = rdkit_n_outer_electrons(atom.atomic_number()).unwrap_or(0);
    let total_valence = assignment.explicit_valence[atom_index] as i32
        + assignment.implicit_hydrogens[atom_index] as i32;
    let num_free = nouter - (total_valence + atom.formal_charge() as i32);
    let norbs = if total_valence + nouter - (atom.formal_charge() as i32) < 8 {
        let radicals = atom.radical_electrons() as i32;
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
                + atom.explicit_hydrogens() as usize
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
    let mut atom_degree = vec![0usize; mol.atoms().len()];
    let mut atom_has_multi = vec![false; mol.atoms().len()];

    for b in mol.bonds() {
        atom_degree[b.begin().index()] += 1;
        atom_degree[b.end().index()] += 1;
        if matches!(
            b.order(),
            BondOrder::Double | BondOrder::Triple | BondOrder::Quadruple | BondOrder::Aromatic
        ) {
            atom_has_multi[b.begin().index()] = true;
            atom_has_multi[b.end().index()] = true;
        }
    }

    let conjugated_bonds = compute_conjugated_bonds(mol, &assignment, &atom_degree);
    let mut atom_has_conjugated_bond = vec![false; mol.atoms().len()];
    for (bi, b) in mol.bonds().iter().enumerate() {
        if conjugated_bonds[bi] {
            atom_has_conjugated_bond[b.begin().index()] = true;
            atom_has_conjugated_bond[b.end().index()] = true;
        }
    }

    let mut atoms = Vec::with_capacity(mol.atoms().len());
    for (i, a) in mol.atoms().iter().enumerate() {
        let implicit_hs = assignment.implicit_hydrogens[i] as i32;
        let explicit_valence = assignment.explicit_valence[i] as i32;
        let num_hs = (a.explicit_hydrogens() as i32 + implicit_hs).max(0) as usize;
        let degree = atom_degree[i] + num_hs;
        let total_valence = explicit_valence + implicit_hs;
        let _ = atom_has_multi[i];
        let hybridization =
            compute_hybridization(mol, &assignment, &atom_degree, &atom_has_conjugated_bond, i);
        atoms.push(OursAtomFeature {
            atomic_num: a.atomic_number(),
            isotope: a.isotope(),
            chirality: ours_chiral_tag_name(a.chiral_tag()),
            degree,
            formal_charge: a.formal_charge(),
            num_hs,
            num_radical_electrons: a.radical_electrons() as usize,
            hybridization,
            is_aromatic: a.is_aromatic(),
            is_in_ring: ring_flags[i],
            explicit_valence,
            implicit_hs,
            total_valence,
        });
    }

    let mut bonds = Vec::with_capacity(mol.bonds().len());
    for (bi, b) in mol.bonds().iter().enumerate() {
        bonds.push(OursBondFeature {
            begin_atom: b.begin().index(),
            end_atom: b.end().index(),
            bond_type: bond_type_name(b.order()),
            stereo: ours_bond_stereo_name(b.stereo()),
            is_conjugated: conjugated_bonds[bi],
        });
    }
    (atoms, bonds)
}

fn extract_possible_stereo_features(
    mol: &Molecule,
) -> (Vec<OursPossibleStereoAtomFeature>, Vec<OursBondFeature>) {
    let (atoms, bonds) = extract_ours_features(mol);
    let assignment = assign_valence(mol, ValenceModel::RdkitLike).unwrap_or_else(|e| {
        panic!(
            "assign_valence failed in possible-stereo test extraction: {:?}",
            e
        );
    });
    let keep_cip_ranks = legacy_possible_stereo_presence_requires_rank_work(mol);
    let cip_ranks = if keep_cip_ranks {
        legacy_rdkit_cip_ranks_for_depict_with_assignment(mol, &assignment)
    } else {
        Vec::new()
    };
    let chirality_possible_flags: Vec<bool> = if cip_ranks.is_empty() {
        vec![false; mol.atoms().len()]
    } else {
        (0..mol.atoms().len())
            .map(|atom_idx| {
                let (legal_center, has_dupes, _nbrs) =
                    legacy_is_atom_potential_chiral_center(mol, &assignment, atom_idx, &cip_ranks);
                legal_center && !has_dupes
            })
            .collect()
    };
    let cip_codes: Vec<Option<String>> = if keep_cip_ranks {
        (0..mol.atoms().len())
            .map(|atom_idx| legacy_cip_code_for_atom(mol, &assignment, &cip_ranks, atom_idx))
            .collect()
    } else {
        vec![None; mol.atoms().len()]
    };
    let cip_ranks = if keep_cip_ranks
        && cip_codes.iter().any(Option::is_some)
        && legacy_assign_bond_stereo_would_leave_unassigned(mol)
    {
        legacy_rdkit_cip_reranks_with_legacy_stereo_and_assignment(
            mol,
            &cip_ranks,
            &cip_codes,
            &assignment,
        )
    } else {
        cip_ranks
    };
    let atoms = atoms
        .into_iter()
        .enumerate()
        .map(|(atom_idx, base)| {
            let atom = &mol.atoms()[atom_idx];
            let cip_rank = if keep_cip_ranks {
                cip_ranks.get(atom_idx).copied()
            } else {
                None
            };
            let chirality_possible = chirality_possible_flags[atom_idx];
            OursPossibleStereoAtomFeature {
                base,
                cip_code: ours_cip_code_name(atom.prop("_CIPCode")),
                cip_rank,
                chirality_possible,
            }
        })
        .collect();
    (atoms, bonds)
}

fn atom_has_directional_bond_to_other_than(
    mol: &Molecule,
    atom_index: usize,
    excluded_neighbor: usize,
) -> bool {
    mol.bonds().iter().any(|bond| {
        let other = if bond.begin().index() == atom_index {
            bond.end().index()
        } else if bond.end().index() == atom_index {
            bond.begin().index()
        } else {
            return false;
        };
        other != excluded_neighbor
            && matches!(
                bond.direction(),
                BondDirection::EndDownRight | BondDirection::EndUpRight
            )
    })
}

fn min_cycle_size_for_bond(mol: &Molecule, bond_index: usize) -> Option<usize> {
    let bond = &mol.bonds()[bond_index];
    let start = bond.begin().index();
    let goal = bond.end().index();
    let mut q = VecDeque::new();
    let mut seen = vec![false; mol.atoms().len()];
    seen[start] = true;
    q.push_back((start, 0usize));
    while let Some((atom_idx, dist)) = q.pop_front() {
        for (idx, other) in mol.bonds().iter().enumerate().filter_map(|(idx, b)| {
            if idx == bond_index {
                return None;
            }
            if b.begin().index() == atom_idx {
                Some((idx, b.end().index()))
            } else if b.end().index() == atom_idx {
                Some((idx, b.begin().index()))
            } else {
                None
            }
        }) {
            let _ = idx;
            if other == goal {
                return Some(dist + 2);
            }
            if !seen[other] {
                seen[other] = true;
                q.push_back((other, dist + 1));
            }
        }
    }
    None
}

fn atom_nonzero_degree_for_legacy_stereo(mol: &Molecule, atom_index: usize) -> usize {
    mol.bonds()
        .iter()
        .filter(|bond| {
            (bond.begin().index() == atom_index || bond.end().index() == atom_index)
                && !(matches!(bond.order(), BondOrder::Null)
                    || matches!(bond.order(), BondOrder::Dative)
                        && bond.begin().index() == atom_index)
        })
        .count()
}

fn has_protium_neighbor_for_legacy_stereo(mol: &Molecule, atom_index: usize) -> bool {
    mol.bonds().iter().any(|bond| {
        let other = if bond.begin().index() == atom_index {
            Some(bond.end().index())
        } else if bond.end().index() == atom_index {
            Some(bond.begin().index())
        } else {
            None
        };
        other
            .map(|idx| {
                let atom = &mol.atoms()[idx];
                atom.atomic_number() == 1 && atom.isotope().is_none()
            })
            .unwrap_or(false)
    })
}

fn bond_is_conjugated_for_legacy_stereo(mol: &Molecule, bond_index: usize) -> bool {
    let bond = &mol.bonds()[bond_index];
    if bond.is_aromatic() || matches!(bond.order(), BondOrder::Double | BondOrder::Triple) {
        return true;
    }
    if !matches!(bond.order(), BondOrder::Single) {
        return false;
    }
    [bond.begin().index(), bond.end().index()]
        .into_iter()
        .any(|center| {
            mol.bonds().iter().enumerate().any(|(other_idx, other)| {
                other_idx != bond_index
                    && (other.begin().index() == center || other.end().index() == center)
                    && (other.is_aromatic()
                        || matches!(other.order(), BondOrder::Double | BondOrder::Triple))
            })
        })
}

fn atom_has_conjugated_bond_for_legacy_stereo(mol: &Molecule, atom_index: usize) -> bool {
    mol.bonds().iter().enumerate().any(|(bond_index, bond)| {
        (bond.begin().index() == atom_index || bond.end().index() == atom_index)
            && bond_is_conjugated_for_legacy_stereo(mol, bond_index)
    })
}

fn is_atom_bridgehead_for_legacy_stereo(mol: &Molecule, atom_index: usize) -> bool {
    if atom_nonzero_degree_for_legacy_stereo(mol, atom_index) < 3 {
        return false;
    }
    let atom_ring_bonds = mol
        .bonds()
        .iter()
        .enumerate()
        .filter(|(bond_index, bond)| {
            (bond.begin().index() == atom_index || bond.end().index() == atom_index)
                && min_cycle_size_for_bond(mol, *bond_index).is_some()
        })
        .count();
    atom_ring_bonds >= 3
}

fn legacy_is_atom_potential_chiral_center(
    mol: &Molecule,
    assignment: &cosmolkit_core::ValenceAssignment,
    atom_index: usize,
    ranks: &[i64],
) -> (bool, bool, Vec<(u32, usize)>) {
    let atom = &mol.atoms()[atom_index];
    let nz_degree = atom_nonzero_degree_for_legacy_stereo(mol, atom_index);
    let total_num_hs = atom.explicit_hydrogens() as usize
        + assignment.implicit_hydrogens[atom_index].max(0) as usize
        + mol
            .bonds()
            .iter()
            .filter(|bond| {
                let other = if bond.begin().index() == atom_index {
                    Some(bond.end().index())
                } else if bond.end().index() == atom_index {
                    Some(bond.begin().index())
                } else {
                    None
                };
                other
                    .map(|idx| mol.atoms()[idx].atomic_number() == 1)
                    .unwrap_or(false)
            })
            .count();
    let tnz_degree = nz_degree + total_num_hs;
    let mut legal_center = true;
    let mut has_dupes = false;
    let mut nbrs = Vec::new();

    if tnz_degree > 4
        || tnz_degree < 3
        || (nz_degree < 3 && !matches!(atom.atomic_number(), 15 | 33))
    {
        legal_center = false;
    } else if nz_degree == 3 {
        if total_num_hs == 1 {
            if has_protium_neighbor_for_legacy_stereo(mol, atom_index) {
                legal_center = false;
            }
        } else {
            legal_center = false;
            match atom.atomic_number() {
                7 => {
                    let in_three_ring = mol.bonds().iter().enumerate().any(|(bond_index, bond)| {
                        (bond.begin().index() == atom_index || bond.end().index() == atom_index)
                            && min_cycle_size_for_bond(mol, bond_index) == Some(3)
                    });
                    if !atom_has_conjugated_bond_for_legacy_stereo(mol, atom_index)
                        && (in_three_ring || is_atom_bridgehead_for_legacy_stereo(mol, atom_index))
                    {
                        legal_center = true;
                    }
                }
                15 | 33 => {
                    legal_center = true;
                }
                16 | 34 => {
                    let ev = assignment.explicit_valence[atom_index];
                    if ev == 4 || (ev == 3 && atom.formal_charge() == 1) {
                        legal_center = true;
                    }
                }
                _ => {}
            }
        }
    }

    if legal_center && !ranks.is_empty() {
        let mut seen_ranks = std::collections::HashMap::<i64, ()>::new();
        for bond in mol.bonds() {
            let other = if bond.begin().index() == atom_index {
                Some(bond.end().index())
            } else if bond.end().index() == atom_index {
                Some(bond.begin().index())
            } else {
                None
            };
            let Some(other_idx) = other else {
                continue;
            };
            nbrs.push((ranks[other_idx] as u32, bond.id().index()));
            if !(matches!(bond.order(), BondOrder::Null)
                || matches!(bond.order(), BondOrder::Dative) && bond.begin().index() == atom_index)
                && seen_ranks.insert(ranks[other_idx], ()).is_some()
            {
                has_dupes = true;
                break;
            }
        }
    }

    (legal_center, has_dupes, nbrs)
}

fn legacy_rdkit_twice_bond_type(order: BondOrder) -> i64 {
    match order {
        BondOrder::Null => 0,
        BondOrder::Single => 2,
        BondOrder::Double => 4,
        BondOrder::Triple => 6,
        BondOrder::Quadruple => 8,
        BondOrder::Aromatic => 3,
        BondOrder::Dative => 2,
        BondOrder::Hydrogen => 0,
        BondOrder::Quintuple
        | BondOrder::Hextuple
        | BondOrder::OneAndHalf
        | BondOrder::TwoAndHalf
        | BondOrder::ThreeAndHalf
        | BondOrder::FourAndHalf
        | BondOrder::FiveAndHalf
        | BondOrder::Ionic
        | BondOrder::DativeOne
        | BondOrder::DativeLeft
        | BondOrder::DativeRight
        | BondOrder::ThreeCenter
        | BondOrder::Other
        | BondOrder::Zero
        | BondOrder::Unspecified => {
            panic!(
                "legacy RDKit depict CIP ranks do not support bond order {:?}",
                order
            )
        }
    }
}

fn legacy_most_common_isotope(atomic_num: u8) -> Option<i64> {
    match atomic_num {
        1 => Some(1),
        5 => Some(11),
        6 => Some(12),
        7 => Some(14),
        8 => Some(16),
        9 => Some(19),
        11 => Some(23),
        12 => Some(24),
        13 => Some(27),
        14 => Some(28),
        15 => Some(31),
        16 => Some(32),
        17 => Some(35),
        19 => Some(39),
        20 => Some(40),
        35 => Some(79),
        53 => Some(127),
        n if n <= 20 => Some((n as i64) * 2),
        _ => Some((atomic_num as i64).max(1) * 2),
    }
}

fn legacy_rdkit_cip_invariants(mol: &Molecule) -> Vec<i64> {
    const MASS_BITS: i64 = 10;
    const MAX_MASS: i64 = 1 << MASS_BITS;

    mol.atoms()
        .iter()
        .map(|atom| {
            let mut invariant = i64::from(atom.atomic_number() % 128);
            let mut mass = 0i64;
            if let Some(isotope) = atom.isotope() {
                mass = i64::from(isotope)
                    - legacy_most_common_isotope(atom.atomic_number()).unwrap_or(0);
                if mass >= 0 {
                    mass += 1;
                }
            }
            mass += MAX_MASS / 2;
            if mass < 0 {
                mass = 0;
            } else {
                mass %= MAX_MASS;
            }
            invariant = (invariant << MASS_BITS) | mass;
            let mapnum = atom
                .atom_map()
                .map(|m| ((m as i64) + 1) % 1024)
                .unwrap_or(0);
            (invariant << 10) | mapnum
        })
        .collect()
}

fn legacy_find_cip_segments(
    sorted: &mut [(usize, i64)],
    cip_entries: &[Vec<i64>],
) -> (Vec<(usize, usize)>, usize) {
    let mut segments = Vec::new();
    if sorted.is_empty() {
        return (segments, 0);
    }

    let mut num_independent = sorted.len();
    let mut running_rank = 0i64;
    sorted[0].1 = running_rank;
    let mut current = 0usize;
    let mut in_equal_section = false;

    for i in 1..sorted.len() {
        if cip_entries[sorted[current].0] == cip_entries[sorted[i].0] {
            sorted[i].1 = running_rank;
            num_independent -= 1;
            if !in_equal_section {
                in_equal_section = true;
                segments.push((i - 1, 0));
            }
        } else {
            running_rank += 1;
            sorted[i].1 = running_rank;
            current = i;

            if in_equal_section {
                if let Some((_, last)) = segments.last_mut() {
                    *last = i;
                }
                in_equal_section = false;
            }
        }
    }

    if in_equal_section {
        if let Some((_, last)) = segments.last_mut() {
            *last = sorted.len() - 1;
        }
    }

    (segments, num_independent)
}

fn legacy_recompute_cip_ranks(sorted: &[(usize, i64)], ranks: &mut [i64]) {
    for &(atom_idx, rank) in sorted {
        ranks[atom_idx] = rank;
    }
}

fn legacy_count_swaps_to_interconvert(probe: &[usize], reference: &[usize]) -> usize {
    let mut probe = probe.to_vec();
    let mut positions: std::collections::HashMap<usize, usize> =
        probe.iter().enumerate().map(|(i, &v)| (v, i)).collect();
    let mut swaps = 0usize;
    for (i, &target) in reference.iter().enumerate() {
        if probe[i] == target {
            continue;
        }
        let j = *positions
            .get(&target)
            .expect("reference element missing from probe permutation");
        let current = probe[i];
        probe.swap(i, j);
        positions.insert(current, j);
        positions.insert(target, i);
        swaps += 1;
    }
    swaps
}

fn legacy_perturbation_order_for_atom(mol: &Molecule, atom_index: usize, probe: &[usize]) -> usize {
    let reference = mol
        .bonds()
        .iter()
        .filter_map(|bond| {
            if bond.begin().index() == atom_index || bond.end().index() == atom_index {
                Some(bond.id().index())
            } else {
                None
            }
        })
        .collect::<Vec<_>>();
    legacy_count_swaps_to_interconvert(probe, &reference)
}

fn legacy_cip_code_for_atom(
    mol: &Molecule,
    assignment: &cosmolkit_core::ValenceAssignment,
    ranks: &[i64],
    atom_index: usize,
) -> Option<String> {
    let atom = &mol.atoms()[atom_index];
    if matches!(atom.chiral_tag(), OursChiralTag::Unspecified) {
        return None;
    }

    let (legal_center, has_dupes, mut nbrs) =
        legacy_is_atom_potential_chiral_center(mol, assignment, atom_index, ranks);
    if !legal_center || has_dupes {
        return None;
    }

    nbrs.sort();
    let probe: Vec<usize> = nbrs.into_iter().map(|(_, bond_index)| bond_index).collect();
    let mut swaps = legacy_perturbation_order_for_atom(mol, atom_index, &probe);
    let total_num_hs = atom.explicit_hydrogens() as usize
        + assignment.implicit_hydrogens[atom_index] as usize
        + mol
            .bonds()
            .iter()
            .filter(|bond| {
                let other = if bond.begin().index() == atom_index {
                    Some(bond.end().index())
                } else if bond.end().index() == atom_index {
                    Some(bond.begin().index())
                } else {
                    None
                };
                other
                    .map(|idx| mol.atoms()[idx].atomic_number() == 1)
                    .unwrap_or(false)
            })
            .count();
    if probe.len() == 3 && total_num_hs == 1 {
        swaps += 1;
    }

    let mut final_tag = atom.chiral_tag();
    if swaps % 2 == 1 {
        final_tag = match final_tag {
            OursChiralTag::TetrahedralCcw => OursChiralTag::TetrahedralCw,
            OursChiralTag::TetrahedralCw => OursChiralTag::TetrahedralCcw,
            OursChiralTag::TrigonalBipyramidal => OursChiralTag::TrigonalBipyramidal,
            OursChiralTag::Unspecified => OursChiralTag::Unspecified,
            _ => final_tag,
        };
    }

    Some(
        if matches!(final_tag, OursChiralTag::TetrahedralCcw) {
            "S"
        } else {
            "R"
        }
        .to_owned(),
    )
}

fn legacy_rdkit_cip_ranks_from_invariants_with_assignment(
    mol: &Molecule,
    invars: Vec<i64>,
    seed_with_invars: bool,
    assignment: &cosmolkit_core::ValenceAssignment,
) -> Vec<i64> {
    let n = mol.atoms().len();
    let coordinated_h_counts = (0..n)
        .map(|index| {
            mol.atoms()[index].explicit_hydrogens() as usize
                + assignment.implicit_hydrogens[index] as usize
        })
        .collect::<Vec<_>>();
    let mut cip_entries: Vec<Vec<i64>> = invars
        .iter()
        .copied()
        .map(|v| {
            let mut entry = Vec::with_capacity(16);
            entry.push(v);
            entry
        })
        .collect();
    let mut sortable: Vec<(usize, i64)> = (0..n).map(|idx| (idx, -1)).collect();

    sortable.sort_by(|a, b| cip_entries[a.0].cmp(&cip_entries[b.0]));
    let (mut needs_sorting, mut num_ranks) = legacy_find_cip_segments(&mut sortable, &cip_entries);
    let mut ranks = vec![0i64; n];
    legacy_recompute_cip_ranks(&sortable, &mut ranks);

    for i in 0..n {
        if seed_with_invars {
            cip_entries[i][0] = invars[i];
        } else {
            cip_entries[i][0] = i64::from(mol.atoms()[i].atomic_number());
            cip_entries[i].push(ranks[i]);
        }
    }
    let cip_rank_index = if seed_with_invars { 1 } else { 2 };

    const K_MAX_BONDS: usize = 16;
    let mut bond_features = vec![(0i64, 0usize); n * K_MAX_BONDS];
    let mut num_neighbors = vec![0usize; n];
    let mut atom_degrees = vec![0usize; n];
    for bond in mol.bonds() {
        atom_degrees[bond.begin().index()] += 1;
        atom_degrees[bond.end().index()] += 1;
    }
    for bond in mol.bonds() {
        let begin_count = if matches!(bond.order(), BondOrder::Double) {
            let nbr = &mol.atoms()[bond.end().index()];
            if nbr.atomic_number() == 15 && matches!(atom_degrees[bond.end().index()], 3 | 4) {
                1
            } else {
                legacy_rdkit_twice_bond_type(bond.order())
            }
        } else {
            legacy_rdkit_twice_bond_type(bond.order())
        };
        let end_count = if matches!(bond.order(), BondOrder::Double) {
            let nbr = &mol.atoms()[bond.begin().index()];
            if nbr.atomic_number() == 15 && matches!(atom_degrees[bond.begin().index()], 3 | 4) {
                1
            } else {
                legacy_rdkit_twice_bond_type(bond.order())
            }
        } else {
            legacy_rdkit_twice_bond_type(bond.order())
        };

        let begin = bond.begin().index();
        let end = bond.end().index();

        let begin_offset = begin * K_MAX_BONDS + num_neighbors[begin];
        if begin_offset < (begin + 1) * K_MAX_BONDS {
            bond_features[begin_offset] = (begin_count, end);
        }
        num_neighbors[begin] += 1;

        let end_offset = end * K_MAX_BONDS + num_neighbors[end];
        if end_offset < (end + 1) * K_MAX_BONDS {
            bond_features[end_offset] = (end_count, begin);
        }
        num_neighbors[end] += 1;
    }

    let max_its = n / 2 + 1;
    let mut num_its = 0usize;
    let mut last_num_ranks: Option<usize> = None;

    while !needs_sorting.is_empty()
        && num_its < max_its
        && last_num_ranks.is_none_or(|last| last < num_ranks)
    {
        for index in 0..n {
            let index_offset = K_MAX_BONDS * index;
            let feature_len = num_neighbors[index].min(K_MAX_BONDS);
            if num_neighbors[index] > 1 {
                bond_features[index_offset..index_offset + feature_len]
                    .sort_by(|a, b| ranks[b.1].cmp(&ranks[a.1]));
            }
            for &(count, nbr_idx) in &bond_features[index_offset..index_offset + feature_len] {
                for _ in 0..count {
                    cip_entries[index].push(ranks[nbr_idx] + 1);
                }
            }
            let new_len = cip_entries[index].len() + coordinated_h_counts[index];
            cip_entries[index].resize(new_len, 0);
        }
        last_num_ranks = Some(num_ranks);
        for &(first, last) in &needs_sorting {
            sortable[first..=last].sort_by(|a, b| cip_entries[a.0].cmp(&cip_entries[b.0]));
        }
        let found = legacy_find_cip_segments(&mut sortable, &cip_entries);
        needs_sorting = found.0;
        num_ranks = found.1;
        legacy_recompute_cip_ranks(&sortable, &mut ranks);

        if Some(num_ranks) != last_num_ranks {
            for i in 0..n {
                cip_entries[i].resize(cip_rank_index + 1, 0);
                cip_entries[i][cip_rank_index] = ranks[i];
            }
        }

        num_its += 1;
    }

    ranks
}

fn legacy_rdkit_cip_ranks_for_depict_with_assignment(
    mol: &Molecule,
    assignment: &cosmolkit_core::ValenceAssignment,
) -> Vec<i64> {
    legacy_rdkit_cip_ranks_from_invariants_with_assignment(
        mol,
        legacy_rdkit_cip_invariants(mol),
        false,
        assignment,
    )
}

fn legacy_rdkit_cip_reranks_with_legacy_stereo_and_assignment(
    mol: &Molecule,
    ranks: &[i64],
    cip_codes: &[Option<String>],
    assignment: &cosmolkit_core::ValenceAssignment,
) -> Vec<i64> {
    let mut factor = 100i64;
    while factor < mol.atoms().len() as i64 {
        factor *= 10;
    }
    let mut double_bond_stereo_contrib = vec![0i64; mol.atoms().len()];
    for bond in mol.bonds() {
        if !matches!(bond.order(), BondOrder::Double) {
            continue;
        }
        let contrib = match bond.stereo() {
            BondStereo::Trans => 1,
            BondStereo::Cis => 2,
            BondStereo::None | BondStereo::Any => 0,
            BondStereo::Z | BondStereo::E | BondStereo::AtropCw | BondStereo::AtropCcw => {
                panic!(
                    "legacy RDKit stereo rerank does not support bond stereo {:?}",
                    bond.stereo()
                )
            }
        };
        double_bond_stereo_contrib[bond.begin().index()] += contrib;
        double_bond_stereo_contrib[bond.end().index()] += contrib;
    }
    let mut invars = vec![0i64; mol.atoms().len()];
    for atom in mol.atoms() {
        let mut invariant = ranks[atom.id().index()] * factor;
        if let Some(code) = &cip_codes[atom.id().index()] {
            if code == "S" {
                invariant += 10;
            } else if code == "R" {
                invariant += 20;
            }
        }
        invariant += double_bond_stereo_contrib[atom.id().index()];
        invars[atom.id().index()] = invariant;
    }
    legacy_rdkit_cip_ranks_from_invariants_with_assignment(mol, invars, true, assignment)
}

fn legacy_assign_bond_stereo_would_leave_unassigned(mol: &Molecule) -> bool {
    for (bond_index, bond) in mol.bonds().iter().enumerate() {
        if !matches!(bond.order(), BondOrder::Double)
            || min_cycle_size_for_bond(mol, bond_index).is_some_and(|size| size < 8)
        {
            continue;
        }
        let begin = bond.begin().index();
        let end = bond.end().index();
        let begin_degree = atom_nonzero_degree_for_legacy_stereo(mol, begin);
        let end_degree = atom_nonzero_degree_for_legacy_stereo(mol, end);
        if !matches!(begin_degree, 2 | 3) || !matches!(end_degree, 2 | 3) {
            continue;
        }
        let begin_has_dir = atom_has_directional_bond_to_other_than(mol, begin, end);
        let end_has_dir = atom_has_directional_bond_to_other_than(mol, end, begin);
        if !(begin_has_dir && end_has_dir) {
            return true;
        }
    }
    false
}

fn legacy_possible_stereo_presence_requires_rank_work(mol: &Molecule) -> bool {
    let assignment = assign_valence(mol, ValenceModel::RdkitLike).unwrap_or_else(|e| {
        panic!(
            "assign_valence failed in possible-stereo presence scan: {:?}",
            e
        );
    });
    let has_stereo_atoms = mol
        .atoms()
        .iter()
        .any(|atom| !matches!(atom.chiral_tag(), OursChiralTag::Unspecified));
    let has_potential_stereo_atoms = mol.atoms().iter().enumerate().any(|(atom_idx, _)| {
        legacy_is_atom_potential_chiral_center(mol, &assignment, atom_idx, &[]).0
    });

    let mut has_stereo_bonds = false;
    let mut has_potential_stereo_bonds = false;
    for bond in mol.bonds() {
        if !matches!(bond.order(), BondOrder::Double) {
            continue;
        }
        let begin = bond.begin().index();
        let end = bond.end().index();
        let is_specified = atom_has_directional_bond_to_other_than(mol, begin, end)
            || atom_has_directional_bond_to_other_than(mol, end, begin);
        if is_specified {
            has_stereo_bonds = true;
        } else if !has_potential_stereo_bonds
            && cosmolkit_core::stereo::should_detect_double_bond_stereo(mol, bond.id())
                .unwrap_or(false)
        {
            has_potential_stereo_bonds = true;
        }
        if has_stereo_bonds && has_potential_stereo_bonds {
            break;
        }
    }

    has_stereo_atoms || has_stereo_bonds || has_potential_stereo_atoms || has_potential_stereo_bonds
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
            a.isotope, e.isotope,
            "isotope mismatch row {} atom {} ({}) [{}]",
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

fn compare_possible_stereo_features(
    ours_atoms: &[OursPossibleStereoAtomFeature],
    ours_bonds: &[OursBondFeature],
    expected: &FeatureSet,
    row_idx: usize,
    smiles: &str,
) {
    let base_atoms: Vec<OursAtomFeature> = ours_atoms.iter().map(|a| a.base.clone()).collect();
    compare_features(
        &base_atoms,
        ours_bonds,
        expected,
        row_idx,
        smiles,
        "possible_stereo",
    );

    for (i, (a, e)) in ours_atoms
        .iter()
        .zip(expected.atom_features.iter())
        .enumerate()
    {
        assert_eq!(
            a.cip_code, e.cip_code,
            "cip_code mismatch row {} atom {} ({}) [possible_stereo]",
            row_idx, i, smiles
        );
        assert_eq!(
            a.cip_rank, e.cip_rank,
            "cip_rank mismatch row {} atom {} ({}) [possible_stereo]",
            row_idx, i, smiles
        );
        assert_eq!(
            a.chirality_possible, e.chirality_possible,
            "chirality_possible mismatch row {} atom {} ({}) [possible_stereo]",
            row_idx, i, smiles
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
                record.direct.is_some()
                    && record.with_hs.is_some()
                    && record.addhs_removehs.is_some()
                    && record.possible_stereo.is_some()
                    && record.chiral_centers.is_some(),
                "rdkit_ok=true requires direct, with_hs, addhs_removehs, possible_stereo, and chiral_centers at row {}",
                idx + 1
            );
            assert!(
                record.error.is_none(),
                "rdkit_ok=true should not carry error at row {}",
                idx + 1
            );
        } else {
            assert!(
                record.direct.is_none()
                    && record.with_hs.is_none()
                    && record.addhs_removehs.is_none()
                    && record.possible_stereo.is_none()
                    && record.chiral_centers.is_none(),
                "rdkit_ok=false should not carry feature/chiral-center sets at row {}",
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
fn graph_feature_golden_records_chiral_centers_with_include_unassigned_branches() {
    let golden = load_golden().expect("should read tests/golden/graph_features.jsonl");
    for (idx, record) in golden.iter().enumerate() {
        if !record.rdkit_ok {
            continue;
        }
        let centers = record
            .chiral_centers
            .as_ref()
            .expect("chiral_centers missing for rdkit_ok row");

        for c in &centers.include_unassigned_false {
            let found = centers
                .include_unassigned_true
                .iter()
                .any(|x| x.atom_idx == c.atom_idx && x.label == c.label);
            assert!(
                found,
                "includeUnassigned=true should contain includeUnassigned=false center at row {} atom {} ({})",
                idx + 1,
                c.atom_idx,
                record.smiles
            );
        }
    }
}

#[test]
fn graph_feature_golden_records_flag_possible_stereo_centers_branch() {
    let golden = load_golden().expect("should read tests/golden/graph_features.jsonl");
    let mut possible_center_rows = 0usize;
    let required_smiles = ["CC(F)(Cl)Br", "FC(Cl)(Br)I", "CC(F)Cl", "CC(Cl)Br"];

    for required in required_smiles {
        assert!(
            golden.iter().any(|record| record.smiles == required),
            "SMILES corpus should contain explicit potential-stereocenter fixture {required}"
        );
    }

    for (idx, record) in golden.iter().enumerate() {
        if !record.rdkit_ok {
            continue;
        }
        let direct = record.direct.as_ref().expect("direct missing");
        let possible = record
            .possible_stereo
            .as_ref()
            .expect("possible_stereo missing for rdkit_ok row");
        assert_eq!(
            direct.atom_features.len(),
            possible.atom_features.len(),
            "possible_stereo atom count mismatch at row {} ({})",
            idx + 1,
            record.smiles
        );
        assert_eq!(
            direct.bond_features.len(),
            possible.bond_features.len(),
            "possible_stereo bond count mismatch at row {} ({})",
            idx + 1,
            record.smiles
        );

        if possible
            .atom_features
            .iter()
            .any(|atom| atom.chirality_possible)
        {
            possible_center_rows += 1;
        }
    }

    assert!(
        possible_center_rows >= 4,
        "SMILES corpus should exercise RDKit AssignStereochemistry(flagPossibleStereoCenters=True)"
    );
}

#[test]
fn graph_feature_golden_records_isotopes_from_smiles() {
    let golden = load_golden().expect("should read tests/golden/graph_features.jsonl");
    let mut isotopes = Vec::new();

    for record in golden.iter().filter(|record| record.rdkit_ok) {
        let direct = record.direct.as_ref().expect("direct missing");
        isotopes.extend(direct.atom_features.iter().filter_map(|atom| atom.isotope));
    }

    assert!(
        isotopes.contains(&13),
        "SMILES corpus should include a carbon-13 isotope"
    );
    assert!(
        isotopes.contains(&2),
        "SMILES corpus should include deuterium"
    );
}

#[test]
fn graph_feature_golden_preserves_rdkit_cip_fields_for_chiral_atoms() {
    let golden = load_golden().expect("should read tests/golden/graph_features.jsonl");
    let mut chiral_atoms = 0usize;
    let mut cip_code_count = 0usize;
    let mut cip_rank_count = 0usize;

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
            if atom.cip_rank.is_some() {
                cip_rank_count += 1;
            }
            if let Some(code) = &atom.cip_code {
                assert!(
                    matches!(
                        code,
                        CipCode::R | CipCode::S | CipCode::LowerR | CipCode::LowerS
                    ),
                    "RDKit _CIPCode must be a known CIP label at row {} atom {} ({})",
                    idx + 1,
                    atom_idx,
                    record.smiles
                );
                cip_code_count += 1;
            }
            assert!(
                !atom.chirality_possible,
                "chiral atom should not carry RDKit _ChiralityPossible in this corpus at row {} atom {} ({})",
                idx + 1,
                atom_idx,
                record.smiles
            );
        }
    }

    assert!(
        chiral_atoms >= 8,
        "SMILES corpus should include multiple RDKit-recognized chiral atoms"
    );
    assert!(
        cip_code_count >= 8,
        "SMILES corpus should include multiple RDKit-assigned CIP codes"
    );
    assert!(
        cip_rank_count >= cip_code_count,
        "atoms with RDKit _CIPCode should also be represented in the ranked CIP field population"
    );
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
                let mol_with_h = mol_direct.with_hydrogens().unwrap_or_else(|e| {
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

                let batch_smiles = vec![record.smiles.clone(), record.smiles.clone()];
                let batch = MoleculeBatch::from_smiles_list(&batch_smiles);
                for (batch_idx, batch_record) in batch.iter().enumerate() {
                    let BatchRecord::Molecule(batch_mol) = batch_record else {
                        panic!(
                            "batch direct graph feature missing at row {} ({}) duplicate {}",
                            idx + 1,
                            record.smiles,
                            batch_idx
                        );
                    };
                    let (batch_atoms, batch_bonds) = extract_ours_features(batch_mol);
                    compare_features(
                        &batch_atoms,
                        &batch_bonds,
                        direct_expected,
                        idx + 1,
                        &record.smiles,
                        "batch_direct",
                    );
                }

                let batch_with_h = batch
                    .with_hydrogens(BatchErrorMode::KeepErrors)
                    .expect("batch add_hydrogens should succeed after scalar branch passed");
                for (batch_idx, batch_record) in batch_with_h.iter().enumerate() {
                    let BatchRecord::Molecule(batch_mol) = batch_record else {
                        panic!(
                            "batch explicit-hydrogen graph feature missing at row {} ({}) duplicate {}",
                            idx + 1,
                            record.smiles,
                            batch_idx
                        );
                    };
                    let (batch_atoms, batch_bonds) = extract_ours_features(batch_mol);
                    compare_features(
                        &batch_atoms,
                        &batch_bonds,
                        with_h_expected,
                        idx + 1,
                        &record.smiles,
                        "batch_with_hs",
                    );
                }
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

#[test]
fn graph_features_match_rdkit_golden_for_addhs_removehs_roundtrip_branch() {
    let golden = load_golden().expect("should read tests/golden/graph_features.jsonl");

    for (idx, record) in golden.iter().enumerate() {
        let ours = Molecule::from_smiles(&record.smiles);
        match (record.rdkit_ok, ours) {
            (true, Ok(mol)) => {
                let mol = mol.with_hydrogens().unwrap_or_else(|e| {
                    panic!(
                        "add_hydrogens failed at row {} ({}): {:?}",
                        idx + 1,
                        record.smiles,
                        e
                    )
                });
                let mol = mol.without_hydrogens().unwrap_or_else(|e| {
                    panic!(
                        "remove_hydrogens failed at row {} ({}): {:?}",
                        idx + 1,
                        record.smiles,
                        e
                    )
                });

                let expected = record
                    .addhs_removehs
                    .as_ref()
                    .expect("addhs_removehs missing");
                let (ours_atoms, ours_bonds) = extract_ours_features(&mol);
                compare_features(
                    &ours_atoms,
                    &ours_bonds,
                    expected,
                    idx + 1,
                    &record.smiles,
                    "addhs_removehs",
                );

                let batch_smiles = vec![record.smiles.clone(), record.smiles.clone()];
                let batch = MoleculeBatch::from_smiles_list(&batch_smiles)
                    .with_hydrogens(BatchErrorMode::KeepErrors)
                    .expect("batch add_hydrogens should succeed after scalar branch passed")
                    .without_hydrogens(BatchErrorMode::KeepErrors)
                    .expect("batch remove_hydrogens should succeed after scalar branch passed");
                for (batch_idx, batch_record) in batch.iter().enumerate() {
                    let BatchRecord::Molecule(batch_mol) = batch_record else {
                        panic!(
                            "batch add/remove hydrogens graph feature missing at row {} ({}) duplicate {}",
                            idx + 1,
                            record.smiles,
                            batch_idx
                        );
                    };
                    let (batch_atoms, batch_bonds) = extract_ours_features(batch_mol);
                    compare_features(
                        &batch_atoms,
                        &batch_bonds,
                        expected,
                        idx + 1,
                        &record.smiles,
                        "batch_addhs_removehs",
                    );
                }
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

#[test]
fn graph_features_match_rdkit_golden_for_flag_possible_stereo_centers_branch() {
    let golden = load_golden().expect("should read tests/golden/graph_features.jsonl");

    for (idx, record) in golden.iter().enumerate() {
        let ours = Molecule::from_smiles(&record.smiles);
        match (record.rdkit_ok, ours) {
            (true, Ok(mol)) => {
                let expected = record
                    .possible_stereo
                    .as_ref()
                    .expect("possible_stereo missing");
                let (ours_atoms, ours_bonds) = extract_possible_stereo_features(&mol);
                compare_possible_stereo_features(
                    &ours_atoms,
                    &ours_bonds,
                    expected,
                    idx + 1,
                    &record.smiles,
                );

                let batch_smiles = vec![record.smiles.clone(), record.smiles.clone()];
                let batch = MoleculeBatch::from_smiles_list(&batch_smiles);
                for (batch_idx, batch_record) in batch.iter().enumerate() {
                    let BatchRecord::Molecule(batch_mol) = batch_record else {
                        panic!(
                            "batch possible-stereo graph feature missing at row {} ({}) duplicate {}",
                            idx + 1,
                            record.smiles,
                            batch_idx
                        );
                    };
                    let (batch_atoms, batch_bonds) = extract_possible_stereo_features(batch_mol);
                    compare_possible_stereo_features(
                        &batch_atoms,
                        &batch_bonds,
                        expected,
                        idx + 1,
                        &format!("{} duplicate {batch_idx}", record.smiles),
                    );
                }
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
