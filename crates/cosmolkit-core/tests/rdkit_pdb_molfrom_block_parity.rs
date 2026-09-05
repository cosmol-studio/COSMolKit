use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use cosmolkit_core::{BondDirection, BondStereo, Molecule};
use serde::Deserialize;

#[derive(Debug, Deserialize)]
struct RdkitPdbOracle {
    pdb_id: String,
    source_url: String,
    atoms: Vec<RdkitAtom>,
    bonds: Vec<RdkitBond>,
    coords: Vec<[f64; 3]>,
}

#[derive(Debug, Deserialize)]
struct RdkitAtom {
    idx: usize,
    atomic_num: u8,
    formal_charge: i8,
    isotope: Option<u16>,
    is_aromatic: bool,
    chiral_tag: String,
    radical_electrons: u8,
    pdb: RdkitPdbResidueInfo,
}

#[derive(Debug, Deserialize)]
struct RdkitPdbResidueInfo {
    atom_name: String,
    serial_number: i32,
    alt_loc: String,
    residue_name: String,
    residue_number: i32,
    chain_id: String,
    insertion_code: String,
    occupancy: f64,
    temp_factor: f64,
    is_hetero_atom: bool,
}

#[derive(Debug, Deserialize)]
struct RdkitBond {
    idx: usize,
    begin: usize,
    end: usize,
    bond_type: String,
    is_aromatic: bool,
    direction: String,
    stereo: String,
    stereo_atoms: Vec<usize>,
}

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn python_path() -> PathBuf {
    repo_root().join(".venv/bin/python")
}

fn cached_pdb_path() -> PathBuf {
    repo_root().join("tmp/7SH6.pdb")
}

fn ensure_7sh6_pdb() -> PathBuf {
    let path = cached_pdb_path();
    if path.exists() {
        return path;
    }

    fs::create_dir_all(path.parent().expect("tmp path should have a parent"))
        .expect("should create tmp cache directory");

    let script = r#"
from pathlib import Path
import sys
import urllib.request

out = Path(sys.argv[1])
url = "https://files.rcsb.org/download/7SH6.pdb"
data = urllib.request.urlopen(url, timeout=60).read()
out.write_bytes(data)
"#;
    let output = Command::new(python_path())
        .arg("-c")
        .arg(script)
        .arg(&path)
        .current_dir(repo_root())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .output()
        .expect("should run project Python to download 7SH6 PDB");
    assert!(
        output.status.success(),
        "failed to download 7SH6 PDB into {}:\nstdout:\n{}\nstderr:\n{}",
        path.display(),
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    path
}

fn rdkit_oracle_for_pdb(path: &Path) -> RdkitPdbOracle {
    let script = r#"
from pathlib import Path
import json
import sys

from rdkit import Chem

pdb_path = Path(sys.argv[1])
text = pdb_path.read_text()
mol = Chem.MolFromPDBBlock(text)
if mol is None:
    raise SystemExit("RDKit MolFromPDBBlock returned None")
conf = mol.GetConformer()

atoms = []
coords = []
for atom in mol.GetAtoms():
    info = atom.GetPDBResidueInfo()
    if info is None:
        raise SystemExit(f"atom {atom.GetIdx()} has no PDB residue info")
    pos = conf.GetAtomPosition(atom.GetIdx())
    isotope = atom.GetIsotope()
    atoms.append({
        "idx": atom.GetIdx(),
        "atomic_num": atom.GetAtomicNum(),
        "formal_charge": atom.GetFormalCharge(),
        "isotope": isotope if isotope != 0 else None,
        "is_aromatic": atom.GetIsAromatic(),
        "chiral_tag": str(atom.GetChiralTag()),
        "radical_electrons": atom.GetNumRadicalElectrons(),
        "pdb": {
            "atom_name": info.GetName(),
            "serial_number": info.GetSerialNumber(),
            "alt_loc": info.GetAltLoc(),
            "residue_name": info.GetResidueName(),
            "residue_number": info.GetResidueNumber(),
            "chain_id": info.GetChainId(),
            "insertion_code": info.GetInsertionCode(),
            "occupancy": info.GetOccupancy(),
            "temp_factor": info.GetTempFactor(),
            "is_hetero_atom": info.GetIsHeteroAtom(),
        },
    })
    coords.append([pos.x, pos.y, pos.z])

bonds = []
for bond in mol.GetBonds():
    bonds.append({
        "idx": bond.GetIdx(),
        "begin": bond.GetBeginAtomIdx(),
        "end": bond.GetEndAtomIdx(),
        "bond_type": str(bond.GetBondType()),
        "is_aromatic": bond.GetIsAromatic(),
        "direction": str(bond.GetBondDir()),
        "stereo": str(bond.GetStereo()),
        "stereo_atoms": list(bond.GetStereoAtoms()),
    })

print(json.dumps({
    "pdb_id": "7SH6",
    "source_url": "https://files.rcsb.org/download/7SH6.pdb",
    "atoms": atoms,
    "bonds": bonds,
    "coords": coords,
}, sort_keys=True))
"#;
    let output = Command::new(python_path())
        .arg("-c")
        .arg(script)
        .arg(path)
        .current_dir(repo_root())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .output()
        .expect("should run project Python RDKit oracle");
    assert!(
        output.status.success(),
        "failed to generate RDKit PDB oracle for {}:\nstdout:\n{}\nstderr:\n{}",
        path.display(),
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    serde_json::from_slice(&output.stdout).unwrap_or_else(|error| {
        panic!(
            "failed to parse RDKit PDB oracle JSON for {}: {error}\nstdout:\n{}",
            path.display(),
            String::from_utf8_lossy(&output.stdout)
        )
    })
}

fn bond_direction_name(direction: BondDirection) -> &'static str {
    match direction {
        BondDirection::None => "NONE",
        BondDirection::BeginWedge => "BEGINWEDGE",
        BondDirection::BeginDash => "BEGINDASH",
        BondDirection::EndUpRight => "ENDUPRIGHT",
        BondDirection::EndDownRight => "ENDDOWNRIGHT",
        BondDirection::EitherDouble => "EITHERDOUBLE",
        BondDirection::Unknown => "UNKNOWN",
    }
}

fn bond_stereo_name(stereo: BondStereo) -> &'static str {
    match stereo {
        BondStereo::None => "STEREONONE",
        BondStereo::Any => "STEREOANY",
        BondStereo::Z => "STEREOZ",
        BondStereo::E => "STEREOE",
        BondStereo::Cis => "STEREOCIS",
        BondStereo::Trans => "STEREOTRANS",
        BondStereo::AtropCw => "STEREOATROPCW",
        BondStereo::AtropCcw => "STEREOATROPCCW",
    }
}

fn assert_close(left: f64, right: f64, label: &str) {
    assert!((left - right).abs() <= 1.0e-5, "{label}: expected {right}, got {left}");
}

#[test]
fn pdb_7sh6_molfrompdbblock_state_matches_rdkit() {
    let pdb_path = ensure_7sh6_pdb();
    let oracle = rdkit_oracle_for_pdb(&pdb_path);
    assert_eq!(oracle.pdb_id, "7SH6");
    assert_eq!(oracle.source_url, "https://files.rcsb.org/download/7SH6.pdb");

    let pdb =
        fs::read_to_string(&pdb_path).unwrap_or_else(|error| panic!("failed to read {}: {error}", pdb_path.display()));
    let mol = Molecule::from_pdb_block(&pdb).unwrap_or_else(|error| panic!("COSMolKit PDB conversion failed: {error}"));

    assert_eq!(mol.num_atoms(), oracle.atoms.len(), "atom count mismatch");
    assert_eq!(mol.num_bonds(), oracle.bonds.len(), "bond count mismatch");

    let coords = mol
        .conformers_3d()
        .first()
        .expect("PDB MolFromPDBBlock result should have a 3D conformer")
        .coordinates();
    assert_eq!(coords.len(), oracle.coords.len(), "coordinate row mismatch");

    for (idx, (atom, expected)) in mol.atoms().iter().zip(&oracle.atoms).enumerate() {
        assert_eq!(expected.idx, idx, "RDKit atom index mismatch at row {idx}");
        assert_eq!(atom.id().index(), idx, "COSMolKit atom index mismatch at row {idx}");
        assert_eq!(
            atom.atomic_number(),
            expected.atomic_num,
            "atomic number mismatch at atom {idx}"
        );
        assert_eq!(
            atom.formal_charge(),
            expected.formal_charge,
            "formal charge mismatch at atom {idx}"
        );
        assert_eq!(atom.isotope(), expected.isotope, "isotope mismatch at atom {idx}");
        assert_eq!(
            atom.is_aromatic(),
            expected.is_aromatic,
            "atom aromaticity mismatch at atom {idx}"
        );
        assert_eq!(
            atom.chiral_tag().rdkit_name(),
            expected.chiral_tag,
            "chiral tag mismatch at atom {idx}"
        );
        assert_eq!(
            atom.radical_electrons(),
            expected.radical_electrons,
            "radical electron mismatch at atom {idx}"
        );

        let info = atom
            .pdb_residue_info()
            .unwrap_or_else(|| panic!("atom {idx} has no PDB residue info"));
        assert_eq!(
            info.atom_name(),
            expected.pdb.atom_name,
            "PDB atom name mismatch at atom {idx}"
        );
        assert_eq!(
            info.serial_number(),
            expected.pdb.serial_number,
            "PDB serial number mismatch at atom {idx}"
        );
        assert_eq!(
            info.alt_loc(),
            expected.pdb.alt_loc,
            "PDB altLoc mismatch at atom {idx}"
        );
        assert_eq!(
            info.residue_name(),
            expected.pdb.residue_name,
            "PDB residue name mismatch at atom {idx}"
        );
        assert_eq!(
            info.residue_number(),
            expected.pdb.residue_number,
            "PDB residue number mismatch at atom {idx}"
        );
        assert_eq!(
            info.chain_id(),
            expected.pdb.chain_id,
            "PDB chain id mismatch at atom {idx}"
        );
        assert_eq!(
            info.insertion_code(),
            expected.pdb.insertion_code,
            "PDB insertion code mismatch at atom {idx}"
        );
        assert_close(
            info.occupancy(),
            expected.pdb.occupancy,
            &format!("PDB occupancy mismatch at atom {idx}"),
        );
        assert_close(
            info.temp_factor(),
            expected.pdb.temp_factor,
            &format!("PDB temp factor mismatch at atom {idx}"),
        );
        assert_eq!(
            info.is_hetero_atom(),
            expected.pdb.is_hetero_atom,
            "PDB hetero flag mismatch at atom {idx}"
        );

        for axis in 0..3 {
            assert_close(
                coords[idx][axis],
                oracle.coords[idx][axis],
                &format!("3D coordinate mismatch at atom {idx} axis {axis}"),
            );
        }
    }

    for (idx, (bond, expected)) in mol.bonds().iter().zip(&oracle.bonds).enumerate() {
        assert_eq!(expected.idx, idx, "RDKit bond index mismatch at row {idx}");
        assert_eq!(
            bond.begin().index(),
            expected.begin,
            "bond begin mismatch at bond {idx}"
        );
        assert_eq!(bond.end().index(), expected.end, "bond end mismatch at bond {idx}");
        assert_eq!(
            bond.order().rdkit_name(),
            expected.bond_type,
            "bond type mismatch at bond {idx}"
        );
        assert_eq!(
            bond.is_aromatic(),
            expected.is_aromatic,
            "bond aromaticity mismatch at bond {idx}"
        );
        assert_eq!(
            bond_direction_name(bond.direction()),
            expected.direction,
            "bond direction mismatch at bond {idx}"
        );
        assert_eq!(
            bond_stereo_name(bond.stereo()),
            expected.stereo,
            "bond stereo mismatch at bond {idx}"
        );
        let stereo_atoms = bond
            .stereo_atoms()
            .map(|atoms| atoms.iter().map(|atom| atom.index()).collect::<Vec<_>>())
            .unwrap_or_default();
        assert_eq!(
            stereo_atoms, expected.stereo_atoms,
            "bond stereo atoms mismatch at bond {idx}"
        );
    }
}
