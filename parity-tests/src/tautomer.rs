//! Complete TAU observations using only the canonical public Molecule surface.
use crate::molecular::Outcome;
use crate::registry::molecule_plan::{TautomerCatalog, TautomerProfile};
use cosmolkit::{Molecule, TautomerParams};
use serde::{Deserialize, Serialize};
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AtomState {
    pub atomic_number: u8,
    pub formal_charge: i16,
    pub explicit_hydrogens: u16,
    pub no_implicit: bool,
    pub isotope: u16,
    pub radical_electrons: u8,
    pub aromatic: bool,
    pub chiral_tag: String,
    pub hybridization: String,
    pub cip_code: Option<String>,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BondState {
    pub begin: usize,
    pub end: usize,
    pub bond_type: String,
    pub aromatic: bool,
    pub conjugated: bool,
    pub direction: String,
    pub stereo: String,
    pub stereo_atoms: Vec<usize>,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MoleculeState {
    pub isomeric_smiles: String,
    pub atoms: Vec<AtomState>,
    pub bonds: Vec<BondState>,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Score {
    pub ring: i32,
    pub substructure: i32,
    pub hetero_hydrogen: i32,
    pub total: i32,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Status {
    Completed,
    MaxTautomersReached,
    MaxTransformsReached,
    Canceled,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Enumeration {
    pub ordered_smiles: Vec<String>,
    pub status: Status,
    pub modified_atoms: Vec<usize>,
    pub modified_bonds: Vec<usize>,
    pub scores: Vec<Score>,
    pub molecule_states: Vec<MoleculeState>,
    pub canonical_smiles: String,
    pub canonical_state: MoleculeState,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Canonicalization {
    pub canonical_smiles: String,
    pub canonical_state: MoleculeState,
    pub canonical_score: Score,
}
pub fn parameters(profile: TautomerProfile) -> Result<TautomerParams, String> {
    let params = match profile.catalog {
        TautomerCatalog::Current => TautomerParams::default(),
        TautomerCatalog::V1 => TautomerParams::v1().map_err(|e| e.to_string())?,
    };
    Ok(params
        .with_max_tautomers(profile.max_tautomers)
        .with_max_transforms(profile.max_transforms)
        .with_remove_sp3_stereo(profile.remove_sp3_stereo)
        .with_remove_bond_stereo(profile.remove_bond_stereo)
        .with_remove_isotopic_hydrogens(profile.remove_isotopic_hydrogens)
        .with_reassign_stereo(profile.reassign_stereo))
}
fn score(molecule: &Molecule) -> Result<Score, String> {
    let mut subject = molecule.clone();
    let score = subject.tautomer_score_().map_err(|e| e.to_string())?;
    Ok(Score {
        ring: score.ring(),
        substructure: score.substructure(),
        hetero_hydrogen: score.hetero_hydrogen(),
        total: score.total(),
    })
}
fn smiles(molecule: &Molecule) -> Result<String, String> {
    String::from_utf8(
        molecule
            .to_smiles()
            .map_err(|e| e.to_string())?
            .into_bytes(),
    )
    .map_err(|e| e.to_string())
}

fn state(molecule: &Molecule) -> Result<MoleculeState, String> {
    Ok(MoleculeState {
        isomeric_smiles: smiles(molecule)?,
        atoms: molecule
            .atoms()
            .iter()
            .map(|a| {
                Ok(AtomState {
                    atomic_number: a.atomic_number(),
                    formal_charge: i16::from(a.formal_charge()),
                    explicit_hydrogens: u16::from(a.explicit_hydrogens()),
                    no_implicit: a.no_implicit(),
                    isotope: a.isotope().unwrap_or(0),
                    radical_electrons: a.radical_electrons(),
                    aromatic: a.is_aromatic(),
                    chiral_tag: a.chiral_tag().rdkit_name().into(),
                    hybridization: a.hybridization().rdkit_name().into(),
                    cip_code: a
                        .prop("_CIPCode")
                        .map(|v| {
                            let text = v.as_string().map_err(|e| e.to_string())?;
                            std::str::from_utf8(text.as_bytes())
                                .map(str::to_owned)
                                .map_err(|e| e.to_string())
                        })
                        .transpose()?,
                })
            })
            .collect::<Result<_, String>>()?,
        bonds: molecule
            .bonds()
            .iter()
            .map(|b| BondState {
                begin: b.begin().index(),
                end: b.end().index(),
                bond_type: b.order().rdkit_name().into(),
                aromatic: b.is_aromatic(),
                conjugated: b.is_conjugated(),
                direction: b.direction().rdkit_name().into(),
                stereo: b
                    .stereo()
                    .rdkit_name()
                    .strip_prefix("STEREO")
                    .expect("source enum prefix")
                    .into(),
                stereo_atoms: b
                    .stereo_atoms()
                    .map(|ids| ids.map(|id| id.index()).to_vec())
                    .unwrap_or_default(),
            })
            .collect(),
    })
}
pub fn enumerate(molecule: &Molecule, profile: TautomerProfile) -> Result<Outcome, String> {
    let before = molecule.clone();
    let mut subject = molecule.clone();
    let mut result = subject
        .enumerate_tautomers_with_params_(&parameters(profile)?)
        .map_err(|e| e.to_string())?;
    let canonical = result.canonical_tautomer().map_err(|e| e.to_string())?;
    let value = Enumeration {
        ordered_smiles: result
            .canonical_smiles()
            .into_iter()
            .map(|text| std::str::from_utf8(text.as_bytes()).map(str::to_owned))
            .collect::<Result<_, _>>()
            .map_err(|e| e.to_string())?,
        status: match result.status() {
            cosmolkit::TautomerEnumerationStatus::Completed => Status::Completed,
            cosmolkit::TautomerEnumerationStatus::MaxTautomersReached => {
                Status::MaxTautomersReached
            }
            cosmolkit::TautomerEnumerationStatus::MaxTransformsReached => {
                Status::MaxTransformsReached
            }
            cosmolkit::TautomerEnumerationStatus::Canceled => Status::Canceled,
        },
        modified_atoms: result
            .modified_atoms()
            .iter()
            .map(|id| id.index())
            .collect(),
        modified_bonds: result
            .modified_bonds()
            .iter()
            .map(|id| id.index())
            .collect(),
        scores: result.iter().map(score).collect::<Result<_, _>>()?,
        molecule_states: result.iter().map(state).collect::<Result<_, _>>()?,
        canonical_smiles: smiles(&canonical)?,
        canonical_state: state(&canonical)?,
    };
    if *molecule != before {
        return Err("TAU value operation mutated its input".into());
    }
    Ok(Outcome::TautomerEnumeration(value))
}
pub fn canonicalize(molecule: &Molecule, profile: TautomerProfile) -> Result<Outcome, String> {
    let before = molecule.clone();
    let mut subject = molecule.clone();
    let canonical = subject
        .canonical_tautomer_with_params_(&parameters(profile)?)
        .map_err(|e| e.to_string())?;
    let value = Canonicalization {
        canonical_smiles: smiles(&canonical)?,
        canonical_state: state(&canonical)?,
        canonical_score: score(&canonical)?,
    };
    if *molecule != before {
        return Err("TAU value operation mutated its input".into());
    }
    Ok(Outcome::TautomerCanonicalization(value))
}
fn valid_state(state: &MoleculeState) -> bool {
    state.bonds.iter().all(|b| {
        b.begin < state.atoms.len()
            && b.end < state.atoms.len()
            && b.stereo_atoms.iter().all(|&a| a < state.atoms.len())
            && (b.stereo_atoms.is_empty() || b.stereo_atoms.len() == 2)
    })
}
pub fn valid_enumeration(value: &Enumeration) -> bool {
    let len = value.ordered_smiles.len();
    len > 0
        && len == value.molecule_states.len()
        && len == value.scores.len()
        && value
            .ordered_smiles
            .windows(2)
            .all(|pair| pair[0] < pair[1])
        && value
            .modified_atoms
            .windows(2)
            .all(|pair| pair[0] < pair[1])
        && value
            .modified_bonds
            .windows(2)
            .all(|pair| pair[0] < pair[1])
        && value.molecule_states.iter().all(valid_state)
        && valid_state(&value.canonical_state)
        && value
            .modified_atoms
            .iter()
            .all(|&id| id < value.canonical_state.atoms.len())
        && value
            .modified_bonds
            .iter()
            .all(|&id| id < value.canonical_state.bonds.len())
        && value
            .scores
            .iter()
            .all(|s| s.total == s.ring + s.substructure + s.hetero_hydrogen)
        && value.canonical_smiles == value.canonical_state.isomeric_smiles
}
pub fn valid_canonicalization(value: &Canonicalization) -> bool {
    valid_state(&value.canonical_state)
        && value.canonical_smiles == value.canonical_state.isomeric_smiles
        && value.canonical_score.total
            == value.canonical_score.ring
                + value.canonical_score.substructure
                + value.canonical_score.hetero_hydrogen
}
