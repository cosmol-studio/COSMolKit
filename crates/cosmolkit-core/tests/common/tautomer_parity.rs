use std::collections::BTreeMap;

use cosmolkit_core::{
    Hybridization, Molecule, SmilesParseParams, TautomerCatalog, TautomerEnumerator, TautomerOptions, score_tautomer,
};
use serde::Deserialize;
use serde_json::{Value, json};

#[derive(Debug, Deserialize)]
pub struct GoldenRecord {
    pub row: usize,
    pub case_id: String,
    pub smiles: String,
    pub expected_canonical_smiles: Option<String>,
    pub sanitize: bool,
    pub remove_hs: bool,
    pub source: String,
    pub parse: Value,
    pub branches: BTreeMap<String, GoldenBranch>,
    pub input_tautomer: Option<GoldenEndpointInput>,
    pub atom_order_permutation: Option<GoldenEndpointInput>,
}

#[derive(Debug, Deserialize)]
pub struct GoldenBranch {
    pub parameters: Value,
    pub ok: bool,
    pub error: Option<Value>,
    pub status: Option<String>,
    pub ordered_smiles: Vec<String>,
    pub modified_atoms: Vec<usize>,
    pub modified_bonds: Vec<usize>,
    pub molecule_states: Vec<Value>,
    pub scores: Vec<Value>,
    pub canonical_smiles: Option<String>,
    pub canonical_state: Option<Value>,
}

#[derive(Debug, Deserialize)]
pub struct GoldenEndpointInput {
    pub smiles: String,
    pub parse: Value,
    pub branches: BTreeMap<String, GoldenEndpointBranch>,
}

#[derive(Debug, Deserialize)]
pub struct GoldenEndpointBranch {
    pub parameters: Value,
    pub ok: bool,
    pub error: Option<Value>,
    pub canonical_smiles: Option<String>,
    pub canonical_state: Option<Value>,
    pub canonical_score: Option<Value>,
}

fn hybridization_name(value: Hybridization) -> &'static str {
    match value {
        Hybridization::Unspecified => "UNSPECIFIED",
        Hybridization::S => "S",
        Hybridization::Sp => "SP",
        Hybridization::Sp2 => "SP2",
        Hybridization::Sp3 => "SP3",
        Hybridization::Sp2d => "SP2D",
        Hybridization::Sp3d => "SP3D",
        Hybridization::Sp3d2 => "SP3D2",
        Hybridization::Other => "OTHER",
    }
}

pub fn molecule_state(molecule: &Molecule) -> Value {
    let atoms = molecule
        .atoms()
        .iter()
        .map(|atom| {
            json!({
                "atomic_number": atom.atomic_number(),
                "formal_charge": atom.formal_charge(),
                "explicit_hydrogens": atom.explicit_hydrogens(),
                "no_implicit": atom.no_implicit(),
                "isotope": atom.isotope().unwrap_or(0),
                "radical_electrons": atom.radical_electrons(),
                "aromatic": atom.is_aromatic(),
                "chiral_tag": atom.chiral_tag().rdkit_name(),
                "hybridization": hybridization_name(atom.hybridization()),
                "cip_code": atom.prop("_CIPCode"),
            })
        })
        .collect::<Vec<_>>();
    let bonds = molecule
        .bonds()
        .iter()
        .map(|bond| {
            let stereo_atoms = bond
                .stereo_atoms()
                .map(|atoms| vec![atoms[0].index(), atoms[1].index()])
                .unwrap_or_default();
            json!({
                "begin": bond.begin().index(),
                "end": bond.end().index(),
                "bond_type": bond.order().rdkit_name(),
                "aromatic": bond.is_aromatic(),
                "conjugated": bond.is_conjugated(),
                "direction": bond.direction().rdkit_name(),
                "stereo": bond.stereo().rdkit_name(),
                "stereo_atoms": stereo_atoms,
            })
        })
        .collect::<Vec<_>>();
    json!({
        "isomeric_smiles": molecule.to_smiles(true).expect("serialize tautomer state"),
        "atoms": atoms,
        "bonds": bonds,
    })
}

pub fn configured_enumerator(parameters: &Value) -> TautomerEnumerator<'static> {
    let options = TautomerOptions::default()
        .with_max_tautomers(parameters["max_tautomers"].as_u64().unwrap() as u32)
        .with_max_transforms(parameters["max_transforms"].as_u64().unwrap() as u32)
        .with_remove_sp3_stereo(parameters["remove_sp3_stereo"].as_bool().unwrap())
        .with_remove_bond_stereo(parameters["remove_bond_stereo"].as_bool().unwrap())
        .with_remove_isotopic_hydrogens(parameters["remove_isotopic_hydrogens"].as_bool().unwrap())
        .with_reassign_stereo(parameters["reassign_stereo"].as_bool().unwrap());
    if parameters["catalog"] == "v1" {
        TautomerEnumerator::from_catalog_and_options(TautomerCatalog::v1().expect("compile V1 catalog"), options)
    } else {
        TautomerEnumerator::from_options(options)
    }
}

pub fn parse_record(record: &GoldenRecord) -> Result<Molecule, String> {
    let parameters = SmilesParseParams::with_sanitize(record.sanitize).with_remove_hs(record.remove_hs);
    Molecule::from_smiles_with_params(&record.smiles, &parameters).map_err(|error| error.to_string())
}

pub fn assert_branch(record: &GoldenRecord, name: &str, expected: &GoldenBranch, molecule: &Molecule) {
    let context = format!(
        "row {} case {} source {} branch {name}",
        record.row, record.case_id, record.source
    );
    assert_eq!(expected.parameters["name"], name, "{context}: profile identity");
    let enumerator = configured_enumerator(&expected.parameters);
    assert_eq!(
        enumerator.max_tautomers(),
        expected.parameters["max_tautomers"].as_u64().unwrap() as u32,
        "{context}: max tautomers option"
    );
    assert_eq!(
        enumerator.max_transforms(),
        expected.parameters["max_transforms"].as_u64().unwrap() as u32,
        "{context}: max transforms option"
    );
    assert_eq!(
        enumerator.remove_sp3_stereo(),
        expected.parameters["remove_sp3_stereo"].as_bool().unwrap(),
        "{context}: remove SP3 stereo option"
    );
    assert_eq!(
        enumerator.remove_bond_stereo(),
        expected.parameters["remove_bond_stereo"].as_bool().unwrap(),
        "{context}: remove bond stereo option"
    );
    assert_eq!(
        enumerator.remove_isotopic_hydrogens(),
        expected.parameters["remove_isotopic_hydrogens"].as_bool().unwrap(),
        "{context}: remove isotopic hydrogens option"
    );
    assert_eq!(
        enumerator.reassign_stereo(),
        expected.parameters["reassign_stereo"].as_bool().unwrap(),
        "{context}: reassign stereo option"
    );
    assert_eq!(
        enumerator.catalog().transforms().len(),
        if expected.parameters["catalog"] == "v1" { 36 } else { 37 },
        "{context}: catalog option"
    );
    let result = enumerator.enumerate(molecule);
    if !expected.ok {
        let error = result.expect_err(&format!("{context}: expected enumeration failure"));
        let expected_type = expected
            .error
            .as_ref()
            .and_then(|value| value["type"].as_str())
            .expect("failed golden branch must record its source error type");
        assert!(
            !error.to_string().is_empty(),
            "{context}: COSMolKit returned an empty error for RDKit {expected_type}"
        );
        return;
    }
    let result = result.unwrap_or_else(|error| panic!("{context}: enumeration failed: {error}"));
    assert_eq!(
        format!("{:?}", result.status()),
        expected.status.as_deref().unwrap(),
        "{context}: status"
    );
    assert_eq!(
        result.canonical_smiles(),
        expected.ordered_smiles,
        "{context}: ordered outputs"
    );
    assert_eq!(
        result.modified_atoms().iter().map(|id| id.index()).collect::<Vec<_>>(),
        expected.modified_atoms,
        "{context}: modified atoms"
    );
    assert_eq!(
        result.modified_bonds().iter().map(|id| id.index()).collect::<Vec<_>>(),
        expected.modified_bonds,
        "{context}: modified bonds"
    );
    let actual_states = result.iter().map(molecule_state).collect::<Vec<_>>();
    assert_eq!(actual_states, expected.molecule_states, "{context}: molecule states");
    let actual_scores = result
        .iter()
        .map(|tautomer| {
            let score = score_tautomer(tautomer).expect("score enumerated tautomer");
            json!({
                "ring": score.ring(),
                "substructure": score.substructure(),
                "hetero_hydrogen": score.hetero_hydrogen(),
                "total": score.total(),
            })
        })
        .collect::<Vec<_>>();
    assert_eq!(actual_scores, expected.scores, "{context}: score components");
    let canonical = enumerator
        .pick_canonical(&result)
        .unwrap_or_else(|error| panic!("{context}: canonical selection failed: {error}"));
    assert_eq!(
        canonical.to_smiles(true).expect("serialize canonical tautomer"),
        expected.canonical_smiles.as_deref().unwrap(),
        "{context}: canonical SMILES"
    );
    assert_eq!(
        molecule_state(&canonical),
        *expected.canonical_state.as_ref().unwrap(),
        "{context}: canonical state"
    );
}

pub fn assert_endpoint_input(row: usize, kind: &str, expected: &GoldenEndpointInput) {
    let context = format!("row {row} {kind} {}", expected.smiles);
    let parsed = Molecule::from_smiles(&expected.smiles);
    if !expected.parse["ok"].as_bool().unwrap() {
        assert!(parsed.is_err(), "{context}: expected parse failure");
        assert!(expected.branches.is_empty(), "{context}: failed parse branches");
        return;
    }
    let molecule = parsed.unwrap_or_else(|error| panic!("{context}: parse failed: {error}"));
    for (name, branch) in &expected.branches {
        assert_eq!(branch.parameters["name"], name.as_str(), "{context}: profile");
        let enumerator = configured_enumerator(&branch.parameters);
        let actual = enumerator.canonicalize(&molecule);
        if !branch.ok {
            let error = actual.expect_err(&format!("{context} {name}: expected failure"));
            let expected_type = branch
                .error
                .as_ref()
                .and_then(|value| value["type"].as_str())
                .expect("failed endpoint branch must record source error type");
            assert!(
                !error.to_string().is_empty(),
                "{context} {name}: empty error for RDKit {expected_type}"
            );
            continue;
        }
        let actual = actual.unwrap_or_else(|error| panic!("{context} {name}: {error}"));
        assert_eq!(
            actual.to_smiles(true).expect("serialize canonical endpoint"),
            branch.canonical_smiles.as_deref().unwrap(),
            "{context} {name}: canonical SMILES"
        );
        assert_eq!(
            molecule_state(&actual),
            *branch.canonical_state.as_ref().unwrap(),
            "{context} {name}: canonical state"
        );
        let score = score_tautomer(&actual).expect("score canonical endpoint");
        assert_eq!(
            json!({
                "ring": score.ring(),
                "substructure": score.substructure(),
                "hetero_hydrogen": score.hetero_hydrogen(),
                "total": score.total(),
            }),
            *branch.canonical_score.as_ref().unwrap(),
            "{context} {name}: canonical score"
        );
    }
}
