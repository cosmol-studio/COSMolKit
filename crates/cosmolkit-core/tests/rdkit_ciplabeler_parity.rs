use std::{
    fs::File,
    io::{BufRead, BufReader},
};

use cosmolkit_core::{
    AtomId, BondId, BondStereo, CipLabelOptions, Molecule, OperationError, SdfReadParams, SmilesParseParams,
};
use serde::Deserialize;
use serde_json::{Value, json};

mod common;
use common::parity_data;

const OUTPUT_NAME: &str = "ciplabeler.jsonl";

#[derive(Debug, Clone, Deserialize, PartialEq, Eq)]
struct PropertyState {
    present: bool,
    value: Value,
    value_type: Option<String>,
    computed: bool,
}

#[derive(Debug, Clone, Deserialize, PartialEq, Eq)]
struct MoleculeState {
    cip_computed: PropertyState,
}

#[derive(Debug, Clone, Deserialize, PartialEq, Eq)]
struct AtomState {
    index: usize,
    chiral_tag: String,
    cip_code: PropertyState,
    cip_neighbor_order: PropertyState,
    cip_rank: PropertyState,
}

#[derive(Debug, Clone, Deserialize, PartialEq, Eq)]
struct BondState {
    index: usize,
    begin: usize,
    end: usize,
    stereo: String,
    stereo_atoms_u32: Vec<u64>,
    cip_code: PropertyState,
    cip_neighbor_order: PropertyState,
}

#[derive(Debug, Clone, Deserialize, PartialEq, Eq)]
struct ObservableState {
    molecule: MoleculeState,
    atoms: Vec<AtomState>,
    bonds: Vec<BondState>,
}

#[derive(Debug, Clone, Deserialize, PartialEq, Eq)]
struct CallResult {
    status: String,
    error_type: Option<String>,
    error_message: Option<String>,
}

#[derive(Debug, Clone, Deserialize)]
struct OracleCall {
    atoms_to_label: Option<Vec<usize>>,
    bonds_to_label: Option<Vec<usize>>,
    max_recursive_iterations: u32,
    result: CallResult,
    state: ObservableState,
}

#[derive(Debug, Clone, Deserialize)]
struct OracleSource {
    fixture: String,
    input_kind: String,
    input: String,
    sanitize: bool,
    remove_hs: bool,
    strict_parsing: bool,
    update_property_cache: bool,
    strict_property_cache: bool,
}

#[derive(Debug, Clone, Deserialize)]
struct OracleRecord {
    schema_version: u32,
    case_id: String,
    source: OracleSource,
    parse_status: String,
    parse_error: Option<String>,
    initial_state: Option<ObservableState>,
    calls: Vec<OracleCall>,
}

fn load_records(profile: &str) -> Vec<OracleRecord> {
    let path = parity_data::expected_path_for_profile("stereo", "rdkit", profile, OUTPUT_NAME);
    let file = File::open(&path).unwrap_or_else(|error| {
        panic!(
            "failed to open {}; regenerate it with `uv run python tools/testdata/rdkit/generate_all.py --only ciplabeler --profile {profile}`: {error}",
            path.display()
        )
    });
    BufReader::new(file)
        .lines()
        .enumerate()
        .map(|(line_index, line)| {
            let line = line
                .unwrap_or_else(|error| panic!("failed to read {} line {}: {error}", path.display(), line_index + 1));
            serde_json::from_str(&line)
                .unwrap_or_else(|error| panic!("failed to parse {} line {}: {error}", path.display(), line_index + 1))
        })
        .collect()
}

fn property_state(value: Option<&str>, computed: bool, value_type: &'static str) -> PropertyState {
    let Some(value) = value else {
        return PropertyState {
            present: false,
            value: Value::Null,
            value_type: None,
            computed: false,
        };
    };
    let value = match value_type {
        "boolean" => json!(match value {
            "0" => false,
            "1" => true,
            other => panic!("invalid stored boolean property `{other}`"),
        }),
        "unsigned" => json!(
            value
                .parse::<u64>()
                .unwrap_or_else(|error| panic!("invalid stored unsigned property `{value}`: {error}"))
        ),
        "unsigned_vector" => serde_json::from_str(value)
            .unwrap_or_else(|error| panic!("invalid stored unsigned-vector property `{value}`: {error}")),
        "string" => json!(value),
        other => panic!("unsupported oracle property type `{other}`"),
    };
    PropertyState {
        present: true,
        value,
        value_type: Some(value_type.to_owned()),
        computed,
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

fn snapshot(molecule: &Molecule) -> ObservableState {
    ObservableState {
        molecule: MoleculeState {
            cip_computed: property_state(
                molecule.prop("_CIPComputed"),
                molecule.is_prop_computed("_CIPComputed"),
                "boolean",
            ),
        },
        atoms: molecule
            .atoms()
            .iter()
            .map(|atom| AtomState {
                index: atom.id().index(),
                chiral_tag: atom.chiral_tag().rdkit_name().to_owned(),
                cip_code: property_state(atom.prop("_CIPCode"), atom.is_prop_computed("_CIPCode"), "string"),
                cip_neighbor_order: property_state(
                    atom.prop("_CIPNeighborOrder"),
                    atom.is_prop_computed("_CIPNeighborOrder"),
                    "unsigned_vector",
                ),
                cip_rank: property_state(atom.prop("_CIPRank"), atom.is_prop_computed("_CIPRank"), "unsigned"),
            })
            .collect(),
        bonds: molecule
            .bonds()
            .iter()
            .map(|bond| BondState {
                index: bond.id().index(),
                begin: bond.begin().index(),
                end: bond.end().index(),
                stereo: bond_stereo_name(bond.stereo()).to_owned(),
                stereo_atoms_u32: bond.stereo_atoms().map_or_else(Vec::new, |atoms| {
                    atoms.into_iter().map(|atom| atom.index() as u64).collect()
                }),
                cip_code: property_state(bond.prop("_CIPCode"), bond.is_prop_computed("_CIPCode"), "string"),
                cip_neighbor_order: property_state(
                    bond.prop("_CIPNeighborOrder"),
                    bond.is_prop_computed("_CIPNeighborOrder"),
                    "unsigned_vector",
                ),
            })
            .collect(),
    }
}

fn parse_source(source: &OracleSource) -> Result<Molecule, String> {
    let mut molecule = match source.input_kind.as_str() {
        "smiles" => Molecule::from_smiles_with_params(
            &source.input,
            &SmilesParseParams::with_sanitize(source.sanitize).with_remove_hs(source.remove_hs),
        )
        .map_err(|error| error.to_string())?,
        "molblock" => Molecule::from_mol_block_with_params(
            &source.input,
            SdfReadParams {
                sanitize: source.sanitize,
                remove_hs: source.remove_hs,
                strict_parsing: source.strict_parsing,
                ..SdfReadParams::default()
            },
        )
        .map_err(|error| error.to_string())?,
        other => return Err(format!("unsupported focused CIP input kind `{other}`")),
    };
    if source.update_property_cache {
        molecule = molecule
            .with_assigned_valence_strict(source.strict_property_cache)
            .map_err(|error| error.to_string())?;
    }
    Ok(molecule)
}

fn options_for(call: &OracleCall) -> CipLabelOptions {
    let mut options = CipLabelOptions::default().with_max_recursive_iterations(call.max_recursive_iterations);
    if let Some(indices) = &call.atoms_to_label {
        options = options.with_atoms(indices.iter().copied().map(AtomId::new));
    }
    if let Some(indices) = &call.bonds_to_label {
        options = options.with_bonds(indices.iter().copied().map(BondId::new));
    }
    options
}

fn call_result(result: Result<(), OperationError>) -> CallResult {
    match result {
        Ok(()) => CallResult {
            status: "ok".to_owned(),
            error_type: None,
            error_message: None,
        },
        Err(OperationError::CipLabeler { source, .. }) => CallResult {
            status: "error".to_owned(),
            error_type: Some("RuntimeError".to_owned()),
            error_message: Some(source.to_string()),
        },
        Err(error) => panic!("modern CIP operation returned a non-CIP error: {error}"),
    }
}

fn assert_state(actual: &ObservableState, expected: &ObservableState, case_id: &str, stage: &str) {
    assert_eq!(
        actual, expected,
        "complete modern CIP observable state mismatch for `{case_id}` at {stage}"
    );
}

fn read_smiles(profile: &str) -> Option<Vec<String>> {
    let path = match profile {
        "ciplabeler_focused" => return None,
        "smiles_small" => parity_data::repo_root().join("testdata/smiles/corpus/smiles_small.smi"),
        "smiles_5000" => parity_data::repo_root().join("testdata/smiles/corpus/smiles_5000.smi"),
        other => panic!("unsupported modern CIP parity profile `{other}`"),
    };
    Some(
        std::fs::read_to_string(&path)
            .unwrap_or_else(|error| panic!("failed to read {}: {error}", path.display()))
            .lines()
            .map(str::trim)
            .filter(|line| !line.is_empty() && !line.starts_with('#'))
            .map(str::to_owned)
            .collect(),
    )
}

fn assert_profile(profile: &str, expected_count: usize) {
    let records = load_records(profile);
    assert_eq!(records.len(), expected_count, "{profile} record count changed");
    let corpus = read_smiles(profile);
    if let Some(corpus) = &corpus {
        assert_eq!(records.len(), corpus.len(), "{profile} corpus coverage");
    }

    for (row, record) in records.into_iter().enumerate() {
        assert_eq!(record.schema_version, 1, "{}", record.case_id);
        assert!(
            !record.source.fixture.is_empty(),
            "{} has no fixture provenance",
            record.case_id
        );
        if let Some(corpus) = &corpus {
            assert_eq!(record.case_id, format!("corpus-{row}"), "{profile} row id");
            assert_eq!(record.source.input, corpus[row], "{profile} row {row}");
            assert_eq!(record.source.input_kind, "smiles", "{profile} row {row}");
        }
        let parsed = parse_source(&record.source);
        match record.parse_status.as_str() {
            "ok" => {
                assert_eq!(record.parse_error, None, "{}", record.case_id);
            }
            "error" => {
                let actual_error = parsed.expect_err("oracle expected a parse error");
                assert_eq!(Some(actual_error), record.parse_error, "{}", record.case_id);
                continue;
            }
            "none" => {
                assert!(parsed.is_err(), "{} unexpectedly parsed", record.case_id);
                assert_eq!(record.parse_error, None, "{}", record.case_id);
                continue;
            }
            other => panic!("{} has invalid parse status `{other}`", record.case_id),
        }

        let mut molecule = parsed.unwrap_or_else(|error| panic!("{} failed to parse: {error}", record.case_id));
        let initial_state = record
            .initial_state
            .as_ref()
            .unwrap_or_else(|| panic!("{} has no initial state", record.case_id));
        assert_state(&snapshot(&molecule), initial_state, &record.case_id, "parse");

        for (call_index, call) in record.calls.iter().enumerate() {
            let actual_result = call_result(molecule.assign_cip_labels_with_options_(options_for(call)));
            assert_eq!(
                actual_result, call.result,
                "modern CIP result mismatch for `{}` call {}",
                record.case_id, call_index
            );
            assert_state(
                &snapshot(&molecule),
                &call.state,
                &record.case_id,
                &format!("call {call_index}"),
            );
        }
    }
}

#[test]
fn focused_modern_ciplabeler_matches_every_pinned_rdkit_observable_field() {
    assert_profile("ciplabeler_focused", 19);
}

#[test]
fn maintained_modern_ciplabeler_profile_matches_every_pinned_rdkit_observable_field() {
    match parity_data::profile_name().as_str() {
        "smiles_small" => assert_profile("smiles_small", 152),
        "smiles_5000" => assert_profile("smiles_5000", 5_000),
        profile => panic!("unsupported modern CIP parity profile `{profile}`"),
    }
}

#[test]
#[ignore = "run explicitly with COSMOLKIT_PARITY_PROFILE=smiles_5000 in release mode"]
fn modern_ciplabeler_matches_all_5000_rows_exactly() {
    assert_eq!(
        parity_data::profile_name(),
        "smiles_5000",
        "set COSMOLKIT_PARITY_PROFILE=smiles_5000"
    );
    assert_profile("smiles_5000", 5_000);
}
