#![cfg(feature = "cap-alignment")]
use std::{
    fs::File,
    io::{BufRead, BufReader},
};

use cosmolkit::{
    AlignmentAtomMap, AlignmentError, AlignmentParameters, AllConformerRmsdParameters,
    BestAlignmentParameters, Conformer3D, ConformerAlignmentParameters, CoordinateRmsdParameters,
    Molecule,
};
use serde::Deserialize;
use serde_json::Value;

const OUTPUT_NAME: &str = "molalign.jsonl";
const NUMERICAL_TOLERANCE: f64 = 1.0e-8;

#[derive(Debug, Clone, Deserialize)]
struct ConformerSource {
    id: usize,
    coordinates: Vec<[f64; 3]>,
}

#[derive(Debug, Clone, Deserialize)]
struct MoleculeSource {
    smiles: String,
    conformers: Vec<ConformerSource>,
}

#[derive(Debug, Clone, Deserialize)]
struct OracleSource {
    probe: Option<MoleculeSource>,
    reference: Option<MoleculeSource>,
    molecule: Option<MoleculeSource>,
    input_smiles: Option<String>,
}

#[derive(Debug, Clone, Deserialize)]
struct ConformerState {
    id: usize,
    is_3d: bool,
    coordinates: Vec<[f64; 3]>,
}

#[derive(Debug, Clone, Deserialize)]
struct ObservableState {
    probe: Option<Vec<ConformerState>>,
    reference: Option<Vec<ConformerState>>,
    molecule: Option<Vec<ConformerState>>,
}

#[derive(Debug, Clone, Deserialize)]
struct OracleRecord {
    schema_version: u32,
    case_id: String,
    call_index: usize,
    operation: String,
    source: OracleSource,
    parameters: Value,
    status: String,
    result: Option<Value>,
    error_kind: Option<String>,
    error_type: Option<String>,
    error_message: Option<String>,
    before: ObservableState,
    after: ObservableState,
}

fn records(profile: &str) -> Vec<OracleRecord> {
    let path = std::env::var_os("ALIGNMENT_ORACLE_ROOT")
        .map(std::path::PathBuf::from)
        .expect("prepared pinned-source oracle root")
        .join(profile)
        .join(OUTPUT_NAME);
    BufReader::new(File::open(&path).expect("open MolAlign oracle"))
        .lines()
        .enumerate()
        .map(|(line_index, line)| {
            serde_json::from_str(&line.expect("read MolAlign oracle row"))
                .unwrap_or_else(|error| panic!("{}:{}: {error}", path.display(), line_index + 1))
        })
        .collect()
}

fn molecule(source: &MoleculeSource) -> Molecule {
    let molecule = Molecule::from_smiles(&source.smiles).expect("parse oracle molecule");
    let coordinates = cosmolkit::CoordinateBlock {
        conformers_3d: source
            .conformers
            .iter()
            .map(|c| Conformer3D::new(c.id, c.coordinates.clone(), true))
            .collect(),
        ..Default::default()
    };
    Molecule::from_parts(
        molecule.topology().clone(),
        coordinates,
        molecule.properties().clone(),
    )
    .expect("checked sparse-ID oracle fixture")
}

fn number(parameters: &Value, name: &str, default: i64) -> i64 {
    parameters
        .get(name)
        .and_then(Value::as_i64)
        .unwrap_or(default)
}

fn flag(parameters: &Value, name: &str, default: bool) -> bool {
    parameters
        .get(name)
        .and_then(Value::as_bool)
        .unwrap_or(default)
}

fn indices(parameters: &Value, name: &str) -> Option<Vec<usize>> {
    parameters.get(name).map(|values| {
        values
            .as_array()
            .expect("index list")
            .iter()
            .map(|value| value.as_u64().expect("unsigned index") as usize)
            .collect()
    })
}

fn weights(parameters: &Value) -> Option<Vec<f64>> {
    parameters.get("weights").map(|values| {
        values
            .as_array()
            .expect("weight list")
            .iter()
            .map(|value| value.as_f64().expect("floating-point weight"))
            .collect()
    })
}

fn atom_map(value: &Value) -> Vec<AlignmentAtomMap> {
    value
        .as_array()
        .expect("atom map")
        .iter()
        .map(|pair| {
            let pair = pair.as_array().expect("atom-map pair");
            AlignmentAtomMap {
                probe_atom: pair[0].as_u64().expect("probe atom") as usize,
                reference_atom: pair[1].as_u64().expect("reference atom") as usize,
            }
        })
        .collect()
}

fn atom_maps(parameters: &Value) -> Vec<Vec<AlignmentAtomMap>> {
    parameters
        .get("atom_maps")
        .map(|maps| {
            maps.as_array()
                .expect("atom maps")
                .iter()
                .map(atom_map)
                .collect()
        })
        .unwrap_or_default()
}

fn alignment_parameters(parameters: &Value) -> AlignmentParameters {
    AlignmentParameters {
        probe_conformer_id: number(parameters, "probe_conformer_id", -1) as i32,
        reference_conformer_id: number(parameters, "reference_conformer_id", -1) as i32,
        atom_map: parameters.get("atom_map").map(atom_map),
        weights: weights(parameters),
        reflect: flag(parameters, "reflect", false),
        max_iterations: number(parameters, "max_iterations", 50) as u32,
    }
}

fn best_parameters(parameters: &Value) -> BestAlignmentParameters {
    BestAlignmentParameters {
        probe_conformer_id: number(parameters, "probe_conformer_id", -1) as i32,
        reference_conformer_id: number(parameters, "reference_conformer_id", -1) as i32,
        atom_maps: atom_maps(parameters),
        weights: weights(parameters),
        reflect: flag(parameters, "reflect", false),
        max_iterations: number(parameters, "max_iterations", 50) as u32,
        max_matches: number(parameters, "max_matches", 1_000_000) as i32,
        symmetrize_conjugated_terminal_groups: flag(
            parameters,
            "symmetrize_conjugated_terminal_groups",
            true,
        ),
        ignore_hydrogens: flag(parameters, "ignore_hydrogens", true),
        num_threads: number(parameters, "num_threads", 1) as i32,
    }
}

fn all_conformer_parameters(parameters: &Value) -> AllConformerRmsdParameters {
    AllConformerRmsdParameters {
        atom_maps: atom_maps(parameters),
        weights: weights(parameters),
        max_matches: number(parameters, "max_matches", 1_000_000) as i32,
        symmetrize_conjugated_terminal_groups: flag(
            parameters,
            "symmetrize_conjugated_terminal_groups",
            true,
        ),
        ignore_hydrogens: flag(parameters, "ignore_hydrogens", true),
        num_threads: number(parameters, "num_threads", 1) as i32,
    }
}

fn coordinate_parameters(parameters: &Value) -> CoordinateRmsdParameters {
    CoordinateRmsdParameters {
        probe_conformer_id: number(parameters, "probe_conformer_id", -1) as i32,
        reference_conformer_id: number(parameters, "reference_conformer_id", -1) as i32,
        atom_maps: atom_maps(parameters),
        weights: weights(parameters),
        max_matches: number(parameters, "max_matches", 1_000_000) as i32,
        symmetrize_conjugated_terminal_groups: flag(
            parameters,
            "symmetrize_conjugated_terminal_groups",
            true,
        ),
    }
}

fn conformer_parameters(parameters: &Value) -> ConformerAlignmentParameters {
    ConformerAlignmentParameters {
        atom_indices: indices(parameters, "atom_indices"),
        conformer_ids: indices(parameters, "conformer_ids"),
        weights: weights(parameters),
        reflect: flag(parameters, "reflect", false),
        max_iterations: number(parameters, "max_iterations", 50) as u32,
    }
}

fn assert_close(actual: f64, expected: f64, context: &str) {
    assert!(
        (actual - expected).abs() <= NUMERICAL_TOLERANCE,
        "{context}: actual={actual:.17}, expected={expected:.17}, tolerance={NUMERICAL_TOLERANCE}"
    );
}

fn assert_matrix(actual: &[[f64; 4]; 4], expected: &Value, context: &str) {
    let rows = expected.as_array().expect("transform rows");
    for row in 0..4 {
        let columns = rows[row].as_array().expect("transform columns");
        for column in 0..4 {
            assert_close(
                actual[row][column],
                columns[column].as_f64().expect("transform value"),
                context,
            );
        }
    }
}

fn snapshot(molecule: &Molecule) -> Vec<ConformerState> {
    molecule
        .conformers_3d()
        .iter()
        .map(|conformer| ConformerState {
            id: conformer.id(),
            is_3d: conformer.is_3d(),
            coordinates: conformer.coordinates().to_vec(),
        })
        .collect()
}

fn assert_state(actual: &[ConformerState], expected: &[ConformerState], context: &str) {
    assert_eq!(actual.len(), expected.len(), "{context}: conformer count");
    for (actual, expected) in actual.iter().zip(expected) {
        assert_eq!(actual.id, expected.id, "{context}: conformer id");
        assert_eq!(actual.is_3d, expected.is_3d, "{context}: dimensionality");
        assert_eq!(
            actual.coordinates.len(),
            expected.coordinates.len(),
            "{context}: coordinate count"
        );
        for (actual, expected) in actual.coordinates.iter().zip(&expected.coordinates) {
            for axis in 0..3 {
                assert_close(actual[axis], expected[axis], context);
            }
        }
    }
}

fn assert_error(actual: AlignmentError, expected: &str, context: &str) {
    let matches = matches!(
        (expected, &actual),
        (
            "conformer_not_found",
            AlignmentError::ConformerNotFound { .. }
        ) | (
            "weight_count_mismatch",
            AlignmentError::WeightCountMismatch { .. }
        ) | ("no_substructure_match", AlignmentError::NoSubstructureMatch)
    );
    assert!(matches, "{context}: expected {expected}, got {actual:?}");
}

fn assert_profile(profile: &str, expected_records: usize) {
    let records = records(profile);
    assert_eq!(records.len(), expected_records, "{profile}: record count");
    for record in records {
        assert_eq!(record.schema_version, 1);
        let context = format!(
            "{} call {} ({})",
            record.case_id, record.call_index, record.operation
        );
        assert!(
            record.error_type.is_none() == (record.status == "ok"),
            "{context}"
        );
        assert!(
            record.error_message.is_none() == (record.status == "ok"),
            "{context}"
        );

        if record.operation == "input_parse" {
            assert_eq!(record.status, "error", "{context}");
            assert_eq!(record.error_kind.as_deref(), Some("input_parse_error"));
            let smiles = record
                .source
                .input_smiles
                .as_deref()
                .expect("input parse record SMILES");
            assert!(
                Molecule::from_smiles(smiles).is_err(),
                "{context}: COSMolKit accepted an input rejected by pinned RDKit"
            );
            continue;
        }

        if let (Some(probe_source), Some(reference_source)) =
            (&record.source.probe, &record.source.reference)
        {
            let probe = molecule(probe_source);
            let reference = molecule(reference_source);
            assert_state(
                &snapshot(&probe),
                record.before.probe.as_deref().expect("probe before state"),
                &context,
            );
            assert_state(
                &snapshot(&reference),
                record
                    .before
                    .reference
                    .as_deref()
                    .expect("reference before state"),
                &context,
            );
            if record.operation == "align_to" {
                let (aligned, actual) = probe
                    .with_alignment_to_with_params(
                        &reference,
                        &alignment_parameters(&record.parameters),
                    )
                    .expect("value-style alignment");
                let expected = record.result.as_ref().expect("oracle result");
                assert_close(
                    actual.rmsd,
                    expected["rmsd"].as_f64().expect("oracle RMSD"),
                    &context,
                );
                assert_matrix(&actual.transform.matrix, &expected["transform"], &context);
                assert_eq!(
                    actual.atom_map,
                    atom_map(&expected["atom_map"]),
                    "{context}: applied atom map"
                );
                assert_state(
                    &snapshot(&aligned),
                    record.after.probe.as_deref().expect("aligned probe state"),
                    &context,
                );
                assert_state(
                    &snapshot(&probe),
                    record.before.probe.as_deref().expect("source probe state"),
                    &context,
                );
                assert_state(
                    &snapshot(&reference),
                    record
                        .after
                        .reference
                        .as_deref()
                        .expect("reference after state"),
                    &context,
                );
                continue;
            }
            let result = match record.operation.as_str() {
                "alignment_transform" => probe.alignment_transform_to_with_params(
                    &reference,
                    &alignment_parameters(&record.parameters),
                ),
                "best_alignment" => probe.best_alignment_to_with_params(
                    &reference,
                    &best_parameters(&record.parameters),
                ),
                "coordinate_rmsd" => {
                    let actual = probe.coordinate_rmsd_to_with_params(
                        &reference,
                        &coordinate_parameters(&record.parameters),
                    );
                    match actual {
                        Ok(rmsd) => {
                            let expected = record.result.as_ref().expect("oracle result");
                            assert_close(
                                rmsd,
                                expected["rmsd"].as_f64().expect("oracle RMSD"),
                                &context,
                            );
                            assert_eq!(record.status, "ok", "{context}");
                        }
                        Err(error) => assert_error(
                            error,
                            record.error_kind.as_deref().expect("oracle error kind"),
                            &context,
                        ),
                    }
                    assert_state(
                        &snapshot(&probe),
                        record.after.probe.as_deref().expect("probe after state"),
                        &context,
                    );
                    assert_state(
                        &snapshot(&reference),
                        record
                            .after
                            .reference
                            .as_deref()
                            .expect("reference after state"),
                        &context,
                    );
                    continue;
                }
                other => panic!("{context}: unsupported operation {other}"),
            };
            match result {
                Ok(actual) => {
                    assert_eq!(record.status, "ok", "{context}");
                    let expected = record.result.as_ref().expect("oracle result");
                    assert_close(
                        actual.rmsd,
                        expected["rmsd"].as_f64().expect("oracle RMSD"),
                        &context,
                    );
                    assert_matrix(&actual.transform.matrix, &expected["transform"], &context);
                    assert_eq!(
                        actual.atom_map,
                        atom_map(&expected["atom_map"]),
                        "{context}: selected atom map"
                    );
                }
                Err(error) => assert_error(
                    error,
                    record.error_kind.as_deref().expect("oracle error kind"),
                    &context,
                ),
            }
            assert_state(
                &snapshot(&probe),
                record.after.probe.as_deref().expect("probe after state"),
                &context,
            );
            assert_state(
                &snapshot(&reference),
                record
                    .after
                    .reference
                    .as_deref()
                    .expect("reference after state"),
                &context,
            );
            continue;
        }

        let source = record.source.molecule.as_ref().expect("molecule source");
        let molecule = molecule(source);
        assert_state(
            &snapshot(&molecule),
            record
                .before
                .molecule
                .as_deref()
                .expect("molecule before state"),
            &context,
        );
        match record.operation.as_str() {
            "all_conformer_best_rms" => {
                let actual = molecule
                    .all_conformer_best_rmsds_with_params(&all_conformer_parameters(
                        &record.parameters,
                    ))
                    .expect("all-conformer RMSD");
                let expected = record.result.as_ref().expect("oracle result");
                let rmsds = expected["rmsds"].as_array().expect("oracle RMSD list");
                let pairs = expected["conformer_pairs"]
                    .as_array()
                    .expect("oracle conformer pairs");
                assert_eq!(actual.len(), rmsds.len(), "{context}");
                for (index, actual) in actual.iter().enumerate() {
                    let pair = pairs[index].as_array().expect("oracle conformer pair");
                    assert_eq!(
                        [actual.probe_conformer_id, actual.reference_conformer_id],
                        [
                            pair[0].as_u64().expect("probe conformer id") as usize,
                            pair[1].as_u64().expect("reference conformer id") as usize,
                        ],
                        "{context}: triangular pair order"
                    );
                    assert_close(
                        actual.rmsd,
                        rmsds[index].as_f64().expect("oracle RMSD"),
                        &context,
                    );
                }
                assert_state(
                    &snapshot(&molecule),
                    record
                        .after
                        .molecule
                        .as_deref()
                        .expect("molecule after state"),
                    &context,
                );
            }
            "align_conformers" => {
                let (aligned, report) = molecule
                    .with_aligned_conformers_with_params(&conformer_parameters(&record.parameters))
                    .expect("value-style conformer alignment");
                let expected = record.result.as_ref().expect("oracle result");
                let rmsds = expected["rmsds"].as_array().expect("oracle RMSD list");
                assert_eq!(report.rmsds.len(), rmsds.len(), "{context}");
                for (actual, expected) in report.rmsds.iter().zip(rmsds) {
                    assert_close(*actual, expected.as_f64().expect("oracle RMSD"), &context);
                }
                assert_state(
                    &snapshot(&aligned),
                    record
                        .after
                        .molecule
                        .as_deref()
                        .expect("molecule after state"),
                    &context,
                );
                assert_state(
                    &snapshot(&molecule),
                    record
                        .before
                        .molecule
                        .as_deref()
                        .expect("source remains unchanged"),
                    &context,
                );
            }
            other => panic!("{context}: unsupported operation {other}"),
        }
    }
}

#[test]
fn focused_molalign_matches_every_pinned_rdkit_observable_field() {
    assert_profile("molalign_focused", 14);
}

#[test]
fn curated_molalign_matches_every_pinned_rdkit_observable_field() {
    assert_profile("smiles_small", 152);
}

#[test]
#[ignore = "run explicitly after preparing the smiles_5000 alignment profile"]
fn molalign_matches_all_5000_rows_and_rotated_operation_branches() {
    assert_profile("smiles_5000", 5_000);
}

// Test-data preparation uses only the public checked constructor and binary API.
// It retains all source conformer IDs and f64 rows; this is not an alignment path.
#[test]
#[ignore = "explicit canonical sparse-ID Python fixture preparation"]
fn prepare_canonical_python_alignment_fixtures() {
    use std::io::Write;
    let root = std::path::PathBuf::from(std::env::var_os("ALIGNMENT_ORACLE_ROOT").unwrap());
    let mut output = std::io::BufWriter::new(
        File::create(root.parent().unwrap().join("canonical-input-pickles.jsonl")).unwrap(),
    );
    let mut count = 0;
    for name in ["molalign_focused", "smiles_small", "smiles_5000"] {
        for record in records(name) {
            let mut inputs = serde_json::Map::new();
            for (key, source) in [
                ("probe", record.source.probe),
                ("reference", record.source.reference),
                ("molecule", record.source.molecule),
            ] {
                if let Some(source) = source {
                    let m = molecule(&source);
                    let bytes = m.to_binary().expect("canonical binary fixture");
                    let restored =
                        Molecule::from_binary(&bytes).expect("canonical binary fixture round trip");
                    assert_state(
                        &snapshot(&restored),
                        &snapshot(&m),
                        &format!("{} fixture round trip", record.case_id),
                    );
                    inputs.insert(key.into(), serde_json::to_value(bytes).unwrap());
                }
            }
            serde_json::to_writer(&mut output,&serde_json::json!({"profile":name,"case_id":record.case_id,"call_index":record.call_index,"inputs":inputs})).unwrap();
            writeln!(&mut output).unwrap();
            count += 1;
        }
    }
    output.flush().unwrap();
    assert_eq!(count, 5166);
}
