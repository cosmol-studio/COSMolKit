//! Public MolAlign comparison adapters, preserving the original field checks.
use cosmolkit::{
    AlignmentAtomMap, AlignmentError, AlignmentParameters, AllConformerRmsdParameters,
    BestAlignmentParameters, Conformer3D, ConformerAlignmentParameters, CoordinateRmsdParameters,
    Molecule, OperationError,
};
use serde::Deserialize;
use serde_json::Value;

const NUMERICAL_TOLERANCE: f64 = 1.0e-8;

#[derive(Clone, Debug, PartialEq, Eq, serde::Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Input {
    pub case: crate::registry::SmilesCase,
    pub row: usize,
    pub preparation: Option<Value>,
}

/// Native preparation supplies the exact source conformers and call arguments,
/// as with forcefield prepared coordinates. Expected results are not inputs.
#[derive(Deserialize)]
struct Call {
    case_id: String,
    call_index: usize,
    operation: String,
    source: OracleSource,
    parameters: Value,
}

fn error_kind(error: &OperationError) -> &'static str {
    let OperationError::Alignment(error) = error else {
        return "unexpected_operation_error";
    };
    match error {
        AlignmentError::ConformerNotFound { .. } => "conformer_not_found",
        AlignmentError::WeightCountMismatch { .. } => "weight_count_mismatch",
        AlignmentError::NoSubstructureMatch => "no_substructure_match",
        // Preserve an actual unexpected failure as a mismatch, not a success
        // or a source-defined error category that it does not implement.
        _ => "unexpected_alignment_error",
    }
}

pub fn run(input: &Input) -> crate::Result<crate::registry::Record> {
    use serde_json::json;
    let prescription = input
        .preparation
        .as_ref()
        .ok_or("missing MolAlign preparation")?;
    let call: Call = serde_json::from_value(prescription.clone()).map_err(|e| e.to_string())?;
    let mut before = ObservableState {
        probe: None,
        reference: None,
        molecule: None,
    };
    let mut after = before.clone();
    let mut source_after = before.clone();
    let (status, result, kind, error_type, error_message) = if call.operation == "input_parse" {
        match Molecule::from_smiles(&input.case.smiles) {
            Ok(_) => ("ok", Some(Value::Null), None, None, None),
            Err(error) => (
                "error",
                None,
                Some("input_parse_error"),
                Some("SmilesParseError"),
                Some(error.to_string()),
            ),
        }
    } else {
        let result: Result<Value, OperationError>;
        if let (Some(probe), Some(reference)) = (&call.source.probe, &call.source.reference) {
            let probe = molecule(probe);
            let reference = molecule(reference);
            before.probe = Some(snapshot(&probe));
            before.reference = Some(snapshot(&reference));
            after = before.clone();
            result = match call.operation.as_str() {
                "alignment_transform" => probe.alignment_transform_to_with_params(&reference, &alignment_parameters(&call.parameters))
                    .map_err(OperationError::Alignment)
                    .map(|r| json!({"rmsd":r.rmsd,"transform":r.transform.matrix,"atom_map":r.atom_map.iter().map(|p|[p.probe_atom,p.reference_atom]).collect::<Vec<_>>()})),
                "best_alignment" => probe.best_alignment_to_with_params(&reference, &best_parameters(&call.parameters))
                    .map_err(OperationError::Alignment)
                    .map(|r| json!({"rmsd":r.rmsd,"transform":r.transform.matrix,"atom_map":r.atom_map.iter().map(|p|[p.probe_atom,p.reference_atom]).collect::<Vec<_>>()})),
                "coordinate_rmsd" => probe.coordinate_rmsd_to_with_params(&reference, &coordinate_parameters(&call.parameters))
                    .map_err(OperationError::Alignment)
                    .map(|r| json!({"rmsd":r})),
                "align_to" => probe.with_alignment_to_with_params(&reference, &alignment_parameters(&call.parameters))
                    .map(|(aligned,r)| {
                        after.probe=Some(snapshot(&aligned));
                        json!({"rmsd":r.rmsd,"transform":r.transform.matrix,"atom_map":r.atom_map.iter().map(|p|[p.probe_atom,p.reference_atom]).collect::<Vec<_>>()})
                    }),
                operation => return Err(format!("unknown pair MolAlign operation: {operation}")),
            };
            source_after.probe = Some(snapshot(&probe));
            source_after.reference = Some(snapshot(&reference));
            after.reference = source_after.reference.clone();
            if call.operation != "align_to" || result.is_err() {
                after.probe = source_after.probe.clone();
            }
        } else {
            let molecule = molecule(
                call.source
                    .molecule
                    .as_ref()
                    .ok_or("missing MolAlign molecule")?,
            );
            before.molecule = Some(snapshot(&molecule));
            after = before.clone();
            result = match call.operation.as_str() {
                "all_conformer_best_rms" => molecule.all_conformer_best_rmsds_with_params(&all_conformer_parameters(&call.parameters))
                    .map_err(OperationError::Alignment)
                    .map(|r| json!({"rmsds":r.iter().map(|r|r.rmsd).collect::<Vec<_>>(),"conformer_pairs":r.iter().map(|r|[r.probe_conformer_id,r.reference_conformer_id]).collect::<Vec<_>>()})),
                "align_conformers" => molecule.with_aligned_conformers_with_params(&conformer_parameters(&call.parameters))
                    .map(|(aligned,r)| {
                        after.molecule=Some(snapshot(&aligned));
                        json!({"rmsds":r.rmsds})
                    }),
                operation => return Err(format!("unknown single MolAlign operation: {operation}")),
            };
            source_after.molecule = Some(snapshot(&molecule));
            if call.operation != "align_conformers" || result.is_err() {
                after.molecule = source_after.molecule.clone();
            }
        }
        match result {
            Ok(result) => ("ok", Some(result), None, None, None),
            Err(error) => (
                "error",
                None,
                Some(error_kind(&error)),
                Some("AlignmentError"),
                Some(error.to_string()),
            ),
        }
    };
    let output = json!({"schema_version":1,"case_id":call.case_id,"call_index":call.call_index,
        "operation":call.operation,"source":prescription["source"],"parameters":call.parameters,
        "status":status,"result":result,"error_kind":kind,"error_type":error_type,"error_message":error_message,
        "before":before,"after":after,"source_after":source_after});
    Ok(crate::registry::Record {
        input: crate::registry::Input::MolAlign(input.clone()),
        output: crate::registry::Value::MolAlign(output),
    })
}

/// Preserve the former comparison: integer IDs, maps and order are exact;
/// every floating result and coordinate uses the original absolute 1e-8.
/// Error categories are exact; language-specific exception text is not parity.
pub fn matches(expected: &Value, actual: &Value) -> bool {
    match (expected, actual) {
        (Value::Object(a), Value::Object(b)) => {
            a.len() == b.len()
                && a.iter().all(|(key, value)| {
                    b.get(key).is_some_and(|other| {
                        if matches!(key.as_str(), "error_type" | "error_message") {
                            (value.is_null() && other.is_null())
                                || (value.is_string() && other.is_string())
                        } else {
                            matches(value, other)
                        }
                    })
                })
        }
        (Value::Array(a), Value::Array(b)) => {
            a.len() == b.len() && a.iter().zip(b).all(|(a, b)| matches(a, b))
        }
        (Value::Number(a), Value::Number(b)) if a.is_f64() && b.is_f64() => {
            (a.as_f64().unwrap() - b.as_f64().unwrap()).abs() <= NUMERICAL_TOLERANCE
        }
        _ => expected == actual,
    }
}

pub fn validate_reference(recipe: &Input, prepared: &Input, output: &Value) -> crate::Result<()> {
    if recipe.case != prepared.case || recipe.row != prepared.row || recipe.preparation.is_some() {
        return Err("MolAlign recipe/case mismatch".into());
    }
    let prescription = prepared
        .preparation
        .as_ref()
        .ok_or("missing MolAlign preparation")?;
    let call: Call = serde_json::from_value(prescription.clone()).map_err(|e| e.to_string())?;
    let record: OracleRecord = serde_json::from_value(output.clone()).map_err(|e| e.to_string())?;
    let operation = if call.operation == "input_parse" {
        "input_parse"
    } else {
        [
            "alignment_transform",
            "best_alignment",
            "coordinate_rmsd",
            "align_to",
            "all_conformer_best_rms",
            "align_conformers",
        ][recipe.row % 6]
    };
    if record.schema_version != 1
        || call.case_id != format!("smiles-{:05}", recipe.row + 1)
        || call.call_index != 0
        || call.operation != operation
        || ["case_id", "call_index", "operation", "source", "parameters"]
            .iter()
            .any(|key| output[*key] != prescription[*key])
        || !matches!(record.status.as_str(), "ok" | "error")
        || (record.status == "ok") != record.error_type.is_none()
        || (record.status == "ok") != record.error_message.is_none()
        || output["source_after"] != output["before"]
    {
        return Err("MolAlign reference identity/status mismatch".into());
    }
    let sources = [
        call.source.probe.as_ref(),
        call.source.reference.as_ref(),
        call.source.molecule.as_ref(),
    ];
    if call.operation == "input_parse" {
        if call.source.input_smiles.as_deref() != Some(recipe.case.smiles.as_str())
            || sources.iter().any(|s| s.is_some())
        {
            return Err("MolAlign rejected input changed".into());
        }
    } else if sources
        .iter()
        .flatten()
        .any(|s| s.smiles != recipe.case.smiles)
        || sources.iter().all(|s| s.is_none())
    {
        return Err("MolAlign prepared molecule identity mismatch".into());
    }
    Ok(())
}

pub fn validate_focused(fixture: &Value, rows: &[Value], count: usize) -> crate::Result<()> {
    let cases = fixture["cases"]
        .as_array()
        .ok_or("missing MolAlign cases")?;
    let mut ordinal = 0;
    for case in cases {
        for (index, call) in case["calls"]
            .as_array()
            .ok_or("missing MolAlign calls")?
            .iter()
            .enumerate()
        {
            let row = rows.get(ordinal).ok_or("missing MolAlign reference row")?;
            let record: OracleRecord =
                serde_json::from_value(row.clone()).map_err(|e| e.to_string())?;
            let expected_source = if case.get("probe").is_some() {
                serde_json::json!({"probe":case["probe"],"reference":case["reference"]})
            } else {
                serde_json::json!({"molecule":case["molecule"]})
            };
            if record.schema_version != 1
                || row["case_id"] != case["case_id"]
                || row["call_index"] != index
                || row["parameters"] != *call
                || row["operation"] != call["operation"]
                || row["source"] != expected_source
                || !matches!(record.status.as_str(), "ok" | "error")
                || (record.status == "ok") != record.error_type.is_none()
                || (record.status == "ok") != record.error_message.is_none()
            {
                return Err(format!(
                    "MolAlign focused reference changed at row {ordinal}"
                ));
            }
            ordinal += 1;
        }
    }
    if ordinal != count || rows.len() != count {
        return Err("MolAlign focused census mismatch".into());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    #[test]
    fn comparator_keeps_exact_maps_order_status_and_source_state() {
        let expected = json!({"rmsd":1.0,"atom_map":[[0,2],[1,3]],
            "source_after":{"coordinates":[[1.0,2.0,3.0]]},"status":"ok",
            "error_kind":null,"error_type":null,"error_message":null});
        let mut changed = expected.clone();
        changed["rmsd"] = json!(1.0 + 0.5e-8);
        assert!(matches(&expected, &changed));
        changed["rmsd"] = json!(1.0 + 2.0e-8);
        assert!(!matches(&expected, &changed));
        for (field, value) in [
            ("atom_map", json!([[1, 3], [0, 2]])),
            ("atom_map", json!([[0, 2]])),
            ("atom_map", json!([[0, 3], [1, 3]])),
            ("atom_map", json!([[0.0, 2], [1, 3]])),
            ("source_after", json!({"coordinates":[[1.0,2.0,3.001]]})),
            ("status", json!("error")),
            ("error_kind", json!("conformer_not_found")),
            ("error_message", json!("unexpected failure")),
        ] {
            let mut changed = expected.clone();
            changed[field] = value;
            assert!(!matches(&expected, &changed), "missed {field} mismatch");
        }
    }

    #[test]
    fn exception_text_is_language_specific_but_category_is_not() {
        let expected = json!({"status":"error","error_kind":"conformer_not_found",
            "error_type":"ValueError","error_message":"Bad Conformer Id"});
        let actual = json!({"status":"error","error_kind":"conformer_not_found",
            "error_type":"AlignmentError","error_message":"conformer 99 not found"});
        assert!(matches(&expected, &actual));
        let mut changed = actual;
        changed["error_kind"] = json!("no_substructure_match");
        assert!(!matches(&expected, &changed));
    }
}

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

#[derive(Debug, Clone, Deserialize, serde::Serialize)]
struct ConformerState {
    id: usize,
    is_3d: bool,
    coordinates: Vec<[f64; 3]>,
}

#[derive(Debug, Clone, Deserialize, serde::Serialize)]
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

pub fn compare_rows(rows: &[Value]) {
    assert!(!rows.is_empty(), "zero MolAlign comparisons");
    for row in rows {
        let record: OracleRecord =
            serde_json::from_value(row.clone()).expect("MolAlign reference schema");
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
