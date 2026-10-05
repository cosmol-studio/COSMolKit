//! Source native DRAW observations imported without expected/checker transformations.
use crate::molecular::Outcome;
use crate::registry::{Input, Value};
use crate::{Corpus, Record, Result, Task, digest, read, registry};
use serde::Deserialize;
use serde_json::{Value as JsonValue, json};
use std::path::PathBuf;
const REFERENCE_SHA: &str = "7656e934c25c7bd1bb3d8fa729f66650c6654ec6d978442806707b2e5f54dad1";
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct NativeRow {
    row: usize,
    smiles: String,
    stage: String,
    svg: String,
    error: String,
}
pub(crate) fn handles(task: &Task) -> bool {
    task.key() == "svg_smiles" && std::env::var_os("DRAW_NATIVE_SVG_REFERENCE").is_some()
}
pub(crate) fn provenance(task: &Task) -> Option<JsonValue> {
    handles(task).then(||json!({"kind":"source_native_observation_import","source_commit":"351f8f378f8ad6bbd517980c38896e66bf907af8","native_immutable_receipt_sha256":"56832144c7353182a7dfb2e909fb1c5d000dda69bbadeac58975b276716f9157","actual_native_run":"native-corpus-svg5000-source-actual","no_freetype":true,"reference_sha256":REFERENCE_SHA,"records":5000,"expected_projection":"none; pinned native RDKit namespace preserved"}))
}
pub(crate) fn generate(task: &Task, cases: &Corpus) -> Result<Vec<Record>> {
    let path = PathBuf::from(
        std::env::var_os("DRAW_NATIVE_SVG_REFERENCE").ok_or("missing native SVG reference")?,
    );
    let bytes = read(&path)?;
    if digest(&bytes) != REFERENCE_SHA {
        return Err("native SVG reference checksum mismatch".into());
    }
    let rows: Vec<NativeRow> = std::str::from_utf8(&bytes)
        .map_err(|e| e.to_string())?
        .lines()
        .map(|s| serde_json::from_str(s).map_err(|e| e.to_string()))
        .collect::<Result<_>>()?;
    let inputs = registry::expand(cases, task);
    if !handles(task) || rows.len() != 5000 || inputs.len() != rows.len() {
        return Err("native SVG import requires its complete 5000 inputs".into());
    }
    inputs
        .into_iter()
        .zip(rows)
        .enumerate()
        .map(|(index, (input, row))| {
            let Input::Molecular { case, .. } = &input else {
                return Err("native SVG requires molecular input".into());
            };
            if row.row != index + 1
                || row.smiles != case.smiles
                || row.stage != "operation"
                || !row.error.is_empty()
            {
                return Err(format!(
                    "native SVG source/input/status mismatch at {}",
                    index + 1
                ));
            }
            let output = Value::Molecular(Outcome::Text(row.svg));
            task.validate_reference(&input, &input, &output)?;
            Ok(Record { input, output })
        })
        .collect()
}
