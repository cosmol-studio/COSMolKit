//! Explicit import of the complete original coordinate generation. Original
//! bytes and manifest remain provenance; no newly authenticated oracle is claimed.
use crate::{
    Corpus, Input, LabeledRecord, Record, Result, Task, digest, encode, read, reference_label,
    registry,
};
use serde_json::{Value, json};
use std::path::PathBuf;

const GENERATION: &str =
    "coordinates_2d_smiles-1a914c53f1834eaf2e6cb077935895b3a53d4e6671573419706febd904aaa589";
const INPUT_SHA: &str = "4155b3b1381d25c3525b0152fcd09393fb310c925246295ef28858d8d9ce4d42";
const REFERENCE_SHA: &str = "77f0989cb73406dfb8411b28d52ce9de3112929bdfa96eb43accd6b60c17f7e5";
const MANIFEST_SHA: &str = "0e0e559e62eb772531b3318a9edba6c22bd0c5e8818ed61941e34131610abe1a";

pub(crate) fn handles(task: &Task) -> bool {
    task.key() == "coordinates_2d_smiles" && std::env::var_os("DRAW_ORIGINAL_REFERENCE").is_some()
}

pub(crate) fn provenance(task: &Task) -> Option<Value> {
    handles(task).then(|| json!({
        "kind": "existing_reference_import", "generation": GENERATION,
        "input_sha256": INPUT_SHA, "reference_sha256": REFERENCE_SHA,
        "original_manifest_sha256": MANIFEST_SHA, "records": 5000,
        "original_manifest": {
            "schema": 3, "task": "coordinates_2d_smiles", "corpus_type": "smiles",
            "generator": "generate_coordinates_2d", "rdkit_version": "2026.03.1",
            "reference_platform": "x86_64-linux",
            "registry_sha256": "60393bea9907629c76cd4133b70b77808947ce4fe7caefb9d93ed3b64e1ad50f",
            "oracle_sha256": "c50ed1125e437eb8f1b2936b25d84c68e24149aa89353ea5ece28eeddfabd765",
            "input_sha256": INPUT_SHA, "reference_sha256": REFERENCE_SHA, "rows": 5000
        },
        "source_identity_limit": "Original manifest records reference version and producer digests, without native source/image identity. Preserved as existing expected data, never promoted to a fresh source-authenticated oracle."
    }))
}

fn verified(path: &std::path::Path, expected: &str) -> Result<Vec<u8>> {
    let bytes = read(path)?;
    if digest(&bytes) != expected {
        return Err(format!(
            "{}: original reference checksum mismatch",
            path.display()
        ));
    }
    Ok(bytes)
}

pub(crate) fn generate(task: &Task, cases: &Corpus) -> Result<Vec<Record>> {
    if !handles(task) {
        return Err("task is outside explicit DRAW original reference import".into());
    }
    let directory = PathBuf::from(
        std::env::var_os("DRAW_ORIGINAL_REFERENCE")
            .ok_or("missing original reference directory")?,
    )
    .join(GENERATION);
    let input = verified(&directory.join("input.json"), INPUT_SHA)?;
    let manifest = verified(&directory.join("manifest.json"), MANIFEST_SHA)?;
    let reference = verified(&directory.join("reference.json"), REFERENCE_SHA)?;
    let expected_inputs = registry::expand(cases, task);
    // Exact serialization preserves original order, every SMILES, case ID and profile.
    if expected_inputs.len() != 5000 || encode(&expected_inputs)? != input {
        return Err(
            "existing coordinate reference requires its exact ordered 5000 input and profile"
                .into(),
        );
    }
    let original_manifest: Value = serde_json::from_slice(&manifest).map_err(|e| e.to_string())?;
    let provenance = provenance(task).ok_or("missing import provenance")?;
    if original_manifest != provenance["original_manifest"] {
        return Err("original coordinate manifest identity mismatch".into());
    }
    let records: Vec<LabeledRecord> =
        serde_json::from_slice(&reference).map_err(|e| e.to_string())?;
    if records.len() != expected_inputs.len() {
        return Err("original coordinate reference row count mismatch".into());
    }
    for (input, row) in expected_inputs.iter().zip(&records) {
        if row.label != reference_label(task, input)? {
            return Err("original coordinate label mismatch".into());
        }
        task.validate_reference(input, &row.record.input, &row.record.output)?;
    }
    Ok(records.into_iter().map(|row| row.record).collect())
}
