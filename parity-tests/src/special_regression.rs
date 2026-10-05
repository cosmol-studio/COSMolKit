//! Fixed, source-authenticated regressions with a separate preparation stage.
//! This module only handles data identity; it never calls chemistry algorithms.
use crate::{Preparation, Result, digest, encode, read, registry, root};
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::{collections::BTreeSet, fs, path::Path, process::Command};

#[derive(Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Manifest {
    schema: u32,
    category: String,
    task: String,
    rdkit_version: String,
    reference_pin_sha256: String,
    generator_sha256: String,
    registry_sha256: String,
    platform: String,
    input_sha256: String,
    reference_sha256: String,
    rows: usize,
}

/// Owned, checked snapshots: comparison cannot reread changed expectations.
pub struct Snapshot {
    pub fixture: Value,
    pub rows: Vec<Value>,
}

fn task(key: &str) -> Result<&'static registry::SpecialRegression> {
    registry::SPECIAL_REGRESSIONS
        .iter()
        .find(|task| task.key == key)
        .ok_or_else(|| format!("unknown special regression: {key}"))
}

fn identity(
    task: &registry::SpecialRegression,
    input: &[u8],
    reference: &[u8],
) -> Result<Manifest> {
    Ok(Manifest {
        schema: 1,
        category: "special_regression".into(),
        task: task.key.into(),
        rdkit_version: registry::RDKIT_VERSION.into(),
        reference_pin_sha256: digest(crate::REFERENCE_PIN.as_bytes()),
        generator_sha256: digest(&read(&root().join(task.generator))?),
        registry_sha256: digest(include_str!("registry.rs").as_bytes()),
        platform: format!("{}-{}", std::env::consts::ARCH, std::env::consts::OS),
        input_sha256: digest(input),
        reference_sha256: digest(reference),
        rows: task.rows,
    })
}

fn validate(input: &[u8], reference: &[u8], count: usize) -> Result<Snapshot> {
    let fixture: Value = serde_json::from_slice(input).map_err(|e| e.to_string())?;
    // These are distinct pinned spellings, not a version-normalization rule:
    // rdBase reports 2026.03.1; Python distribution/fixture reports 2026.3.1.
    let pin: Value = serde_json::from_str(crate::REFERENCE_PIN).map_err(|e| e.to_string())?;
    if fixture["schema_version"] != 1
        || fixture["reference"]["version"] != pin["python_distribution_version"]
    {
        return Err("special regression fixture schema/reference mismatch".into());
    }
    let mut ids = Vec::new();
    for field in ["cases", "octahedral_switch_cases"] {
        let cases = fixture[field]
            .as_array()
            .ok_or_else(|| format!("missing fixture table: {field}"))?;
        for case in cases {
            ids.push(
                case["case_id"]
                    .as_str()
                    .ok_or("missing fixture case id")?
                    .to_owned(),
            );
        }
    }
    if ids.len() != count || ids.iter().collect::<BTreeSet<_>>().len() != count {
        return Err("special regression case census/uniqueness mismatch".into());
    }
    let reference = std::str::from_utf8(reference).map_err(|e| e.to_string())?;
    let rows: Vec<Value> = reference
        .lines()
        .map(|line| serde_json::from_str(line).map_err(|e| e.to_string()))
        .collect::<Result<_>>()?;
    if rows.len() != count {
        return Err("special regression reference row count mismatch".into());
    }
    for (id, row) in ids.iter().zip(&rows) {
        if row["case_id"].as_str() != Some(id)
            || !matches!(row["status"].as_str(), Some("ok" | "error"))
            || !row["before"].is_object()
            || !row["after"].is_object()
            || !row["environment"].is_object()
            || !row["conf_id"].is_i64()
            || !row["replace_existing_tags"].is_boolean()
        {
            return Err(format!(
                "special regression reference identity/schema: {id}"
            ));
        }
    }
    Ok(Snapshot { fixture, rows })
}

/// Read-only global preflight of the complete fixed selection; no CK calls.
pub fn preflight(key: &str, data: &Path) -> Result<Snapshot> {
    let task = task(key)?;
    let input = read(&root().join(task.fixture))?;
    let dir = data.join(format!("special-regression-{}", task.key));
    let load = || -> Result<Snapshot> {
        let published_input = read(&dir.join("input.json"))?;
        let reference = read(&dir.join("reference.jsonl"))?;
        let manifest: Manifest = serde_json::from_slice(&read(&dir.join("manifest.json"))?)
            .map_err(|e| e.to_string())?;
        if published_input != input || manifest != identity(task, &input, &reference)? {
            return Err("stale/corrupt special regression identity or checksums".into());
        }
        validate(&input, &reference, task.rows)
    };
    load().map_err(|error| {
        format!(
            "{error}\n0 Rust operation calls; prepare first: cargo run -p cosmolkit-parity-tests --release -- prepare --special-regression {key} --data {}",
            data.display()
        )
    })
}

/// Complete selected special lane, checked before any comparison begins.
pub fn preflight_selection(keys: &[String], data: &Path) -> Result<usize> {
    let mut seen = BTreeSet::new();
    let mut rows = 0;
    for key in keys {
        if !seen.insert(key) {
            return Err(format!("duplicate special regression: {key}"));
        }
        rows += preflight(key, data)?.rows.len();
    }
    Ok(rows)
}

/// Prepare with the existing pinned generator; tests never call this function.
pub fn prepare(key: &str, data: &Path, python: &Path) -> Result<Preparation> {
    let task = task(key)?;
    fs::create_dir_all(data).map_err(|e| e.to_string())?;
    let lock = fs::OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .truncate(false)
        .open(data.join(".prepare.lock"))
        .map_err(|e| e.to_string())?;
    lock.lock().map_err(|e| e.to_string())?;
    if preflight(key, data).is_ok() {
        return Ok(Preparation {
            reused_tasks: 1,
            generated_tasks: 0,
            rows: task.rows,
        });
    }
    let input = read(&root().join(task.fixture))?;
    let generator = read(&root().join(task.generator))?;
    let temporary = tempfile::Builder::new()
        .prefix(".prepare-special-")
        .tempdir_in(data)
        .map_err(|e| e.to_string())?;
    // The legacy generator selects this fixed matrix by the exact output name,
    // without reading/expanding its SMILES option. Retain that owner unchanged.
    let output_path = temporary.path().join(task.output);
    let output_path = fs::canonicalize(temporary.path())
        .map_err(|e| e.to_string())?
        .join(output_path.file_name().unwrap());
    let output = Command::new(python)
        .arg(root().join(task.generator))
        .arg("--output")
        .arg(&output_path)
        .current_dir(root())
        .output()
        .map_err(|e| e.to_string())?;
    fs::write(temporary.path().join("generator.stdout"), &output.stdout)
        .map_err(|e| e.to_string())?;
    fs::write(temporary.path().join("generator.stderr"), &output.stderr)
        .map_err(|e| e.to_string())?;
    if !output.status.success() {
        let evidence = temporary.keep();
        return Err(format!(
            "special regression generator failed: {}; evidence {}\n{}",
            output.status,
            evidence.display(),
            String::from_utf8_lossy(&output.stderr)
        ));
    }
    let reference = read(&output_path)?;
    let checked: Result<()> = (|| {
        if input != read(&root().join(task.fixture))?
            || generator != read(&root().join(task.generator))?
        {
            return Err("generator/fixture changed during preparation".into());
        }
        validate(&input, &reference, task.rows)?;
        fs::write(temporary.path().join("input.json"), &input).map_err(|e| e.to_string())?;
        fs::write(temporary.path().join("reference.jsonl"), &reference)
            .map_err(|e| e.to_string())?;
        fs::write(
            temporary.path().join("manifest.json"),
            encode(&identity(task, &input, &reference)?)?,
        )
        .map_err(|e| e.to_string())?;
        Ok(())
    })();
    if let Err(error) = checked {
        return Err(format!("{error}; evidence {}", temporary.keep().display()));
    }
    let destination = data.join(format!("special-regression-{}", task.key));
    let backup = destination.exists().then(|| {
        data.join(format!(
            ".invalid-special-{}-{}",
            task.key,
            temporary.path().file_name().unwrap().to_string_lossy()
        ))
    });
    if let Some(backup) = &backup {
        fs::rename(&destination, backup).map_err(|e| e.to_string())?;
    }
    if let Err(error) = fs::rename(temporary.path(), &destination) {
        if let Some(backup) = backup {
            fs::rename(backup, destination).map_err(|e| e.to_string())?;
        }
        return Err(error.to_string());
    }
    preflight(key, data)?;
    Ok(Preparation {
        reused_tasks: 0,
        generated_tasks: 1,
        rows: task.rows,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixed_pair() -> (Vec<u8>, Vec<u8>) {
        let input = serde_json::json!({
            "schema_version": 1,
            "reference": {"version": "2026.3.1"},
            "cases": [{"case_id": "fixed_ok"}],
            "octahedral_switch_cases": [{"case_id": "fixed_error"}],
        });
        let rows = ["fixed_ok", "fixed_error"].map(|id| {
            serde_json::json!({
                "case_id": id, "status": if id == "fixed_ok" {"ok"} else {"error"},
                "before": {}, "after": {}, "environment": {},
                "conf_id": -1, "replace_existing_tags": true,
            })
            .to_string()
        });
        (
            serde_json::to_vec(&input).unwrap(),
            rows.join("\n").into_bytes(),
        )
    }

    #[test]
    fn special_regression_checks_complete_census_order_and_errors() {
        let (input, reference) = fixed_pair();
        let ready = validate(&input, &reference, 2).unwrap();
        assert_eq!(ready.rows.len(), 2);
        assert_eq!(ready.rows[1]["status"], "error");
        assert!(validate(&input, &reference, 3).is_err());
        let reversed = std::str::from_utf8(&reference)
            .unwrap()
            .lines()
            .rev()
            .collect::<Vec<_>>()
            .join("\n");
        assert!(validate(&input, reversed.as_bytes(), 2).is_err());
        let mut invalid = ready.rows;
        invalid[0]["case_id"] = Value::String("fixed_error".into());
        let duplicate = invalid
            .iter()
            .map(Value::to_string)
            .collect::<Vec<_>>()
            .join("\n");
        assert!(validate(&input, duplicate.as_bytes(), 2).is_err());
    }

    #[test]
    fn special_regression_missing_data_is_read_only_and_actionable() {
        let dir = tempfile::tempdir().unwrap();
        let error = preflight("assign_chiral_tags_from_structure", dir.path())
            .err()
            .unwrap();
        assert!(error.contains("0 Rust operation calls"));
        assert!(error.contains("--special-regression assign_chiral_tags_from_structure"));
        assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 0);
        assert!(task("unknown").is_err());
    }

    #[test]
    fn special_regression_identity_binds_inputs_outputs_and_generator() {
        let task = task("assign_chiral_tags_from_structure").unwrap();
        let (input, reference) = fixed_pair();
        let baseline = identity(task, &input, &reference).unwrap();
        assert_ne!(baseline, identity(task, b"changed", &reference).unwrap());
        assert_ne!(baseline, identity(task, &input, b"changed").unwrap());
        assert_eq!(baseline.rows, 77);
        assert_eq!(baseline.category, "special_regression");
        assert_eq!(baseline.generator_sha256.len(), 64);
    }
}
