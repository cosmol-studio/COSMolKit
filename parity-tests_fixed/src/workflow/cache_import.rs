//! One-time, user-approved import of the audited recovery recipe. This never
//! changes comparison validation or derives a reference result from CK.
use super::*;

const RECIPES: [(&str, &str); 2] = [
    (
        "08562057eb9addad658dd6f6ce34b76ac324c41bc94e95b77abc7b65f2605424",
        "a86630e02b9dadae5988e04c45ffaf58f6797337cea0fd0ef368852a69a171c5",
    ),
    (
        "6700a227438eaa20a45f0cdc1937dc5de1f6e71d2a5120df0e9de19bc62585e7",
        "44c33cd3fe302488c87fe12af230ba1edeb051605f5d1bb11430efc9b35e4bfb",
    ),
];

fn old_key(key: &str) -> &str {
    match key {
        "fingerprint_layered_smiles" => "layered_fingerprint_smiles",
        "fingerprint_maccs_smiles" => "maccs_fingerprint_smiles",
        "fingerprint_morgan_count_smiles" => "morgan_count_fingerprint_smiles",
        "fingerprint_morgan_smiles" => "morgan_fingerprint_smiles",
        "fingerprint_morgan_sparse_count_smiles" => "morgan_sparse_count_fingerprint_smiles",
        "fingerprint_morgan_sparse_smiles" => "morgan_sparse_fingerprint_smiles",
        "fingerprint_pattern_smiles" => "pattern_fingerprint_smiles",
        "fingerprint_topological_smiles" => "topological_fingerprint_smiles",
        "remove_hs_smiles" => "remove_hydrogens_smiles",
        _ => key,
    }
}

fn convert_input(value: &mut Value) -> Result<()> {
    match value {
        Value::Array(items) => {
            for item in items {
                convert_input(item)?;
            }
        }
        Value::Object(fields) => {
            for (old, new) in [
                ("do_isomeric_smiles", "isomeric_smiles"),
                ("do_kekule", "kekule"),
                ("remove_hydrogens", "remove_hs"),
            ] {
                if fields.contains_key(old) && fields.contains_key(new) {
                    return Err(format!("ambiguous recovery input: {old} and {new}"));
                }
                if let Some(value) = fields.remove(old) {
                    fields.insert(new.into(), value);
                }
            }
            for value in fields.values_mut() {
                convert_input(value)?;
            }
        }
        _ => {}
    }
    Ok(())
}

fn convert_row(mut row: Value) -> Result<Value> {
    convert_input(row.get_mut("input").ok_or("recovery row missing input")?)?;
    Ok(row)
}

pub(super) fn import(
    plan: &Plan,
    spec: &Spec,
    input: &Value,
    checkout: &Path,
    folder: &Path,
) -> Result<Option<usize>> {
    if !matches!(plan.selection, Selection::Corpus(_)) || reference::uses_gemmi(spec) {
        return Ok(None);
    }
    let Selection::Corpus(name) = &plan.selection else {
        unreachable!()
    };
    let source = checkout
        .join("parity-tests_fixed/expected/corpus")
        .join(name)
        .join(old_key(spec.key()));
    if !source.join("manifest.json").exists() {
        return Ok(None);
    }
    let original_manifest = read(&source.join("manifest.json"))?;
    let manifest: Manifest =
        serde_json::from_slice(&original_manifest).map_err(|e| e.to_string())?;
    let target_digest = reference::source_digest(spec)?;
    if !RECIPES.contains(&(manifest.generator_sha256.as_str(), target_digest.as_str())) {
        eprintln!("  recovery not reused: recipe pair has not been audited");
        return Ok(None);
    }
    if reference::source_digest_at(spec, checkout)? != manifest.generator_sha256 {
        return Err(format!(
            "{}: recovery generator does not match its manifest",
            spec.key()
        ));
    }
    let current = identity_with_output_digest(&plan.selection, spec, input, String::new(), 0)?;
    if manifest.schema != current.schema
        || manifest.selection != current.selection
        || manifest.task != old_key(spec.key())
        || manifest.reference_identity != current.reference_identity
        || manifest.platform != current.platform
    {
        return Err(format!(
            "{}: recovery reference identity mismatch",
            spec.key()
        ));
    }
    let mut stored: Value =
        serde_json::from_slice(&read(&source.join("input.json"))?).map_err(|e| e.to_string())?;
    if encoded_digest(&stored)? != manifest.input_sha256 {
        return Err(format!("{}: corrupt recovery input", spec.key()));
    }
    convert_input(&mut stored)?;
    if stored != *input {
        eprintln!("  recovery not reused: actual parameters differ after approved key conversion");
        return Ok(None);
    }
    drop(stored);
    let rows = OwnedRows::copy_from(&source.join("reference.jsonl"))?;
    if rows.sha256 != manifest.output_sha256 || rows.len() != manifest.rows {
        return Err(format!("{}: corrupt recovery output", spec.key()));
    }
    let temporary = tempfile::Builder::new()
        .prefix(".import-")
        .tempdir_in(folder)
        .map_err(|e| e.to_string())?;
    let mut writer = RowWriter::new_in(temporary.path())?;
    for row in rows.iter()? {
        writer.push(convert_row(row?)?)?;
    }
    let converted = writer.finish()?;
    validate_rows(spec, input, converted.iter()?, converted.len())?;
    let current = identity_with_output_digest(
        &plan.selection,
        spec,
        input,
        converted.sha256.clone(),
        converted.len(),
    )?;
    if current.generator_sha256 != target_digest {
        return Err("reference generator changed during import".into());
    }
    let count = converted.len();
    converted.persist(&temporary.path().join("reference.jsonl"))?;
    let mut input_file = BufWriter::new(
        fs::File::create(temporary.path().join("input.json")).map_err(|e| e.to_string())?,
    );
    serde_json::to_writer(&mut input_file, input).map_err(|e| e.to_string())?;
    input_file.flush().map_err(|e| e.to_string())?;
    drop(input_file);
    fs::write(temporary.path().join("manifest.json"), encode(&current)?)
        .map_err(|e| e.to_string())?;
    fs::write(
        temporary.path().join("import.json"),
        encode(&json!({
            "source": source, "original_manifest": manifest,
            "original_manifest_sha256": digest(&original_manifest),
            "conversion": "input-only parameter key renames; output values unchanged",
        "audited_recipe_pair": [manifest.generator_sha256, target_digest],
        }))?,
    )
    .map_err(|e| e.to_string())?;
    let destination = folder.join(spec.key());
    let backup = folder.join(format!(
        ".invalid-{}-{}",
        spec.key(),
        temporary.path().file_name().unwrap().to_string_lossy()
    ));
    let exists = destination.exists();
    if exists {
        fs::rename(&destination, &backup).map_err(|e| e.to_string())?;
    }
    if let Err(error) = fs::rename(temporary.path(), &destination) {
        if exists {
            fs::rename(backup, destination)
                .map_err(|e| format!("publish {error}; rollback {e}"))?;
        }
        return Err(error.to_string());
    }
    Ok(Some(count))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn input_conversion_preserves_values_and_rejects_collisions() {
        let mut value = json!([{"do_isomeric_smiles":true, "do_kekule":false,
            "nested":{"remove_hydrogens":true}, "text":"remove_hydrogens"}]);
        convert_input(&mut value).unwrap();
        assert_eq!(
            value,
            json!([{"isomeric_smiles":true, "kekule":false,
            "nested":{"remove_hs":true}, "text":"remove_hydrogens"}])
        );
        assert!(convert_input(&mut json!({"remove_hydrogens":false,"remove_hs":true})).is_err());
    }
    #[test]
    fn native_output_is_never_converted_or_rounded() {
        let row = json!({"input":{"remove_hydrogens":true}, "output":{
            "remove_hydrogens":false,"integer":u64::MAX,
            "floats":[-0.0,0.9389657496748851_f64]}});
        let output = row["output"].clone();
        let converted = convert_row(row).unwrap();
        assert_eq!(
            encode(&converted["output"]).unwrap(),
            encode(&output).unwrap()
        );
        assert_eq!(converted["input"], json!({"remove_hs":true}));
        let parsed: Value = serde_json::from_slice(&encode(&converted).unwrap()).unwrap();
        assert_eq!(
            parsed["output"]["floats"][0].as_f64().unwrap().to_bits(),
            (-0.0_f64).to_bits()
        );
        assert_eq!(
            parsed["output"]["floats"][1].as_f64().unwrap().to_bits(),
            0.9389657496748851_f64.to_bits()
        );
        assert_eq!(parsed["output"]["integer"].as_u64(), Some(u64::MAX));
    }
    #[test]
    fn renamed_tasks_do_not_waive_different_parameters() {
        assert_eq!(old_key("remove_hs_smiles"), "remove_hydrogens_smiles");
        let mut old = json!({"remove_hydrogens":true});
        convert_input(&mut old).unwrap();
        assert_ne!(old, json!({"remove_hs":false}));
        assert!(convert_row(json!({"output":null})).is_err());
    }
}
