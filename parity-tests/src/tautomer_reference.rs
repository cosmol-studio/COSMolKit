//! Import freshly authenticated native TAU observations during preparation only.
//! Streaming source reads retain all rows; no CK operation generates expectations.
use crate::molecular::{Outcome, Stage};
use crate::registry::{
    self, Corpus, Input, Operation, Record, Task, Value,
    molecule_plan::{Profile, TaskId, TautomerProfile},
};
use crate::{Result, root};
use serde_json::{Value as Json, json};
use sha2::{Digest, Sha256};
use std::io::{BufRead, BufReader, Read};
const INPUT: &str = "testdata/smiles/corpus/smiles_5000.smi";
const INPUT_SHA: &str = "a4d579cd72621af27772256bb23ba796452276bb924fd20aac83625ffa67d849";
const DIRECTORY: &str = "testdata/tautomer/expected/rdkit/smiles_5000";
const ENUM_SHA: &str = "4a1aee883d98122c21e838d2f1bcf3a4ec4d5c65cdb62819b62f8195aa8d9d82";
const CANON_SHA: &str = "9d32b78df71a44c09a9d301606d8eac449d7b64f11c5904cd9a73fdb89887a08";
const MANIFEST_SHA: &str = "cd96f380b449cd057f42461770c14dec96c99f9f5f1eb79f1b69a1b12e78ce48";
const SOURCE_PIN: &str = "351f8f378f8ad6bbd517980c38896e66bf907af8";
pub fn handles(task: &Task) -> bool {
    matches!(
        task.operation,
        Operation::Molecular(TaskId::TautomerEnumeration | TaskId::TautomerCanonicalization)
    )
}
fn source(task: &Task) -> (&'static str, &'static str) {
    match task.operation {
        Operation::Molecular(TaskId::TautomerEnumeration) => ("tautomer.jsonl", ENUM_SHA),
        Operation::Molecular(TaskId::TautomerCanonicalization) => {
            ("canonicalization.jsonl", CANON_SHA)
        }
        _ => unreachable!("selected TAU task"),
    }
}
pub fn provenance(task: &Task) -> Option<Json> {
    handles(task).then(||{let (file,sha)=source(task);json!({"kind":"authenticated_native_reference_import","path":format!("{DIRECTORY}/{file}"),"sha256":sha,"input_sha256":INPUT_SHA,"source_revision":SOURCE_PIN,"reference_version":registry::RDKIT_VERSION,"records":5000,"branches":10000,"source_identity_manifest":format!("{DIRECTORY}/manifest.json"),"source_identity_manifest_sha256":MANIFEST_SHA})})
}
fn verified(path: &std::path::Path, expected: &str) -> Result<()> {
    let mut input = std::fs::File::open(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let mut hasher = Sha256::new();
    let mut buffer = vec![0u8; 1024 * 1024];
    loop {
        let count = input.read(&mut buffer).map_err(|e| e.to_string())?;
        if count == 0 {
            break;
        }
        hasher.update(&buffer[..count]);
    }
    if hasher
        .finalize()
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect::<String>()
        != expected
    {
        return Err(format!(
            "{}: authenticated reference checksum mismatch",
            path.display()
        ));
    }
    Ok(())
}
fn observation(profile: Profile, branch: &Json) -> Result<Outcome> {
    if branch["ok"] != true {
        return Ok(Outcome::Error {
            stage: Stage::Operation,
            detail: serde_json::to_string(&branch["error"]).map_err(|e| e.to_string())?,
        });
    }
    if !branch["error"].is_null() {
        return Err("successful TAU reference contains error".into());
    }
    let mut observed = branch
        .as_object()
        .ok_or("TAU branch is not an object")?
        .clone();
    observed.remove("parameters");
    observed.remove("ok");
    observed.remove("error");
    match profile {
        Profile::TautomerEnumeration { .. } => serde_json::from_value(Json::Object(observed))
            .map(Outcome::TautomerEnumeration)
            .map_err(|e| e.to_string()),
        Profile::TautomerCanonicalization { .. } => serde_json::from_value(Json::Object(observed))
            .map(Outcome::TautomerCanonicalization)
            .map_err(|e| e.to_string()),
        _ => Err("profile outside TAU reference import".into()),
    }
}
pub fn generate(task: &Task, cases: &Corpus) -> Result<Vec<Record>> {
    if !handles(task) {
        return Err("task outside TAU reference import".into());
    }
    verified(&root().join(INPUT), INPUT_SHA)?;
    let original = crate::molecular::read_corpus(&root().join(INPUT))?;
    if cases.molecules != original || original.len() != 5000 {
        return Err("TAU reference requires its exact ordered 5000 input".into());
    }
    verified(&root().join(DIRECTORY).join("manifest.json"), MANIFEST_SHA)?;
    let manifest: Json =
        serde_json::from_slice(&crate::read(&root().join(DIRECTORY).join("manifest.json"))?)
            .map_err(|e| e.to_string())?;
    if manifest["source_revision"] != SOURCE_PIN
        || manifest["input_sha256"] != INPUT_SHA
        || manifest["enumeration_sha256"] != ENUM_SHA
        || manifest["canonicalization_sha256"] != CANON_SHA
        || manifest["rows"] != 5000
    {
        return Err("TAU source identity manifest mismatch".into());
    }
    // Its full native-build proof is retained in the manifest; its content digest
    // is also covered by the preparation identity below, without mutating goldens.
    let (name, hash) = source(task);
    let path = root().join(DIRECTORY).join(name);
    verified(&path, hash)?;
    let Operation::Molecular(id) = task.operation else {
        unreachable!()
    };
    let profiles = id.profiles();
    let mut records = Vec::with_capacity(10000);
    let mut rows = 0;
    let reader = BufReader::new(std::fs::File::open(path).map_err(|e| e.to_string())?);
    for line in reader.lines() {
        let row: Json =
            serde_json::from_str(&line.map_err(|e| e.to_string())?).map_err(|e| e.to_string())?;
        let case = original
            .get(rows)
            .ok_or("TAU reference contains extra row")?;
        if row["row"] != rows
            || row["case_id"] != format!("smiles_5000:{rows}")
            || row["smiles"] != case.smiles
            || row["schema_version"] != 1
            || row["sanitize"] != true
            || row["remove_hs"] != true
        {
            return Err(format!(
                "TAU reference input identity mismatch at row {rows}"
            ));
        }
        for &profile in &profiles {
            let parameters = match profile {
                Profile::TautomerEnumeration { parameters }
                | Profile::TautomerCanonicalization { parameters } => parameters,
                _ => unreachable!(),
            };
            let output = if row["parse"]["ok"] != true {
                Outcome::Error {
                    stage: Stage::Parse,
                    detail: serde_json::to_string(&row["parse"]["error"])
                        .map_err(|e| e.to_string())?,
                }
            } else {
                let branch = &row["branches"][parameters.branch()];
                let mut raw = branch["parameters"]
                    .as_object()
                    .ok_or("missing TAU branch parameters")?
                    .clone();
                if raw.remove("name") != Some(Json::String(parameters.branch().into())) {
                    return Err("TAU branch name mismatch".into());
                }
                let native: TautomerProfile =
                    serde_json::from_value(Json::Object(raw)).map_err(|e| e.to_string())?;
                if native != parameters {
                    return Err("TAU branch frozen parameters differ".into());
                }
                observation(profile, branch)?
            };
            records.push(Record {
                input: Input::Molecular {
                    case: case.clone(),
                    profile,
                },
                output: Value::Molecular(output),
            });
        }
        rows += 1;
    }
    if rows != 5000 || records.len() != 10000 {
        return Err("TAU reference must retain 5000 rows/10000 branches".into());
    }
    Ok(records)
}
