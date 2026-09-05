use std::fs;
use std::path::Path;
use std::process::Command;

use serde_json::Value;

const SCHEMA_VERSION: &str = "cosmolkit-inchi-official-c-v1";
const UPSTREAM_TAG: &str = "v1.07.5";
const UPSTREAM_COMMIT: &str = "11a87982bb518f57ac013f0b258c283655e1ea1d";
const SOURCE_TREE_SHA256: &str = "4f903503c370a6dd21de0a2e6beb77e076a1cd6063b34b6bde4b57d91f9c4edd";
const API_VERSION: &str = "1.07.5";

fn repository_root() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(Path::parent)
        .expect("cosmolkit-inchi must be located under crates/")
}

fn read_json(path: &Path) -> Value {
    let bytes = fs::read(path).unwrap_or_else(|error| panic!("{}: {error}", path.display()));
    serde_json::from_slice(&bytes)
        .unwrap_or_else(|error| panic!("{} is not valid JSON: {error}", path.display()))
}

fn validate_reference_metadata() -> Vec<u8> {
    let root = repository_root();
    let schema_path = root.join("testdata/inchi/reference/official_c_oracle.schema.json");
    let golden_path = root.join("testdata/inchi/reference/official_c_oracle_version.jsonl");

    let schema = read_json(&schema_path);
    assert_eq!(
        schema["properties"]["schema_version"]["const"],
        SCHEMA_VERSION
    );
    assert_eq!(
        schema["properties"]["oracle"]["properties"]["tag"]["const"],
        UPSTREAM_TAG
    );
    assert_eq!(
        schema["properties"]["oracle"]["properties"]["commit"]["const"],
        UPSTREAM_COMMIT
    );
    assert_eq!(
        schema["properties"]["oracle"]["properties"]["tree_sha256"]["const"],
        SOURCE_TREE_SHA256
    );
    assert_eq!(
        schema["properties"]["oracle"]["properties"]["api_version"]["const"],
        API_VERSION
    );

    let golden_bytes =
        fs::read(&golden_path).unwrap_or_else(|error| panic!("{}: {error}", golden_path.display()));
    let mut golden_lines = golden_bytes.split(|byte| *byte == b'\n');
    let golden_line = golden_lines.next().expect("version golden is empty");
    assert!(
        golden_lines.all(|line| line.is_empty()),
        "version golden must contain exactly one JSONL record"
    );
    let golden: Value = serde_json::from_slice(golden_line)
        .unwrap_or_else(|error| panic!("{} is not valid JSONL: {error}", golden_path.display()));
    assert_eq!(golden["schema_version"], SCHEMA_VERSION);
    assert_eq!(golden["oracle"]["tag"], UPSTREAM_TAG);
    assert_eq!(golden["oracle"]["commit"], UPSTREAM_COMMIT);
    assert_eq!(golden["oracle"]["tree_sha256"], SOURCE_TREE_SHA256);
    assert_eq!(golden["oracle"]["api_version"], API_VERSION);
    assert_eq!(golden["case_id"], "oracle-version");
    assert_eq!(golden["operation"], "version");
    assert_eq!(golden["output"]["status"], 0);
    assert_eq!(golden["output"]["value"], API_VERSION);

    golden_bytes
}

#[test]
fn official_c_oracle_reference_metadata_is_pinned() {
    let _ = validate_reference_metadata();
}

#[test]
#[ignore = "requires the pinned vendored official InChI source; run explicitly with --ignored"]
fn official_c_oracle_version_and_provenance() {
    let root = repository_root();
    let runner_path = root.join("tools/oracles/official_inchi/run.sh");
    let golden_bytes = validate_reference_metadata();

    let output = Command::new("sh")
        .arg(&runner_path)
        .arg("--version-record")
        .current_dir(root)
        .output()
        .unwrap_or_else(|error| panic!("failed to start {}: {error}", runner_path.display()));
    assert!(
        output.status.success(),
        "official C oracle failed with {}:\n{}",
        output.status,
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(output.stdout, golden_bytes);
}
