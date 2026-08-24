use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex, OnceLock};

use serde::Deserialize;
use sha2::{Digest, Sha256};

#[derive(Debug, Deserialize)]
struct ExpectedManifest {
    schema_version: u32,
    family: String,
    domain: String,
    profile: String,
    reference_implementation: ReferenceImplementation,
    reference_runtime: ReferenceRuntime,
    input: ManifestFile,
    outputs: Vec<ManifestOutput>,
}

#[derive(Debug, Deserialize)]
struct ReferenceImplementation {
    name: String,
    version: String,
}

#[derive(Debug, Deserialize)]
struct ReferenceRuntime {
    implementation: String,
    version: String,
}

#[derive(Debug, Deserialize)]
struct ManifestFile {
    path: String,
    sha256: String,
}

#[derive(Debug, Deserialize)]
struct ManifestOutput {
    path: String,
    output_schema_version: u32,
    generator: Vec<ManifestFile>,
    inputs: Vec<ManifestFile>,
    options: ManifestOptions,
    platform: ManifestPlatform,
    sha256: String,
    records: usize,
}

#[derive(Debug, Deserialize)]
struct ManifestOptions {
    profile: String,
    arguments: Vec<String>,
}

#[derive(Debug, Deserialize)]
struct ManifestPlatform {
    system: String,
    machine: String,
}

#[derive(Debug, Deserialize)]
struct ReferenceIdentity {
    name: String,
    version: String,
    runtime: ReferenceRuntime,
}

type ValidationResult = Result<(), String>;
type ValidationCell = Arc<OnceLock<ValidationResult>>;

static OUTPUT_VALIDATIONS: OnceLock<Mutex<HashMap<PathBuf, ValidationCell>>> = OnceLock::new();
static FILE_CHECKSUMS: OnceLock<Mutex<HashMap<PathBuf, Arc<OnceLock<Result<String, String>>>>>> =
    OnceLock::new();
#[cfg(test)]
static OUTPUT_SCAN_COUNTS: OnceLock<Mutex<HashMap<PathBuf, usize>>> = OnceLock::new();

pub fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..")
}

pub fn profile_name() -> String {
    match std::env::var("COSMOLKIT_PARITY_PROFILE").as_deref() {
        Ok("small" | "smiles_small") | Err(_) => "smiles_small".to_string(),
        Ok("strict" | "5000" | "smiles_5000") => "smiles_5000".to_string(),
        Ok(other) => panic!(
            "unknown COSMOLKIT_PARITY_PROFILE '{other}'; expected smiles_small or smiles_5000"
        ),
    }
}

pub fn smiles_path() -> PathBuf {
    if let Ok(path) = std::env::var("COSMOLKIT_PARITY_SMILES") {
        return PathBuf::from(path);
    }
    match profile_name().as_str() {
        "small" | "smiles_small" => repo_root().join("testdata/smiles/corpus/smiles_small.smi"),
        "strict" | "5000" | "smiles_5000" => {
            repo_root().join("testdata/smiles/corpus/smiles_5000.smi")
        }
        other => panic!(
            "unknown COSMOLKIT_PARITY_PROFILE '{other}'; expected smiles_small or smiles_5000"
        ),
    }
}

pub fn expected_family_dir(domain: &str, reference: &str) -> PathBuf {
    repo_root()
        .join("testdata")
        .join(domain)
        .join("expected")
        .join(reference)
        .join(profile_name())
}

pub fn rdkit_expected_domain(file_name: &str) -> &'static str {
    match file_name {
        "inchi.jsonl" => "inchi",
        "molblock_v2000_kekulized.jsonl"
        | "molblock_v2000_minimal.jsonl"
        | "molfile_read.jsonl" => "molblock",
        "sdf_read.jsonl" | "sdf_write.jsonl" => "sdf",
        "morgan_fingerprint.jsonl"
        | "atom_pair_fingerprint.jsonl"
        | "maccs_fingerprint.jsonl"
        | "rdkit_topological_fingerprint.jsonl"
        | "topological_torsion_fingerprint.jsonl"
        | "avalon_fingerprint.jsonl" => "fingerprint",
        "conformer_generation.jsonl"
        | "conformer_generation_library.jsonl"
        | "confseq_embed_template.jsonl" => "conformer",
        "forcefield_params.jsonl" | "mmff_builtin.jsonl" => "forcefield",
        "forcefield_coverage.jsonl" => "forcefield_coverage",
        "smiles_writer.jsonl" | "isomeric_smiles.jsonl" => "smiles",
        "svg_drawer.jsonl" | "prepared_draw_molecule.jsonl" => "depiction",
        "graph_features.jsonl" => "graph",
        "molecular_descriptors.jsonl" => "descriptors",
        "tetrahedral_stereo_geometry.jsonl" | "assign_atom_chiral_tags_from_structure.jsonl" => {
            "stereo"
        }
        "dg_bounds_matrix.jsonl" => "distgeom",
        "mol2_read.jsonl" => "mol2",
        "xyz_read.jsonl" => "xyz",
        "kekulize_clear_flags_false.jsonl" => "kekulize",
        "delete_substructs.jsonl" | "delete_substructs_onlyfrags_chirality.jsonl" => "substructure",
        "rdkit_builtin_fixture_migration.jsonl" => "rdkit_builtin",
        other => panic!(
            "unknown RDKit expected output '{other}'; add an explicit domain mapping before use"
        ),
    }
}

pub fn golden_path(file_name: &str) -> PathBuf {
    expected_path(rdkit_expected_domain(file_name), "rdkit", file_name)
}

pub fn expected_path(domain: &str, reference: &str, file_name: &str) -> PathBuf {
    expected_path_for_profile(domain, reference, &profile_name(), file_name)
}

pub fn expected_path_for_profile(
    domain: &str,
    reference: &str,
    profile: &str,
    file_name: &str,
) -> PathBuf {
    let family_dir = repo_root()
        .join("testdata")
        .join(domain)
        .join("expected")
        .join(reference)
        .join(profile);
    validate_expected_output_cached(&family_dir, domain, reference, file_name).unwrap_or_else(
        |error| {
            panic!(
                "invalid {reference} {domain} expected data: {error}; prepare it with `{}`",
                prepare_command(reference, domain, profile)
            )
        },
    );
    family_dir.join(file_name)
}

pub fn count_smiles_rows() -> usize {
    count_data_rows(&smiles_path()).unwrap_or_else(|error| panic!("{error}"))
}

pub fn rdkit_prepare_command(suite: &str) -> String {
    format!(
        ".venv/bin/python tools/testdata/rdkit/generate_all.py --python .venv/bin/python --profile {} --suite {suite} --jobs 4",
        profile_name()
    )
}

fn prepare_command(reference: &str, domain: &str, profile: &str) -> String {
    match reference {
        "rdkit" => rdkit_prepare_command(domain),
        "gemmi" => format!(
            ".venv/bin/python tools/testdata/gemmi/generate_all.py --profile {profile} --suite mmcif_writer --jobs 4"
        ),
        other => format!(
            ".venv/bin/python tools/testdata/{other}/generate_all.py --profile {profile} --suite {domain}"
        ),
    }
}

pub fn regenerate_command() -> String {
    rdkit_prepare_command("all")
}

/// Performs the preparation-time audit of every input and output in a family.
/// Ordinary tests should call [`expected_path`], which validates only the
/// output consumed by that test.
pub fn validate_expected_family(
    family_dir: &Path,
    expected_domain: &str,
    expected_reference: &str,
) -> Result<(), String> {
    let expected_profile = family_dir
        .file_name()
        .and_then(|name| name.to_str())
        .ok_or_else(|| {
            format!(
                "expected-data directory has no profile name: {}",
                family_dir.display()
            )
        })?;
    let manifest_path = family_dir.join("manifest.json");
    let manifest_bytes = std::fs::read(&manifest_path)
        .map_err(|error| format!("failed to read {}: {error}", manifest_path.display()))?;
    let manifest: ExpectedManifest = serde_json::from_slice(&manifest_bytes)
        .map_err(|error| format!("failed to parse {}: {error}", manifest_path.display()))?;

    validate_manifest_identity(
        &manifest,
        &manifest_path,
        expected_domain,
        expected_reference,
        expected_profile,
    )?;

    validate_current_input(&manifest.input, expected_profile)?;
    for output in &manifest.outputs {
        validate_output_identity(output, &manifest.input, expected_profile)?;
        validate_output_file(family_dir, output)?;
    }
    Ok(())
}

fn validate_expected_output_cached(
    family_dir: &Path,
    expected_domain: &str,
    expected_reference: &str,
    file_name: &str,
) -> ValidationResult {
    let key = family_dir.join(file_name);
    let cell = {
        let mut validations = OUTPUT_VALIDATIONS
            .get_or_init(|| Mutex::new(HashMap::new()))
            .lock()
            .map_err(|_| "expected-output validation cache lock is poisoned".to_string())?;
        validations
            .entry(key)
            .or_insert_with(|| Arc::new(OnceLock::new()))
            .clone()
    };
    cell.get_or_init(|| {
        validate_expected_output(family_dir, expected_domain, expected_reference, file_name)
    })
    .clone()
}

fn validate_expected_output(
    family_dir: &Path,
    expected_domain: &str,
    expected_reference: &str,
    file_name: &str,
) -> Result<(), String> {
    let expected_profile = family_dir
        .file_name()
        .and_then(|name| name.to_str())
        .ok_or_else(|| {
            format!(
                "expected-data directory has no profile name: {}",
                family_dir.display()
            )
        })?;
    let manifest_path = family_dir.join("manifest.json");
    let manifest_bytes = std::fs::read(&manifest_path)
        .map_err(|error| format!("failed to read {}: {error}", manifest_path.display()))?;
    let manifest: ExpectedManifest = serde_json::from_slice(&manifest_bytes)
        .map_err(|error| format!("failed to parse {}: {error}", manifest_path.display()))?;

    validate_manifest_identity(
        &manifest,
        &manifest_path,
        expected_domain,
        expected_reference,
        expected_profile,
    )?;
    validate_current_input(&manifest.input, expected_profile)?;

    let requested = Path::new(file_name);
    if requested.is_absolute()
        || requested
            .components()
            .any(|component| !matches!(component, std::path::Component::Normal(_)))
    {
        return Err(format!(
            "expected output path must be relative and normalized: {file_name}"
        ));
    }

    let output = manifest
        .outputs
        .iter()
        .find(|output| Path::new(&output.path) == requested)
        .ok_or_else(|| {
            format!(
                "{} does not declare requested output {file_name}",
                manifest_path.display()
            )
        })?;
    validate_output_identity(output, &manifest.input, expected_profile)?;
    validate_output_file(family_dir, output)
}

fn validate_manifest_identity(
    manifest: &ExpectedManifest,
    manifest_path: &Path,
    expected_domain: &str,
    expected_reference: &str,
    expected_profile: &str,
) -> Result<(), String> {
    if manifest.schema_version != 1 {
        return Err(format!(
            "{} has unsupported schema version {}",
            manifest_path.display(),
            manifest.schema_version
        ));
    }
    let reference_identity = reference_identity(expected_reference)?;
    if manifest.family != expected_reference
        || manifest.domain != expected_domain
        || manifest.profile != expected_profile
        || manifest.reference_implementation.name != reference_identity.name
        || manifest.reference_implementation.version != reference_identity.version
        || manifest.reference_runtime.implementation != reference_identity.runtime.implementation
        || manifest.reference_runtime.version != reference_identity.runtime.version
    {
        return Err(format!(
            "{} identity does not match family={expected_reference} domain={expected_domain} profile={}",
            manifest_path.display(),
            expected_profile
        ));
    }
    Ok(())
}

fn reference_identity(reference: &str) -> Result<ReferenceIdentity, String> {
    let path = repo_root()
        .join("testdata/reference")
        .join(format!("{reference}.json"));
    let bytes = std::fs::read(&path)
        .map_err(|error| format!("failed to read {}: {error}", path.display()))?;
    serde_json::from_slice(&bytes)
        .map_err(|error| format!("failed to parse {}: {error}", path.display()))
}

fn validate_current_input(input: &ManifestFile, expected_profile: &str) -> Result<(), String> {
    let current_path = match expected_profile {
        "bio_mmcif_writer" => repo_root().join("testdata/bio/gemmi_mmcif_writer_profile.json"),
        profile if profile == profile_name() => smiles_path(),
        "atom_pair_focused" => repo_root()
            .join("testdata/fingerprint/fixtures/rdkit/atom_pair_fingerprint_focused.smi"),
        "ciplabeler_focused" => {
            repo_root().join("testdata/stereo/fixtures/ciplabeler_focused.json")
        }
        "smiles_small" => repo_root().join("testdata/smiles/corpus/smiles_small.smi"),
        "smiles_5000" => repo_root().join("testdata/smiles/corpus/smiles_5000.smi"),
        other => return Err(format!("unknown expected-data profile '{other}'")),
    };
    let current = normalize_existing_path(&current_path)?;
    let declared = normalize_existing_path(&identity_path(input)?)?;
    if current != declared {
        return Err(format!(
            "manifest input {} does not match active corpus {}",
            declared.display(),
            current.display()
        ));
    }
    verify_checksum_cached(&declared, &input.sha256)
}

fn validate_output_identity(
    output: &ManifestOutput,
    primary_input: &ManifestFile,
    expected_profile: &str,
) -> Result<(), String> {
    if output.output_schema_version != 1 {
        return Err(format!(
            "{} has unsupported output schema version {}",
            output.path, output.output_schema_version
        ));
    }
    if output.generator.is_empty() {
        return Err(format!("{} declares no generator source", output.path));
    }
    for identity in &output.generator {
        verify_identity_file(identity)?;
    }
    if output.inputs.is_empty() {
        return Err(format!("{} declares no generation input", output.path));
    }
    for identity in &output.inputs {
        verify_identity_file(identity)?;
    }
    if !output
        .inputs
        .iter()
        .any(|input| input.path == primary_input.path && input.sha256 == primary_input.sha256)
    {
        return Err(format!(
            "{} does not include the active corpus identity",
            output.path
        ));
    }
    let expected_argument_prefix = [
        "--input".to_string(),
        primary_input.path.clone(),
        "--output".to_string(),
        output.path.clone(),
    ];
    let deterministic_shard_arguments_are_valid = match output.options.arguments.as_slice() {
        arguments if arguments == expected_argument_prefix => true,
        [prefix @ .., flag, count]
            if prefix == expected_argument_prefix
                && flag == "--shards"
                && count.parse::<usize>().is_ok_and(|count| count > 0) =>
        {
            true
        }
        _ => false,
    };
    if output.options.profile != expected_profile || !deterministic_shard_arguments_are_valid {
        return Err(format!(
            "{} generation options do not match the active profile",
            output.path
        ));
    }
    let (system, machine) = current_platform_identity();
    if output.platform.system != system || output.platform.machine != machine {
        return Err(format!(
            "{} platform identity is {}/{}, expected {system}/{machine}",
            output.path, output.platform.system, output.platform.machine
        ));
    }
    Ok(())
}

fn verify_identity_file(identity: &ManifestFile) -> Result<(), String> {
    let path = normalize_existing_path(&identity_path(identity)?)?;
    verify_checksum_cached(&path, &identity.sha256)
}

fn identity_path(identity: &ManifestFile) -> Result<PathBuf, String> {
    let path = Path::new(&identity.path);
    if path.is_absolute() {
        return Ok(path.to_path_buf());
    }
    if path
        .components()
        .any(|component| !matches!(component, std::path::Component::Normal(_)))
    {
        return Err(format!(
            "manifest identity path must be relative and normalized: {}",
            identity.path
        ));
    }
    Ok(repo_root().join(path))
}

fn normalize_existing_path(path: &Path) -> Result<PathBuf, String> {
    path.canonicalize()
        .map_err(|error| format!("failed to resolve {}: {error}", path.display()))
}

fn current_platform_identity() -> (&'static str, &'static str) {
    let system = match std::env::consts::OS {
        "linux" => "Linux",
        "macos" => "Darwin",
        "windows" => "Windows",
        other => other,
    };
    let machine = match (std::env::consts::OS, std::env::consts::ARCH) {
        ("macos", "aarch64") => "arm64",
        ("windows", "x86_64") => "AMD64",
        (_, other) => other,
    };
    (system, machine)
}

fn validate_output_file(family_dir: &Path, output: &ManifestOutput) -> Result<(), String> {
    let output_path = family_dir.join(&output.path);
    #[cfg(test)]
    {
        let mut counts = OUTPUT_SCAN_COUNTS
            .get_or_init(|| Mutex::new(HashMap::new()))
            .lock()
            .map_err(|_| "output scan-count lock is poisoned".to_string())?;
        *counts.entry(output_path.clone()).or_default() += 1;
    }
    let (checksum, records) = checksum_and_data_rows(&output_path)?;
    if checksum != output.sha256 {
        return Err(format!(
            "{} checksum is {checksum}, expected {}",
            output_path.display(),
            output.sha256
        ));
    }
    if records != output.records {
        return Err(format!(
            "{} record count is {records}, expected {}",
            output_path.display(),
            output.records
        ));
    }
    Ok(())
}

fn verify_checksum_cached(path: &Path, expected: &str) -> Result<(), String> {
    let cell = {
        let mut checksums = FILE_CHECKSUMS
            .get_or_init(|| Mutex::new(HashMap::new()))
            .lock()
            .map_err(|_| "file-checksum cache lock is poisoned".to_string())?;
        checksums
            .entry(path.to_path_buf())
            .or_insert_with(|| Arc::new(OnceLock::new()))
            .clone()
    };
    let actual = cell.get_or_init(|| sha256_file(path)).clone()?;
    if actual != expected {
        return Err(format!(
            "{} checksum is {actual}, expected {expected}",
            path.display()
        ));
    }
    Ok(())
}

fn lowercase_hex(bytes: &[u8]) -> String {
    const DIGITS: &[u8; 16] = b"0123456789abcdef";
    let mut encoded = String::with_capacity(bytes.len() * 2);
    for &byte in bytes {
        encoded.push(char::from(DIGITS[usize::from(byte >> 4)]));
        encoded.push(char::from(DIGITS[usize::from(byte & 0x0f)]));
    }
    encoded
}

fn checksum_and_data_rows(path: &Path) -> Result<(String, usize), String> {
    let file =
        File::open(path).map_err(|error| format!("failed to open {}: {error}", path.display()))?;
    let mut reader = BufReader::new(file);
    let mut digest = Sha256::new();
    let mut records = 0;
    let mut line = Vec::new();
    loop {
        line.clear();
        let count = reader
            .read_until(b'\n', &mut line)
            .map_err(|error| format!("failed to read {}: {error}", path.display()))?;
        if count == 0 {
            break;
        }
        digest.update(&line);
        let first_content = line
            .iter()
            .copied()
            .find(|byte| !byte.is_ascii_whitespace());
        if first_content.is_some_and(|byte| byte != b'#') {
            records += 1;
        }
    }
    Ok((lowercase_hex(digest.finalize().as_ref()), records))
}

fn sha256_file(path: &Path) -> Result<String, String> {
    let mut file =
        File::open(path).map_err(|error| format!("failed to open {}: {error}", path.display()))?;
    let mut digest = Sha256::new();
    let mut buffer = [0_u8; 64 * 1024];
    loop {
        let count = file
            .read(&mut buffer)
            .map_err(|error| format!("failed to read {}: {error}", path.display()))?;
        if count == 0 {
            break;
        }
        digest.update(&buffer[..count]);
    }
    Ok(lowercase_hex(digest.finalize().as_ref()))
}

fn count_data_rows(path: &Path) -> Result<usize, String> {
    let file =
        File::open(path).map_err(|error| format!("failed to open {}: {error}", path.display()))?;
    let mut count = 0;
    for line in BufReader::new(file).lines() {
        let line = line.map_err(|error| format!("failed to read {}: {error}", path.display()))?;
        let line = line.trim();
        if !line.is_empty() && !line.starts_with('#') {
            count += 1;
        }
    }
    Ok(count)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sha256_file_uses_the_standard_lowercase_hex_representation() {
        let path = std::env::temp_dir().join(format!(
            "cosmolkit-test-support-sha256-{}",
            std::process::id()
        ));
        std::fs::write(&path, b"abc").expect("checksum fixture should be written");
        assert_eq!(
            sha256_file(&path).expect("fixture checksum should be available"),
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );
        std::fs::remove_file(path).expect("checksum fixture should be removed");
    }

    #[test]
    fn expected_data_json_floats_round_trip_to_exact_binary64_values() {
        for (literal, expected_bits) in [
            ("-1.3006687716310001", 0xbff4_cf8a_0ed1_5695),
            ("-0.39018908733716867", 0xbfd8_f8db_a657_a16b),
            ("1.3682138422726557", 0x3ff5_e434_32a7_edcf),
            ("1.8158026507898606", 0x3ffd_0d87_1492_1ef7),
        ] {
            let standard = literal
                .parse::<f64>()
                .expect("fixture float should parse with the standard library");
            let json =
                serde_json::from_str::<f64>(literal).expect("fixture float should parse as JSON");

            assert_eq!(standard.to_bits(), expected_bits, "literal {literal}");
            assert_eq!(json.to_bits(), expected_bits, "literal {literal}");
        }
    }

    fn test_manifest(
        domain: &str,
        output_name: &str,
        output_path: &Path,
        generator_sha256: String,
    ) -> serde_json::Value {
        let corpus = repo_root().join("testdata/smiles/corpus/smiles_small.smi");
        let corpus_sha256 = sha256_file(&corpus).expect("corpus checksum should be available");
        let (system, machine) = current_platform_identity();
        let reference = reference_identity("rdkit").expect("RDKit identity should be available");
        serde_json::json!({
            "schema_version": 1,
            "family": "rdkit",
            "domain": domain,
            "profile": profile_name(),
            "reference_implementation": {
                "name": reference.name,
                "version": reference.version
            },
            "reference_runtime": {
                "implementation": reference.runtime.implementation,
                "version": reference.runtime.version
            },
            "input": {
                "path": "testdata/smiles/corpus/smiles_small.smi",
                "sha256": corpus_sha256
            },
            "outputs": [{
                "path": output_name,
                "output_schema_version": 1,
                "generator": [{
                    "path": "testdata/smiles/corpus/smiles_small.smi",
                    "sha256": generator_sha256
                }],
                "inputs": [{
                    "path": "testdata/smiles/corpus/smiles_small.smi",
                    "sha256": sha256_file(&corpus).expect("corpus checksum should be available")
                }],
                "options": {
                    "profile": profile_name(),
                    "arguments": [
                        "--input",
                        "testdata/smiles/corpus/smiles_small.smi",
                        "--output",
                        output_name
                    ]
                },
                "platform": {"system": system, "machine": machine},
                "sha256": sha256_file(output_path)
                    .expect("output checksum should be available"),
                "records": count_data_rows(output_path)
                    .expect("output rows should be countable")
            }]
        })
    }

    fn temporary_family_dir(label: &str) -> (PathBuf, PathBuf) {
        let unique = format!(
            "cosmolkit-test-support-{label}-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .expect("system clock should follow UNIX epoch")
                .as_nanos()
        );
        let temporary_root = std::env::temp_dir().join(unique);
        let family_dir = temporary_root.join(profile_name());
        std::fs::create_dir_all(&family_dir).expect("temporary family should be created");
        (temporary_root, family_dir)
    }

    #[test]
    fn profile_paths_are_repository_owned() {
        assert!(smiles_path().starts_with(repo_root().join("testdata")));
        assert!(!smiles_path().starts_with(repo_root().join("third_party")));
    }

    #[test]
    fn output_lookup_does_not_validate_unrelated_outputs() {
        let (temporary_root, family_dir) = temporary_family_dir("unrelated-output");

        let requested_path = family_dir.join("requested.jsonl");
        std::fs::write(&requested_path, b"{\"case\":1}\n")
            .expect("requested output should be written");
        let unrelated_path = family_dir.join("unrelated.jsonl");
        std::fs::write(&unrelated_path, b"{\"case\":2}\n")
            .expect("unrelated output should be written");
        let corpus = repo_root().join("testdata/smiles/corpus/smiles_small.smi");
        let mut manifest = test_manifest(
            "test_domain",
            "requested.jsonl",
            &requested_path,
            sha256_file(&corpus).expect("corpus checksum should be available"),
        );
        manifest["outputs"]
            .as_array_mut()
            .expect("outputs should be an array")
            .push(serde_json::json!({
                "path": "unrelated.jsonl",
                "output_schema_version": 1,
                "generator": [],
                "inputs": [],
                "options": {"profile": profile_name(), "arguments": []},
                "platform": {"system": "invalid", "machine": "invalid"},
                "sha256": "invalid-on-purpose",
                "records": 1
            }));
        std::fs::write(
            family_dir.join("manifest.json"),
            serde_json::to_vec(&manifest).expect("manifest should serialize"),
        )
        .expect("manifest should be written");

        validate_expected_output(&family_dir, "test_domain", "rdkit", "requested.jsonl")
            .expect("requested output should validate independently");
        assert!(validate_expected_family(&family_dir, "test_domain", "rdkit").is_err());

        std::fs::remove_dir_all(temporary_root).expect("temporary family should be removed");
    }

    #[test]
    fn output_lookup_rejects_stale_generator_identity() {
        let (temporary_root, family_dir) = temporary_family_dir("stale-generator");
        let output_path = family_dir.join("requested.jsonl");
        std::fs::write(&output_path, b"{\"case\":1}\n").expect("output should be written");
        let manifest = test_manifest(
            "test_domain",
            "requested.jsonl",
            &output_path,
            "stale-generator-checksum".to_string(),
        );
        std::fs::write(
            family_dir.join("manifest.json"),
            serde_json::to_vec(&manifest).expect("manifest should serialize"),
        )
        .expect("manifest should be written");

        let error =
            validate_expected_output(&family_dir, "test_domain", "rdkit", "requested.jsonl")
                .expect_err("stale generator identity must be rejected");
        assert!(error.contains("checksum"), "unexpected error: {error}");
        std::fs::remove_dir_all(temporary_root).expect("temporary family should be removed");
    }

    #[test]
    fn cached_output_lookup_scans_requested_output_once() {
        let (temporary_root, family_dir) = temporary_family_dir("cached-output");
        let output_path = family_dir.join("requested.jsonl");
        std::fs::write(&output_path, b"{\"case\":1}\n\n# comment\n")
            .expect("output should be written");
        let corpus = repo_root().join("testdata/smiles/corpus/smiles_small.smi");
        let manifest = test_manifest(
            "test_domain",
            "requested.jsonl",
            &output_path,
            sha256_file(&corpus).expect("corpus checksum should be available"),
        );
        std::fs::write(
            family_dir.join("manifest.json"),
            serde_json::to_vec(&manifest).expect("manifest should serialize"),
        )
        .expect("manifest should be written");

        validate_expected_output_cached(&family_dir, "test_domain", "rdkit", "requested.jsonl")
            .expect("first validation should pass");
        validate_expected_output_cached(&family_dir, "test_domain", "rdkit", "requested.jsonl")
            .expect("cached validation should pass");
        let scans = OUTPUT_SCAN_COUNTS
            .get()
            .expect("scan counts should exist")
            .lock()
            .expect("scan-count lock should not be poisoned");
        assert_eq!(scans.get(&output_path), Some(&1));
        drop(scans);
        std::fs::remove_dir_all(temporary_root).expect("temporary family should be removed");
    }
}
