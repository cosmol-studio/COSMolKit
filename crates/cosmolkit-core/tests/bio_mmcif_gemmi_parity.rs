use std::path::{Path, PathBuf};

use cosmolkit_core::{BioCoorFormat, BioStructure};
use cosmolkit_test_support::{expected_path_for_profile, repo_root};
use serde::Deserialize;

const PROFILE: &str = "bio_mmcif_writer";

#[derive(Debug, Deserialize)]
struct WriterProfile {
    schema_version: u32,
    profile: String,
    gemmi_commit: String,
    cases: Vec<WriterCase>,
}

#[derive(Debug, Deserialize)]
struct WriterCase {
    case_id: String,
    input: String,
    input_format: String,
    arguments: Vec<String>,
    output: String,
}

fn load_profile() -> WriterProfile {
    let path = repo_root().join("testdata/bio/gemmi_mmcif_writer_profile.json");
    let bytes = std::fs::read(&path).unwrap_or_else(|error| panic!("failed to read {}: {error}", path.display()));
    let profile: WriterProfile =
        serde_json::from_slice(&bytes).unwrap_or_else(|error| panic!("failed to parse {}: {error}", path.display()));
    assert_eq!(profile.schema_version, 1);
    assert_eq!(profile.profile, PROFILE);
    assert_eq!(profile.gemmi_commit, "5cc1c23c6007e0e6cbd69289c6f7c0bff50e943e");
    profile
}

fn read_case(case: &WriterCase) -> BioStructure {
    let input = repo_root().join(&case.input);
    match case.input_format.as_str() {
        "pdb" => BioStructure::from_pdb(&input),
        "mmcif" => BioStructure::from_mmcif(&input),
        other => panic!("{} declares unsupported input format {other}", case.case_id),
    }
    .unwrap_or_else(|error| panic!("failed to read {}: {error}", input.display()))
}

fn assert_semantic_state_eq(actual: &BioStructure, expected: &BioStructure, case_id: &str) {
    assert_eq!(actual.name(), expected.name(), "{case_id}: structure name");
    assert_eq!(actual.atoms(), expected.atoms(), "{case_id}: atoms");
    assert_eq!(actual.residues(), expected.residues(), "{case_id}: residues");
    assert_eq!(actual.chains(), expected.chains(), "{case_id}: chains");
    assert_eq!(actual.entities(), expected.entities(), "{case_id}: entities");
    assert_eq!(actual.models(), expected.models(), "{case_id}: models");
    assert_eq!(actual.coordinates(), expected.coordinates(), "{case_id}: coordinates");
    assert_eq!(actual.metadata(), expected.metadata(), "{case_id}: metadata");
    assert_eq!(actual.crystal(), expected.crystal(), "{case_id}: crystal");
    assert_eq!(actual.resolution(), expected.resolution(), "{case_id}: resolution");
    assert_eq!(actual.ter_status(), expected.ter_status(), "{case_id}: TER status");
    assert_eq!(
        actual.mod_residues(),
        expected.mod_residues(),
        "{case_id}: modified residues"
    );
    assert_eq!(actual.connections(), expected.connections(), "{case_id}: connections");
    assert_eq!(actual.cispeps(), expected.cispeps(), "{case_id}: cis peptides");
    assert_eq!(actual.assemblies(), expected.assemblies(), "{case_id}: assemblies");
    assert_eq!(actual.has_origx(), expected.has_origx(), "{case_id}: ORIGX presence");
    assert_eq!(actual.origx(), expected.origx(), "{case_id}: ORIGX transform");
    assert_eq!(
        actual.ncs_operators(),
        expected.ncs_operators(),
        "{case_id}: NCS operators"
    );
    assert_eq!(
        actual.ncs_oper_identity_id(),
        expected.ncs_oper_identity_id(),
        "{case_id}: NCS identity"
    );
}

fn read_generated_mmcif(text: &str, path: &Path) -> BioStructure {
    BioStructure::from_str_with_format(text, &path.to_string_lossy(), BioCoorFormat::Mmcif)
        .unwrap_or_else(|error| panic!("failed to reparse {}: {error}", path.display()))
}

fn assert_exact_output(actual: &str, expected: &str, case_id: &str) {
    if actual == expected {
        return;
    }
    let offset = actual
        .as_bytes()
        .iter()
        .zip(expected.as_bytes())
        .position(|(actual, expected)| actual != expected)
        .unwrap_or_else(|| actual.len().min(expected.len()));
    let line = actual.as_bytes()[..offset.min(actual.len())]
        .iter()
        .filter(|&&byte| byte == b'\n')
        .count()
        + 1;
    let line_start = actual.as_bytes()[..offset.min(actual.len())]
        .iter()
        .rposition(|&byte| byte == b'\n')
        .map_or(0, |index| index + 1);
    let column = offset.saturating_sub(line_start) + 1;
    let context_start = offset.saturating_sub(80);
    let actual_end = (offset + 160).min(actual.len());
    let expected_end = (offset + 160).min(expected.len());
    panic!(
        "{case_id}: exact Gemmi output differs at byte {offset}, line {line}, column {column}\nactual context: {:?}\nGemmi context: {:?}\nactual bytes: {}; Gemmi bytes: {}",
        &actual[context_start.min(actual.len())..actual_end],
        &expected[context_start.min(expected.len())..expected_end],
        actual.len(),
        expected.len()
    );
}

#[test]
fn mmcif_writer_matches_every_prepared_pinned_gemmi_fixture_exactly_and_semantically() {
    let profile = load_profile();
    assert!(!profile.cases.is_empty(), "writer profile must declare cases");

    for case in &profile.cases {
        assert!(
            case.arguments.is_empty(),
            "{}: default writer case expected",
            case.case_id
        );
        let source = read_case(case);
        let actual = source
            .to_mmcif()
            .unwrap_or_else(|error| panic!("{} failed to write: {error}", case.case_id));
        let expected_path = expected_path_for_profile("bio", "gemmi", PROFILE, &case.output);
        let expected = std::fs::read_to_string(&expected_path)
            .unwrap_or_else(|error| panic!("failed to read {}: {error}", expected_path.display()));

        assert_exact_output(&actual, &expected, &case.case_id);

        let actual_path = PathBuf::from(format!("{}.cosmolkit.cif", case.case_id));
        let actual_reparsed = read_generated_mmcif(&actual, &actual_path);
        let expected_reparsed = read_generated_mmcif(&expected, &expected_path);
        assert_semantic_state_eq(&actual_reparsed, &expected_reparsed, &case.case_id);
    }
}
