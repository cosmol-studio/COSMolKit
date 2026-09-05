use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

use cosmolkit_core::fingerprint::{
    AdditionalOutput, FingerprintFuncArguments, TopologicalTorsionArguments,
    TypedFingerprintGenerator, generatorFromJSON, generatorToJSON, getTopologicalTorsionGenerator,
};
use cosmolkit_core::{
    AtomId, Molecule, TopologicalTorsionLegacyKind, TopologicalTorsionLegacyParams,
    TopologicalTorsionLegacyResult, get_atom_code, get_topological_torsion_code,
    topological_torsion_legacy_fingerprint,
};
use rayon::prelude::*;
use serde::Deserialize;
use serde_json::Value;

mod common;
use common::parity_data;

const PROFILE_NAMES: [&str; 12] = [
    "all_provenance",
    "chirality",
    "count_simulation_provenance",
    "custom_atom_invariants",
    "default",
    "extra_bits",
    "folding_collisions",
    "from_first_atom",
    "ignore_first_atom",
    "shortest_paths",
    "three_atom_custom_count_bounds",
    "unfolded_no_count_simulation",
];

const CORPUS_PROFILE_NAMES: [&str; 9] = [
    "all_provenance",
    "chirality",
    "count_simulation_provenance",
    "custom_atom_invariants",
    "default",
    "extra_bits_folding",
    "five_atom_shortest_paths",
    "three_atom_custom_count_bounds",
    "unfolded_no_count_simulation",
];

#[derive(Debug, Deserialize)]
struct GoldenRecord {
    row: usize,
    smiles: String,
    rdkit_ok: bool,
    profiles: BTreeMap<String, GoldenProfile>,
    legacy: Option<GoldenLegacy>,
    helpers: Option<GoldenHelpers>,
    errors: Option<BTreeMap<String, GoldenError>>,
    error: Option<String>,
}

#[derive(Debug, Deserialize)]
struct GoldenProfile {
    parameters: Value,
    info_string: String,
    json: Value,
    json_restored_count: GoldenCount,
    sparse_count: GoldenCount,
    sparse_bit: GoldenBits,
    count: GoldenCount,
    bit: GoldenBits,
    bulk: GoldenBulk,
    additional_output: Option<GoldenOutputs>,
}

#[derive(Debug, Deserialize)]
struct GoldenBulk {
    sparse_count: GoldenCount,
    sparse_bit: GoldenBits,
    count: GoldenCount,
    bit: GoldenBits,
}

#[derive(Debug, Deserialize)]
struct GoldenCount {
    size: u64,
    nonzero_elements: BTreeMap<String, i32>,
}

#[derive(Debug, Deserialize)]
struct GoldenBits {
    size: u64,
    on_bits: Vec<u64>,
}

#[derive(Debug, Deserialize)]
struct GoldenOutputs {
    sparse_count: GoldenAdditionalOutput,
    sparse_bit: GoldenAdditionalOutput,
    count: GoldenAdditionalOutput,
    bit: GoldenAdditionalOutput,
}

#[derive(Debug, Deserialize)]
struct GoldenAdditionalOutput {
    atom_to_bits: Vec<Vec<u64>>,
    bit_info_map: BTreeMap<String, Vec<[u32; 2]>>,
    bit_paths: BTreeMap<String, Vec<Vec<usize>>>,
    atom_counts: Vec<u32>,
    atoms_per_bit: BTreeMap<String, Vec<Vec<usize>>>,
}

#[derive(Debug, Deserialize)]
struct GoldenLegacy {
    unfolded: GoldenCount,
    unfolded_chiral: GoldenCount,
    unfolded_custom: GoldenCount,
    hashed: GoldenCount,
    hashed_rooted: GoldenCount,
    hashed_ignored: GoldenCount,
    bit_vectors: BTreeMap<String, GoldenBits>,
}

#[derive(Debug, Deserialize)]
struct GoldenHelpers {
    atom_codes: Vec<GoldenAtomCode>,
    paths: Vec<GoldenPath>,
    ids: Vec<u64>,
}

#[derive(Debug, Deserialize)]
struct GoldenAtomCode {
    code: u32,
    chiral_code: u32,
}

#[derive(Debug, Deserialize)]
struct GoldenPath {
    path: Vec<usize>,
    score: u64,
}

#[derive(Debug, Deserialize)]
struct GoldenError {
    #[serde(rename = "type")]
    error_type: Option<String>,
    message: Option<String>,
}

fn golden_reader() -> BufReader<File> {
    let path = parity_data::golden_path("topological_torsion_fingerprint.jsonl");
    assert!(
        path.exists(),
        "missing RDKit Topological Torsion golden: {}. Prepare it with \
         .venv/bin/python tools/testdata/rdkit/generate_all.py --python \
         .venv/bin/python --profile {} --suite fingerprint --jobs 4",
        path.display(),
        parity_data::profile_name(),
    );
    BufReader::new(
        File::open(&path)
            .unwrap_or_else(|error| panic!("failed to open {}: {error}", path.display())),
    )
}

fn corpus_smiles() -> Vec<String> {
    let path = parity_data::smiles_path();
    std::fs::read_to_string(&path)
        .unwrap_or_else(|error| panic!("failed to read {}: {error}", path.display()))
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(str::to_owned)
        .collect()
}

fn value_bool(parameters: &Value, name: &str, default: bool) -> bool {
    parameters
        .get(name)
        .and_then(Value::as_bool)
        .unwrap_or(default)
}

fn value_u32(parameters: &Value, name: &str, default: u32) -> u32 {
    parameters
        .get(name)
        .and_then(Value::as_u64)
        .map(|value| u32::try_from(value).expect("profile integer fits u32"))
        .unwrap_or(default)
}

fn profile_generator(
    profile: &GoldenProfile,
) -> cosmolkit_core::TopologicalTorsionFingerprintGenerator {
    let parameters = &profile.parameters;
    let count_bounds = parameters
        .get("countBounds")
        .and_then(Value::as_array)
        .map(|values| {
            values
                .iter()
                .map(|value| value.as_u64().expect("count bound") as u32)
                .collect()
        })
        .unwrap_or_else(|| vec![1, 2, 4, 8]);
    let mut arguments = TopologicalTorsionArguments::new(
        value_bool(parameters, "includeChirality", false),
        value_u32(parameters, "torsionAtomCount", 4),
        value_bool(parameters, "countSimulation", true),
        count_bounds,
        value_u32(parameters, "fpSize", 2048),
    )
    .expect("valid oracle profile arguments");
    arguments.df_only_shortest_paths = value_bool(parameters, "onlyShortestPaths", false);
    arguments.fingerprint_arguments.d_num_bits_per_feature =
        value_u32(parameters, "numBitsPerFeature", 1);
    getTopologicalTorsionGenerator(&arguments, None, true).expect("valid oracle generator")
}

fn profile_call_arguments(
    molecule: &Molecule,
    profile: &GoldenProfile,
    provenance: bool,
) -> FingerprintFuncArguments {
    let parameters = &profile.parameters;
    let mut arguments = FingerprintFuncArguments {
        from_atoms: (parameters.get("fromAtoms").and_then(Value::as_str) == Some("first")
            && molecule.num_atoms() > 0)
            .then(|| vec![0]),
        ignore_atoms: (parameters.get("ignoreAtoms").and_then(Value::as_str) == Some("first")
            && molecule.num_atoms() > 0)
            .then(|| vec![0]),
        custom_atom_invariants: (parameters
            .get("customAtomInvariants")
            .and_then(Value::as_str)
            == Some("index_plus_17"))
        .then(|| {
            (0..molecule.num_atoms())
                .map(|index| index as u32 + 17)
                .collect()
        }),
        ..Default::default()
    };
    if provenance {
        let mut output = AdditionalOutput::new();
        output.allocate_atom_to_bits();
        output.allocate_bit_info_map();
        output.allocate_bit_paths();
        output.allocate_atom_counts();
        output.allocate_atoms_per_bit();
        arguments.additional_output = Some(output);
    }
    arguments
}

fn expected_counts(golden: &GoldenCount) -> BTreeMap<u64, i32> {
    golden
        .nonzero_elements
        .iter()
        .map(|(bit, count)| (bit.parse().expect("u64 bit id"), *count))
        .collect()
}

fn expected_paths(golden: &BTreeMap<String, Vec<Vec<usize>>>) -> BTreeMap<u64, Vec<Vec<usize>>> {
    golden
        .iter()
        .map(|(bit, paths)| (bit.parse().expect("u64 bit id"), paths.clone()))
        .collect()
}

fn expected_bit_info(golden: &BTreeMap<String, Vec<[u32; 2]>>) -> BTreeMap<u64, Vec<(u32, u32)>> {
    golden
        .iter()
        .map(|(bit, entries)| {
            (
                bit.parse().expect("u64 bit id"),
                entries.iter().map(|entry| (entry[0], entry[1])).collect(),
            )
        })
        .collect()
}

fn assert_output(context: &str, actual: &AdditionalOutput, expected: &GoldenAdditionalOutput) {
    assert_eq!(
        actual.atom_to_bits.as_ref(),
        Some(&expected.atom_to_bits),
        "{context}"
    );
    assert_eq!(
        actual.atom_counts.as_ref(),
        Some(&expected.atom_counts),
        "{context}"
    );
    assert_eq!(
        actual.bit_info_map.as_ref(),
        Some(&expected_bit_info(&expected.bit_info_map)),
        "{context}"
    );
    assert_eq!(
        actual.bit_paths.as_ref(),
        Some(&expected_paths(&expected.bit_paths)),
        "{context}"
    );
    assert_eq!(
        actual.atoms_per_bit.as_ref(),
        Some(&expected_paths(&expected.atoms_per_bit)),
        "{context}"
    );
}

fn assert_modern_profile(row: usize, name: &str, molecule: &Molecule, expected: &GoldenProfile) {
    let context = format!("row {row} profile {name}");
    let generator = profile_generator(expected);
    let info_string = format!(
        "{} --- {} --- {} --- {} --- No bond invariants generator",
        generator
            .fingerprint_arguments
            .fingerprint_arguments
            .common_arguments_string(),
        generator.fingerprint_arguments.infoString(),
        generator.atom_environment_generator.infoString(),
        generator.atom_invariants_generator.infoString(),
    );
    assert_eq!(info_string, expected.info_string, "{context}");
    let actual_json: Value = serde_json::from_str(&generatorToJSON(
        &TypedFingerprintGenerator::TopologicalTorsion(generator.clone()),
    ))
    .expect("Rust generator JSON");
    assert_eq!(actual_json, expected.json, "{context}: JSON mismatch");

    let provenance = expected.additional_output.is_some();
    let mut sparse_count_args = profile_call_arguments(molecule, expected, provenance);
    let sparse_count = generator
        .getSparseCountFingerprint(molecule, &mut sparse_count_args)
        .unwrap_or_else(|error| panic!("{context}: sparse count failed: {error}"));
    assert_eq!(sparse_count.size(), expected.sparse_count.size, "{context}");
    assert_eq!(
        sparse_count.nonzero_elements(),
        &expected_counts(&expected.sparse_count),
        "{context}"
    );

    let mut sparse_bit_args = profile_call_arguments(molecule, expected, provenance);
    let sparse_bit = generator
        .getSparseFingerprint(molecule, &mut sparse_bit_args)
        .unwrap_or_else(|error| panic!("{context}: sparse bit failed: {error}"));
    assert_eq!(sparse_bit.size(), expected.sparse_bit.size, "{context}");
    assert_eq!(
        sparse_bit.on_bits().iter().copied().collect::<Vec<_>>(),
        expected.sparse_bit.on_bits,
        "{context}"
    );

    let mut count_args = profile_call_arguments(molecule, expected, provenance);
    let count = generator
        .getCountFingerprint(molecule, &mut count_args)
        .unwrap_or_else(|error| panic!("{context}: count failed: {error}"));
    assert_eq!(count.size(), expected.count.size, "{context}");
    assert_eq!(
        count.nonzero_elements(),
        &expected_counts(&expected.count),
        "{context}"
    );

    let mut bit_args = profile_call_arguments(molecule, expected, provenance);
    let bit = generator
        .getFingerprint(molecule, &mut bit_args)
        .unwrap_or_else(|error| panic!("{context}: bit failed: {error}"));
    assert_eq!(bit.n_bits() as u64, expected.bit.size, "{context}");
    assert_eq!(
        bit.on_bits()
            .into_iter()
            .map(|bit| bit as u64)
            .collect::<Vec<_>>(),
        expected.bit.on_bits,
        "{context}"
    );

    if let Some(outputs) = &expected.additional_output {
        assert_output(
            &format!("{context}: sparse_count provenance"),
            sparse_count_args
                .additional_output
                .as_ref()
                .expect("output"),
            &outputs.sparse_count,
        );
        assert_output(
            &format!("{context}: sparse_bit provenance"),
            sparse_bit_args.additional_output.as_ref().expect("output"),
            &outputs.sparse_bit,
        );
        assert_output(
            &format!("{context}: count provenance"),
            count_args.additional_output.as_ref().expect("output"),
            &outputs.count,
        );
        assert_output(
            &format!("{context}: bit provenance"),
            bit_args.additional_output.as_ref().expect("output"),
            &outputs.bit,
        );
    }

    let references = [Some(molecule)];
    let bulk_sparse_count = generator
        .getSparseCountFingerprints(&references, 2)
        .expect("bulk sparse count")
        .pop()
        .flatten()
        .expect("bulk item");
    assert_eq!(
        bulk_sparse_count.size(),
        expected.bulk.sparse_count.size,
        "{context}"
    );
    assert_eq!(
        bulk_sparse_count.nonzero_elements(),
        &expected_counts(&expected.bulk.sparse_count),
        "{context}"
    );
    let bulk_sparse_bit = generator
        .getSparseFingerprints(&references, 2)
        .expect("bulk sparse bit")
        .pop()
        .flatten()
        .expect("bulk item");
    assert_eq!(
        bulk_sparse_bit.size(),
        expected.bulk.sparse_bit.size,
        "{context}"
    );
    assert_eq!(
        bulk_sparse_bit
            .on_bits()
            .iter()
            .copied()
            .collect::<Vec<_>>(),
        expected.bulk.sparse_bit.on_bits,
        "{context}"
    );
    let bulk_count = generator
        .getCountFingerprints(&references, 2)
        .expect("bulk count")
        .pop()
        .flatten()
        .expect("bulk item");
    assert_eq!(bulk_count.size(), expected.bulk.count.size, "{context}");
    assert_eq!(
        bulk_count.nonzero_elements(),
        &expected_counts(&expected.bulk.count),
        "{context}"
    );
    let bulk_bit = generator
        .getFingerprints(&references, 2)
        .expect("bulk bit")
        .pop()
        .flatten()
        .expect("bulk item");
    assert_eq!(
        bulk_bit.n_bits() as u64,
        expected.bulk.bit.size,
        "{context}"
    );
    assert_eq!(
        bulk_bit
            .on_bits()
            .into_iter()
            .map(|bit| bit as u64)
            .collect::<Vec<_>>(),
        expected.bulk.bit.on_bits,
        "{context}"
    );

    let restored = generatorFromJSON(&expected.json.to_string()).expect("restore oracle JSON");
    let TypedFingerprintGenerator::TopologicalTorsion(restored) = restored else {
        panic!("{context}: restored wrong generator family")
    };
    let mut restored_args = profile_call_arguments(molecule, expected, false);
    let restored_count = restored
        .getCountFingerprint(molecule, &mut restored_args)
        .expect("restored count");
    assert_eq!(
        restored_count.nonzero_elements(),
        &expected_counts(&expected.json_restored_count),
        "{context}: restored JSON output"
    );
}

fn assert_legacy(row: usize, molecule: &Molecule, expected: &GoldenLegacy) {
    let context = format!("row {row} legacy");
    let cases = [
        (
            "unfolded",
            TopologicalTorsionLegacyParams::default(),
            &expected.unfolded,
        ),
        (
            "unfolded_chiral",
            TopologicalTorsionLegacyParams {
                include_chirality: true,
                ..Default::default()
            },
            &expected.unfolded_chiral,
        ),
        (
            "unfolded_custom",
            TopologicalTorsionLegacyParams {
                atom_invariants: Some(
                    (0..molecule.num_atoms())
                        .map(|index| index as u32 + 17)
                        .collect(),
                ),
                ..Default::default()
            },
            &expected.unfolded_custom,
        ),
        (
            "hashed",
            TopologicalTorsionLegacyParams {
                kind: TopologicalTorsionLegacyKind::HashedCount,
                n_bits: 1000,
                ..Default::default()
            },
            &expected.hashed,
        ),
        (
            "hashed_rooted",
            TopologicalTorsionLegacyParams {
                kind: TopologicalTorsionLegacyKind::HashedCount,
                n_bits: 1000,
                from_atoms: (molecule.num_atoms() > 0).then(|| vec![0]),
                ..Default::default()
            },
            &expected.hashed_rooted,
        ),
        (
            "hashed_ignored",
            TopologicalTorsionLegacyParams {
                kind: TopologicalTorsionLegacyKind::HashedCount,
                n_bits: 1000,
                ignore_atoms: (molecule.num_atoms() > 0).then(|| vec![0]),
                ..Default::default()
            },
            &expected.hashed_ignored,
        ),
    ];
    for (name, params, golden) in cases {
        let result = topological_torsion_legacy_fingerprint(molecule, &params)
            .unwrap_or_else(|error| panic!("{context} {name}: {error}"));
        let TopologicalTorsionLegacyResult::SparseCount(actual) = result else {
            panic!("{context} {name}: wrong vector type")
        };
        assert_eq!(actual.size(), golden.size, "{context} {name}");
        assert_eq!(
            actual.nonzero_elements(),
            &expected_counts(golden),
            "{context} {name}"
        );
    }
    for (bits_per_entry, golden) in &expected.bit_vectors {
        let result = topological_torsion_legacy_fingerprint(
            molecule,
            &TopologicalTorsionLegacyParams {
                kind: TopologicalTorsionLegacyKind::HashedBit,
                n_bits: 256,
                n_bits_per_entry: bits_per_entry.parse().expect("bits per entry"),
                ..Default::default()
            },
        )
        .unwrap_or_else(|error| panic!("{context} bit {bits_per_entry}: {error}"));
        let TopologicalTorsionLegacyResult::Bit(actual) = result else {
            panic!("{context} bit {bits_per_entry}: wrong vector type")
        };
        assert_eq!(actual.n_bits() as u64, golden.size, "{context}");
        assert_eq!(
            actual
                .on_bits()
                .into_iter()
                .map(|bit| bit as u64)
                .collect::<Vec<_>>(),
            golden.on_bits,
            "{context} bit {bits_per_entry}"
        );
    }
}

fn assert_helpers(row: usize, molecule: &Molecule, expected: &GoldenHelpers) {
    assert_eq!(expected.atom_codes.len(), molecule.num_atoms(), "row {row}");
    for (index, golden) in expected.atom_codes.iter().enumerate() {
        assert_eq!(
            get_atom_code(molecule, AtomId::new(index), 0, false).expect("atom code"),
            golden.code,
            "row {row} atom {index}"
        );
        assert_eq!(
            get_atom_code(molecule, AtomId::new(index), 0, true).expect("chiral atom code"),
            golden.chiral_code,
            "row {row} chiral atom {index}"
        );
    }
    for golden in &expected.paths {
        let codes = golden
            .path
            .iter()
            .enumerate()
            .map(|(position, index)| {
                let subtract = if position == 0 || position + 1 == golden.path.len() {
                    1
                } else {
                    2
                };
                get_atom_code(molecule, AtomId::new(*index), subtract, false).expect("path code")
            })
            .collect::<Vec<_>>();
        assert_eq!(
            get_topological_torsion_code(&codes, false).expect("torsion score"),
            golden.score,
            "row {row} path {:?}",
            golden.path
        );
    }
    let unfolded = topological_torsion_legacy_fingerprint(
        molecule,
        &TopologicalTorsionLegacyParams::default(),
    )
    .expect("unfolded ids");
    let TopologicalTorsionLegacyResult::SparseCount(unfolded) = unfolded else {
        panic!("row {row}: wrong ids vector")
    };
    let ids = unfolded
        .nonzero_elements()
        .iter()
        .flat_map(|(bit, count)| std::iter::repeat_n(*bit, *count as usize))
        .collect::<Vec<_>>();
    assert_eq!(ids, expected.ids, "row {row}: helper ids");
}

fn assert_error_profile(row: usize, molecule: &Molecule, errors: &BTreeMap<String, GoldenError>) {
    assert_eq!(
        errors.keys().map(String::as_str).collect::<Vec<_>>(),
        vec![
            "empty_count_bounds",
            "invalid_json",
            "oversized_torsion",
            "short_custom_invariants",
            "short_path",
        ],
        "row {row}: error profile set changed"
    );
    assert_eq!(
        errors["invalid_json"].error_type.as_deref(),
        Some("RuntimeError")
    );
    assert!(
        errors["invalid_json"]
            .message
            .as_deref()
            .is_some_and(|value| value.contains("expected 'null'"))
    );
    assert_eq!(
        errors["short_path"].error_type.as_deref(),
        Some("IndexError")
    );
    assert_eq!(
        errors["short_path"].message.as_deref(),
        Some("list index out of range")
    );
    assert!(errors["oversized_torsion"].error_type.is_none());
    assert!(errors["empty_count_bounds"].error_type.is_none());
    assert!(errors["short_custom_invariants"].error_type.is_none());

    assert!(
        generatorFromJSON("not json").is_err(),
        "row {row}: invalid JSON"
    );
    let oversized = TopologicalTorsionArguments::new(false, 8, true, vec![1, 2, 4, 8], 2048)
        .expect("RDKit-compatible argument construction succeeds before result-size use");
    assert!(
        getTopologicalTorsionGenerator(&oversized, None, true).is_err(),
        "row {row}: source undefined oversized shift must fail closed"
    );
    let generator =
        getTopologicalTorsionGenerator(&TopologicalTorsionArguments::default(), None, true)
            .expect("default generator");
    let mut arguments = FingerprintFuncArguments {
        custom_atom_invariants: Some(vec![1]),
        ..Default::default()
    };
    let custom_invariant_result = generator.getCountFingerprint(molecule, &mut arguments);
    if molecule.num_atoms() > 1 {
        assert!(
            custom_invariant_result.is_err(),
            "row {row}: genuinely short custom invariant vector must fail closed"
        );
    } else {
        assert!(
            custom_invariant_result.is_ok(),
            "row {row}: one invariant is valid for a one-atom molecule"
        );
    }
}

fn assert_active_profile_exactly(expected_profiles: &[&str]) {
    let corpus = corpus_smiles();
    let mut outcomes = golden_reader()
        .lines()
        .enumerate()
        .par_bridge()
        .map(|(row, content)| {
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                let content = content.unwrap_or_else(|error| {
                    panic!("failed to read golden line {}: {error}", row + 1)
                });
                let record: GoldenRecord = serde_json::from_str(&content).unwrap_or_else(|error| {
                    panic!("failed to parse golden line {}: {error}", row + 1)
                });
                let smiles = corpus
                    .get(row)
                    .unwrap_or_else(|| panic!("golden has extra row {row}"));
                assert_eq!(record.row, row, "row index changed");
                assert_eq!(record.smiles, *smiles, "row {row}: SMILES changed");
                if !record.rdkit_ok {
                    assert!(
                        record.error.is_some(),
                        "row {row}: missing RDKit parse error"
                    );
                    assert!(
                        record.profiles.is_empty(),
                        "row {row}: failed row has profiles"
                    );
                    assert!(
                        Molecule::from_smiles(smiles).is_err(),
                        "row {row}: parse mismatch"
                    );
                    return;
                }
                assert!(record.error.is_none(), "row {row}: success row has error");
                assert_eq!(
                    record
                        .profiles
                        .keys()
                        .map(String::as_str)
                        .collect::<Vec<_>>(),
                    expected_profiles,
                    "row {row}: profile set changed"
                );
                let molecule = Molecule::from_smiles(smiles)
                    .unwrap_or_else(|error| panic!("row {row} parse failed: {error}"));
                for (name, profile) in &record.profiles {
                    assert_modern_profile(row, name, &molecule, profile);
                }
                assert_legacy(row, &molecule, record.legacy.as_ref().expect("legacy"));
                assert_helpers(row, &molecule, record.helpers.as_ref().expect("helpers"));
                assert_error_profile(row, &molecule, record.errors.as_ref().expect("errors"));
            }))
            .err()
            .map(|payload| {
                let message = payload
                    .downcast_ref::<String>()
                    .cloned()
                    .or_else(|| {
                        payload
                            .downcast_ref::<&str>()
                            .map(|value| (*value).to_owned())
                    })
                    .unwrap_or_else(|| "non-string panic".to_owned());
                (row, message)
            })
            .map_or_else(|| (row, None), |failure| (row, Some(failure.1)))
        })
        .collect::<Vec<_>>();
    outcomes.sort_by_key(|(row, _)| *row);
    assert_eq!(
        outcomes.len(),
        corpus.len(),
        "golden must consume every corpus row"
    );
    assert!(
        outcomes
            .iter()
            .enumerate()
            .all(|(expected_row, (actual_row, _))| expected_row == *actual_row),
        "golden row sequence must be complete and contiguous"
    );
    let failures = outcomes
        .into_iter()
        .filter_map(|(row, failure)| failure.map(|message| (row, message)))
        .collect::<Vec<_>>();
    assert!(
        failures.is_empty(),
        "Topological Torsion exact focused parity failures: {failures:#?}"
    );
}

#[test]
fn topological_torsion_matches_every_generated_focused_profile_exactly() {
    match parity_data::profile_name().as_str() {
        "smiles_small" => assert_active_profile_exactly(&PROFILE_NAMES),
        "smiles_5000" => assert_active_profile_exactly(&CORPUS_PROFILE_NAMES),
        profile => panic!("unsupported Topological Torsion parity profile {profile}"),
    }
}

#[test]
#[ignore = "run explicitly with COSMOLKIT_PARITY_PROFILE=smiles_5000 in release mode"]
fn topological_torsion_matches_all_5000_rows_and_profiles_exactly() {
    assert_eq!(
        parity_data::profile_name(),
        "smiles_5000",
        "set COSMOLKIT_PARITY_PROFILE=smiles_5000"
    );
    assert_active_profile_exactly(&CORPUS_PROFILE_NAMES);
}
